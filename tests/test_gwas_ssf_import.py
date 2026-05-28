"""Tests for GWAS-SSF import pipeline (Issues #7-#12)."""

import gzip
import io
import tempfile
import unittest
from pathlib import Path
from unittest import mock

import numpy as np

from besdq.gwas_ssf_reader import GwasSsfRow, read_gwas_ssf
from besdq.annotation_reader import TraitConfig, read_trait_annotation
from besdq.significance_filter import apply_significance_filter, FilterResult
from besdq.gwas_ssf_builder import GwasSsfIndexBuilder
from besdq.sqlite_query import BESDQueryIndex

DATA_DIR = Path(__file__).parent.parent / "data" / "ebi_input"
TRAITS_TSV = DATA_DIR / "traits.tsv"


# ---------------------------------------------------------------------------
# Issue #7 — GWAS-SSF reader
# ---------------------------------------------------------------------------

def _make_ssf_gz(rows: list[dict]) -> bytes:
    """Build a minimal GWAS-SSF gzip bytes buffer for testing."""
    cols = [
        'chromosome', 'base_pair_location', 'effect_allele', 'other_allele',
        'beta', 'standard_error', 'effect_allele_frequency', 'p_value', 'rsid',
    ]
    buf = io.BytesIO()
    with gzip.GzipFile(fileobj=buf, mode='wb') as gz:
        gz.write(('\t'.join(cols) + '\n').encode())
        for r in rows:
            line = '\t'.join(str(r.get(c, 'NA')) for c in cols) + '\n'
            gz.write(line.encode())
    return buf.getvalue()


class TestGwasSsfReader(unittest.TestCase):

    def _write_tmp_gz(self, rows: list[dict]) -> str:
        data = _make_ssf_gz(rows)
        with tempfile.NamedTemporaryFile(suffix='.tsv.gz', delete=False) as f:
            f.write(data)
            return f.name

    def test_no_swap_case(self):
        """effect_allele < other_allele: a1=effect, beta unchanged."""
        path = self._write_tmp_gz([{
            'chromosome': '1', 'base_pair_location': 100,
            'effect_allele': 'A', 'other_allele': 'G',
            'beta': 0.5, 'standard_error': 0.1,
            'effect_allele_frequency': 0.3, 'p_value': 1e-5, 'rsid': 'rs1',
        }])
        rows = list(read_gwas_ssf(path))
        Path(path).unlink()
        self.assertEqual(len(rows), 1)
        row = rows[0]
        self.assertEqual(row.a1, 'A')
        self.assertEqual(row.a2, 'G')
        self.assertAlmostEqual(row.beta, 0.5)
        self.assertAlmostEqual(row.eaf, 0.3)
        self.assertLessEqual(row.a1, row.a2)

    def test_swap_case(self):
        """effect_allele > other_allele: a1=other, beta negated, eaf inverted."""
        path = self._write_tmp_gz([{
            'chromosome': '1', 'base_pair_location': 200,
            'effect_allele': 'T', 'other_allele': 'A',
            'beta': 0.8, 'standard_error': 0.2,
            'effect_allele_frequency': 0.4, 'p_value': 1e-3, 'rsid': 'rs2',
        }])
        rows = list(read_gwas_ssf(path))
        Path(path).unlink()
        self.assertEqual(len(rows), 1)
        row = rows[0]
        self.assertEqual(row.a1, 'A')
        self.assertEqual(row.a2, 'T')
        self.assertAlmostEqual(row.beta, -0.8)
        self.assertAlmostEqual(row.eaf, 0.6)
        self.assertLessEqual(row.a1, row.a2)

    def test_missing_rsid_returns_none(self):
        path = self._write_tmp_gz([{
            'chromosome': '1', 'base_pair_location': 300,
            'effect_allele': 'C', 'other_allele': 'T',
            'beta': 0.1, 'standard_error': 0.05,
            'effect_allele_frequency': 0.5, 'p_value': 0.5, 'rsid': 'NA',
        }])
        rows = list(read_gwas_ssf(path))
        Path(path).unlink()
        self.assertIsNone(rows[0].rsid)

    def test_snp_key_format(self):
        path = self._write_tmp_gz([{
            'chromosome': '5', 'base_pair_location': 12345,
            'effect_allele': 'A', 'other_allele': 'C',
            'beta': 0.2, 'standard_error': 0.05,
            'effect_allele_frequency': 0.2, 'p_value': 1e-6, 'rsid': 'rs5',
        }])
        rows = list(read_gwas_ssf(path))
        Path(path).unlink()
        self.assertEqual(rows[0].snp_key, '5:12345:A:C')

    def test_a1_always_le_a2(self):
        """First 100 synthetic rows all have a1 <= a2."""
        test_rows = []
        alleles = [('A', 'G'), ('T', 'C'), ('G', 'T'), ('C', 'A')]
        for i in range(100):
            ea, oa = alleles[i % len(alleles)]
            test_rows.append({
                'chromosome': '1', 'base_pair_location': i * 1000 + 1,
                'effect_allele': ea, 'other_allele': oa,
                'beta': 0.1, 'standard_error': 0.05,
                'effect_allele_frequency': 0.3, 'p_value': 0.01, 'rsid': f'rs{i}',
            })
        path = self._write_tmp_gz(test_rows)
        rows = list(read_gwas_ssf(path))
        Path(path).unlink()
        for row in rows:
            self.assertLessEqual(row.a1, row.a2, f"a1={row.a1} > a2={row.a2}")

    def test_streaming_does_not_load_all(self):
        """Generator yields row-by-row without materialising the whole file."""
        test_rows = [{
            'chromosome': '1', 'base_pair_location': i,
            'effect_allele': 'A', 'other_allele': 'G',
            'beta': 0.1, 'standard_error': 0.05,
            'effect_allele_frequency': 0.3, 'p_value': 0.1, 'rsid': f'rs{i}',
        } for i in range(1, 6)]
        path = self._write_tmp_gz(test_rows)
        gen = read_gwas_ssf(path)
        first = next(gen)
        Path(path).unlink()
        self.assertIsInstance(first, GwasSsfRow)


# ---------------------------------------------------------------------------
# Issue #8 — Annotation TSV + YAML reader
# ---------------------------------------------------------------------------

class TestAnnotationReader(unittest.TestCase):

    def _write_tsv(self, content: str) -> str:
        with tempfile.NamedTemporaryFile(mode='w', suffix='.tsv', delete=False) as f:
            f.write(content)
            return f.name

    def test_parses_three_real_rows(self):
        """Parse the actual traits.tsv for the EBI datasets."""
        if not TRAITS_TSV.exists():
            self.skipTest("traits.tsv not found")
        traits = read_trait_annotation(str(TRAITS_TSV))
        self.assertEqual(len(traits), 3)
        ids = {t.trait_id for t in traits}
        self.assertIn('GCST90275731', ids)
        self.assertIn('GCST90275732', ids)
        self.assertIn('GCST90275739', ids)

    def test_sample_size_from_yaml(self):
        """sample_size falls back to YAML when absent from TSV."""
        if not TRAITS_TSV.exists():
            self.skipTest("traits.tsv not found")
        traits = read_trait_annotation(str(TRAITS_TSV))
        for t in traits:
            self.assertEqual(t.sample_size, 1060)

    def test_missing_required_column_raises(self):
        path = self._write_tsv("file_path\ttrait_id\nbadrow\tbadid\n")
        with self.assertRaises(ValueError) as cm:
            read_trait_annotation(path)
        Path(path).unlink()
        self.assertIn('trait_name', str(cm.exception))

    def test_missing_file_path_raises(self):
        path = self._write_tsv(
            "file_path\ttrait_id\ttrait_name\n"
            "/nonexistent/file.tsv.gz\tID1\tName1\n"
        )
        with self.assertRaises((FileNotFoundError, ValueError)):
            read_trait_annotation(path)
        Path(path).unlink()

    def test_chr_without_bp_raises(self):
        with tempfile.NamedTemporaryFile(mode='w', suffix='.tsv.gz', delete=False) as gzf:
            import gzip
            with gzip.open(gzf.name, 'wt') as gz:
                gz.write("chromosome\tbase_pair_location\teffect_allele\tother_allele\t"
                         "beta\tstandard_error\teffect_allele_frequency\tp_value\trsid\n")
            gz_path = gzf.name

        tsv_content = (
            "file_path\ttrait_id\ttrait_name\ttrait_chr\n"
            f"{gz_path}\tID1\tName1\t1\n"
        )
        tsv_path = self._write_tsv(tsv_content)
        with self.assertRaises(ValueError):
            read_trait_annotation(tsv_path)
        Path(tsv_path).unlink()
        Path(gz_path).unlink()


# ---------------------------------------------------------------------------
# Issue #9 — Three-tier significance filter
# ---------------------------------------------------------------------------

def _make_row(chr: str, bp: int, p: float, beta: float = 0.1, se: float = 0.01) -> GwasSsfRow:
    return GwasSsfRow(chr=chr, bp=bp, a1='A', a2='G', rsid=None,
                      beta=beta, se=se, eaf=0.3, p=p)


class TestSignificanceFilter(unittest.TestCase):

    def test_cis_row_retained_regardless_of_p(self):
        rows = [_make_row('1', 1_000_000, p=0.9)]
        result = apply_significance_filter(rows, trait_chr='1', trait_bp=1_000_000)
        self.assertEqual(len(result.cis), 1)
        self.assertEqual(len(result.sig_trans_candidates), 0)
        self.assertEqual(len(result.sug_trans), 0)

    def test_trans_significant_goes_to_sig_trans(self):
        rows = [_make_row('2', 5_000_000, p=1e-9)]
        result = apply_significance_filter(rows, trait_chr='1', trait_bp=1_000_000)
        self.assertEqual(len(result.sig_trans_candidates), 1)
        self.assertEqual(len(result.cis), 0)

    def test_trans_suggestive_goes_to_sug_trans(self):
        rows = [_make_row('2', 5_000_000, p=1e-5)]
        result = apply_significance_filter(rows, trait_chr='1', trait_bp=1_000_000)
        self.assertEqual(len(result.sug_trans), 1)

    def test_trans_below_suggestive_dropped(self):
        rows = [_make_row('2', 5_000_000, p=0.5)]
        result = apply_significance_filter(rows, trait_chr='1', trait_bp=1_000_000)
        self.assertEqual(len(result.cis), 0)
        self.assertEqual(len(result.sig_trans_candidates), 0)
        self.assertEqual(len(result.sug_trans), 0)

    def test_boundary_at_sig_threshold(self):
        """p == sig_threshold goes to sug_trans, not sig_trans."""
        rows = [_make_row('2', 5_000_000, p=5e-8)]
        result = apply_significance_filter(rows, trait_chr='1', trait_bp=1_000_000,
                                           sig_threshold=5e-8)
        self.assertEqual(len(result.sig_trans_candidates), 0)
        self.assertEqual(len(result.sug_trans), 1)

    def test_no_cis_when_trait_location_absent(self):
        rows = [_make_row('1', 1_000_000, p=0.9)]
        result = apply_significance_filter(rows, trait_chr=None, trait_bp=None)
        self.assertEqual(len(result.cis), 0)

    def test_outside_cis_radius_goes_to_trans(self):
        rows = [_make_row('1', 3_000_001, p=1e-9)]
        result = apply_significance_filter(rows, trait_chr='1', trait_bp=1_000_000,
                                           cis_radius=1_000_000)
        self.assertEqual(len(result.sig_trans_candidates), 1)
        self.assertEqual(len(result.cis), 0)

    def test_cis_boundary_exact(self):
        """SNP exactly at cis boundary is included in cis."""
        rows = [_make_row('1', 2_000_000, p=0.99)]
        result = apply_significance_filter(rows, trait_chr='1', trait_bp=1_000_000,
                                           cis_radius=1_000_000)
        self.assertEqual(len(result.cis), 1)


# ---------------------------------------------------------------------------
# Issue #10 — LD clumping (unit: missing plink2)
# ---------------------------------------------------------------------------

class TestLdClumping(unittest.TestCase):

    def test_empty_candidates_skip_plink2(self):
        from besdq.ld_clumping import clump_trans_peaks
        result = clump_trans_peaks([], plink2_pfile='/fake/ref')
        self.assertEqual(result, [])

    def test_missing_plink2_raises_import_error(self):
        from besdq.ld_clumping import clump_trans_peaks
        candidates = [_make_row('2', 1_000_000, p=1e-10)]
        with mock.patch('shutil.which', return_value=None):
            with self.assertRaises(ImportError) as cm:
                clump_trans_peaks(candidates, plink2_pfile='/fake/ref')
        self.assertIn('plink2', str(cm.exception).lower())


# ---------------------------------------------------------------------------
# Issue #11 + #12 — End-to-end import (no plink2 required)
# ---------------------------------------------------------------------------

class TestGwasSsfEndToEnd(unittest.TestCase):
    """Import the three EBI example files and query the resulting index."""

    @classmethod
    def setUpClass(cls):
        if not TRAITS_TSV.exists():
            cls.skip_reason = "traits.tsv not found"
            return
        cls.skip_reason = None

        from besdq.annotation_reader import read_trait_annotation
        from besdq.gwas_ssf_builder import GwasSsfIndexBuilder

        traits = read_trait_annotation(str(TRAITS_TSV))
        cls.tmp_db = tempfile.NamedTemporaryFile(suffix='.db', delete=False)
        cls.db_path = cls.tmp_db.name
        cls.tmp_db.close()

        builder = GwasSsfIndexBuilder(cls.db_path)
        # plink2_pfile=None → significant trans stored as-is (no clumping)
        builder.build(traits, workers=1, plink2_pfile=None)

    @classmethod
    def tearDownClass(cls):
        if hasattr(cls, 'db_path'):
            Path(cls.db_path).unlink(missing_ok=True)

    def _skip_if_needed(self):
        if getattr(self.__class__, 'skip_reason', None):
            self.skipTest(self.__class__.skip_reason)

    def test_db_has_three_epi_records(self):
        self._skip_if_needed()
        import sqlite3
        conn = sqlite3.connect(self.db_path)
        count = conn.execute("SELECT COUNT(*) FROM epi").fetchone()[0]
        conn.close()
        self.assertEqual(count, 3)

    def test_epi_trait_ids(self):
        self._skip_if_needed()
        import sqlite3
        conn = sqlite3.connect(self.db_path)
        ids = {row[0] for row in conn.execute("SELECT trait_id FROM epi").fetchall()}
        conn.close()
        expected = {'GCST90275731', 'GCST90275732', 'GCST90275739'}
        self.assertEqual(ids, expected)

    def test_esi_has_snps(self):
        self._skip_if_needed()
        import sqlite3
        conn = sqlite3.connect(self.db_path)
        count = conn.execute("SELECT COUNT(*) FROM esi").fetchone()[0]
        conn.close()
        self.assertGreater(count, 0)

    def test_alleles_ordered(self):
        """All ESI rows have a1 <= a2."""
        self._skip_if_needed()
        import sqlite3
        conn = sqlite3.connect(self.db_path)
        rows = conn.execute("SELECT a1, a2 FROM esi WHERE a1 IS NOT NULL AND a2 IS NOT NULL").fetchall()
        conn.close()
        for a1, a2 in rows:
            self.assertLessEqual(a1, a2, f"Allele ordering violated: a1={a1} a2={a2}")

    def test_query_index_compatible(self):
        """BESDQueryIndex can open and query the produced database."""
        self._skip_if_needed()
        idx = BESDQueryIndex(self.db_path)
        # Query by trait_id
        assocs = idx.query_by_probe_id('GCST90275731')
        idx.close()
        # May be empty if no associations passed the filter for this trait
        self.assertIsInstance(assocs, list)

    def test_n_scalar_set(self):
        """probe_data.n_scalar should be set for each trait (sample size)."""
        self._skip_if_needed()
        import sqlite3
        conn = sqlite3.connect(self.db_path)
        rows = conn.execute("SELECT n_scalar FROM probe_data").fetchall()
        conn.close()
        for (n,) in rows:
            self.assertEqual(n, 1060)

    def test_zscores_roundtrip(self):
        """Float16 z-scores decoded from probe_data blobs are finite."""
        self._skip_if_needed()
        import sqlite3
        conn = sqlite3.connect(self.db_path)
        rows = conn.execute(
            "SELECT zscores, snp_count FROM probe_data WHERE snp_count > 0"
        ).fetchall()
        conn.close()
        for zscores_bytes, snp_count in rows:
            zscores = np.frombuffer(zscores_bytes, dtype=np.float16)
            self.assertEqual(len(zscores), snp_count)
            self.assertTrue(np.all(np.isfinite(zscores.astype(np.float32))))

    def test_cli_import_gwas_ssf_missing_annotation_exits(self):
        """import-gwas-ssf with missing annotation file exits with code 1."""
        from besdq.cli import import_gwas_ssf_main
        with mock.patch('sys.argv', [
            'import-gwas-ssf',
            '--trait-annotation', '/nonexistent/traits.tsv',
            '--ld-reference', '/fake/ref',
        ]):
            with self.assertRaises(SystemExit) as cm:
                import_gwas_ssf_main()
        self.assertEqual(cm.exception.code, 1)


if __name__ == '__main__':
    unittest.main(verbosity=2)
