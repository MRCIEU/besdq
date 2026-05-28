"""Unit tests for BESD query functionality."""

import unittest
import tempfile
from pathlib import Path
from unittest import mock
import io
import contextlib
import numpy as np
from besdq import BESDQueryEngine, BESDIndexBuilder, BESDQueryIndex
from besdq.besd_reader import IndexReader, calculate_p_value as reader_calculate_p_value
from besdq.sqlite_query import calculate_p_value as sqlite_calculate_p_value
from besdq.stats import _reconstruct_beta_se, _derive_af, _compute_n_from_data
from besdq import cli


# Path to test data
DATA_DIR = Path(__file__).parent.parent / "data"
WESTRA_BESD = str(DATA_DIR / "westra_eqtl_hg19")

# Known Westra dataset sample size
WESTRA_N = 1980


class TestStatsPureFunctions(unittest.TestCase):
    """Tests for _reconstruct_beta_se, _derive_af, _compute_n_from_data."""

    def test_reconstruct_scalar_n_round_trip(self):
        """Round-trip: reconstruct beta/SE from z-scores with scalar n."""
        n = 1980
        af = np.array([0.3, 0.5, 0.2])
        var_y = 1.0
        ses_orig = np.sqrt(var_y / (n * 2.0 * af * (1.0 - af)))
        betas_orig = np.array([-0.5, 0.3, 1.2])
        z = betas_orig / ses_orig

        betas_rec, ses_rec = _reconstruct_beta_se(z, af, n=n, var_y=var_y)

        np.testing.assert_allclose(ses_rec, ses_orig, rtol=1e-10)
        np.testing.assert_allclose(betas_rec, betas_orig, rtol=1e-10)

    def test_reconstruct_array_n_round_trip(self):
        """Round-trip with array n (VectorN-style)."""
        n = np.array([1000.0, 1500.0, 2000.0])
        af = np.array([0.4, 0.25, 0.45])
        ses_orig = np.sqrt(1.0 / (n * 2.0 * af * (1.0 - af)))
        betas_orig = np.array([0.1, -0.2, 0.05])
        z = betas_orig / ses_orig

        betas_rec, ses_rec = _reconstruct_beta_se(z, af, n=n)

        np.testing.assert_allclose(ses_rec, ses_orig, rtol=1e-10)
        np.testing.assert_allclose(betas_rec, betas_orig, rtol=1e-10)

    def test_reconstruct_nonunit_var_y(self):
        """SE scales by sqrt(var_y) relative to var_y=1 reconstruction."""
        n, af, var_y = 2000, 0.3, 2.5
        z = np.array([3.0, -5.0])
        _, ses_vy1 = _reconstruct_beta_se(z, np.array([af, af]), n=n, var_y=1.0)
        _, ses_vy = _reconstruct_beta_se(z, np.array([af, af]), n=n, var_y=var_y)

        np.testing.assert_allclose(ses_vy, ses_vy1 * np.sqrt(var_y), rtol=1e-10)

    def test_derive_af_round_trip(self):
        """_derive_af inverts _reconstruct_beta_se: SE recovered exactly."""
        n = 1980
        af_orig = np.array([0.15, 0.4, 0.48])
        var_y = 1.0
        ses = np.sqrt(var_y / (n * 2.0 * af_orig * (1.0 - af_orig)))

        af_derived = _derive_af(ses, n, var_y=var_y)
        np.testing.assert_allclose(af_derived, af_orig, rtol=1e-9)

        _, ses_rec = _reconstruct_beta_se(np.zeros(3), af_derived, n=n, var_y=var_y)
        np.testing.assert_allclose(ses_rec, ses, rtol=1e-9)

    def test_derive_af_zero_se_returns_nan(self):
        """se=0 should return NaN without raising."""
        result = _derive_af(np.array([0.0, 0.05]), n=1000)
        self.assertTrue(np.isnan(result[0]))
        self.assertTrue(np.isfinite(result[1]))

    def test_derive_af_negative_discriminant_returns_nan(self):
        """se too large (discriminant < 0) should return NaN without raising."""
        # se extremely large → 2*var_y/(n*se²) ≈ 0, inner ≈ 1 → ok
        # se very small → 2*var_y/(n*se²) >> 1 → inner < 0 → NaN
        result = _derive_af(np.array([1e-6]), n=10)
        self.assertTrue(np.isnan(result[0]))

    def test_reconstruct_float16_z_precision(self):
        """Float16 z-score quantisation gives < 0.1% relative error in beta."""
        n, af, var_y = 1980, 0.2, 1.0
        se_orig = float(np.sqrt(var_y / (n * 2.0 * af * (1.0 - af))))
        beta_orig = -0.5
        z_exact = beta_orig / se_orig
        z_f16 = float(np.float16(z_exact))

        betas_rec, _ = _reconstruct_beta_se(
            np.array([z_f16]), np.array([af]), n=n, var_y=var_y
        )
        rel_err = abs(float(betas_rec[0]) - beta_orig) / abs(beta_orig)
        self.assertLess(rel_err, 0.001)


class TestBESDQueryEngine(unittest.TestCase):
    """Test BESD query engine with westra_eqtl_hg19 dataset."""

    @classmethod
    def setUpClass(cls):
        cls.engine = BESDQueryEngine(WESTRA_BESD)

    def test_load_data(self):
        self.assertEqual(len(self.engine.snps), 506049)
        self.assertEqual(len(self.engine.probes), 5966)
        self.assertEqual(self.engine.besd.format_type, '3F')

    def test_single_position_query(self):
        associations = self.engine.query_cis_window(
            snp_chr='1', snp_start_kb=1191.87, snp_end_kb=1191.87,
            probe_chr='1', probe_start_kb=1140.818, probe_end_kb=1140.818,
        )
        self.assertEqual(len(associations), 1)
        self.assertEqual(associations[0]['snp_id'], 'rs3818646')
        self.assertEqual(associations[0]['probe_id'], 'ILMN_2349633')

    def test_range_query(self):
        associations = self.engine.query_cis_window(
            snp_chr='1', snp_start_kb=100, snp_end_kb=2000,
            probe_chr='1', probe_start_kb=1000, probe_end_kb=2000,
        )
        self.assertEqual(len(associations), 175)

    def test_p_value_calculation(self):
        associations = self.engine.query_cis_window(
            snp_chr='1', snp_start_kb=1191.87, snp_end_kb=1191.87,
            probe_chr='1', probe_start_kb=1140.818, probe_end_kb=1140.818,
        )
        self.assertEqual(len(associations), 1)
        assoc = associations[0]
        self.assertAlmostEqual(assoc['beta'], -0.436080, places=4)
        self.assertAlmostEqual(assoc['se'], 0.040022, places=4)
        self.assertGreaterEqual(assoc['pval'], 0)

    def test_beta_and_se_stored(self):
        associations = self.engine.query_cis_window(
            snp_chr='1', snp_start_kb=100, snp_end_kb=2000,
            probe_chr='1', probe_start_kb=1000, probe_end_kb=2000,
        )
        for assoc in associations:
            self.assertIn('beta', assoc)
            self.assertIn('se', assoc)
            self.assertIsInstance(assoc['beta'], float)
            self.assertIsInstance(assoc['se'], float)
            self.assertGreater(assoc['se'], 0)

    def test_snp_indexing(self):
        associations = self.engine.query_cis_window(
            snp_chr='1', snp_start_kb=1191.87, snp_end_kb=1191.87,
            probe_chr='1', probe_start_kb=1140.818, probe_end_kb=1140.818,
        )
        assoc = associations[0]
        self.assertEqual(assoc['snp_id'], 'rs3818646')
        self.assertEqual(assoc['snp_chr'], '1')
        self.assertEqual(assoc['snp_bp'], 1191870)
        self.assertEqual(assoc['a1'], 'T')
        self.assertEqual(assoc['a2'], 'C')

    def test_probe_indexing(self):
        associations = self.engine.query_cis_window(
            snp_chr='1', snp_start_kb=1191.87, snp_end_kb=1191.87,
            probe_chr='1', probe_start_kb=1140.818, probe_end_kb=1140.818,
        )
        assoc = associations[0]
        self.assertEqual(assoc['probe_id'], 'ILMN_2349633')
        self.assertEqual(assoc['probe_chr'], '1')
        self.assertEqual(assoc['probe_bp'], 1140818)
        self.assertEqual(assoc['gene'], 'TNFRSF18')

    def test_chromosome_filtering(self):
        associations = self.engine.query_cis_window(
            snp_chr='1', snp_start_kb=100, snp_end_kb=2000,
            probe_chr='1', probe_start_kb=1000, probe_end_kb=2000,
        )
        for assoc in associations:
            self.assertEqual(assoc['snp_chr'], '1')
            self.assertEqual(assoc['probe_chr'], '1')

    def test_position_filtering(self):
        associations = self.engine.query_cis_window(
            snp_chr='1', snp_start_kb=100, snp_end_kb=2000,
            probe_chr='1', probe_start_kb=1000, probe_end_kb=2000,
        )
        for assoc in associations:
            self.assertGreaterEqual(assoc['snp_bp'], 100000)
            self.assertLessEqual(assoc['snp_bp'], 2000000)
            self.assertGreaterEqual(assoc['probe_bp'], 1000000)
            self.assertLessEqual(assoc['probe_bp'], 2000000)

    def test_empty_query(self):
        associations = self.engine.query_cis_window(
            snp_chr='1', snp_start_kb=100, snp_end_kb=110,
            probe_chr='1', probe_start_kb=1900, probe_end_kb=2000,
        )
        self.assertEqual(len(associations), 0)

    def test_multi_snp_query(self):
        assocs1 = self.engine.query_by_snp_id('rs3818646')
        assocs2 = self.engine.query_by_snp_id('rs7515488')
        self.assertEqual(len(assocs1), 5)
        self.assertEqual(len(assocs2), 6)

    def test_multi_probe_query(self):
        assocs1 = self.engine.query_by_probe_id('ILMN_2349633')
        assocs2 = self.engine.query_by_probe_id('ILMN_2112256')
        self.assertEqual(len(assocs1), 20)
        self.assertEqual(len(assocs2), 20)
        probe_ids = {a['probe_id'] for a in assocs1 + assocs2}
        self.assertEqual(len(probe_ids), 2)

    def test_gene_query(self):
        associations = self.engine.query_by_gene('TNFRSF18')
        self.assertEqual(len(associations), 20)
        for assoc in associations:
            self.assertEqual(assoc['gene'], 'TNFRSF18')


class TestBESDQueryIndex(unittest.TestCase):
    """Test SQLite-indexed BESD query engine (z-score encoding)."""

    @classmethod
    def setUpClass(cls):
        cls.temp_db = tempfile.NamedTemporaryFile(suffix='.db', delete=False)
        cls.db_path = cls.temp_db.name
        cls.temp_db.close()

        # Westra has no AF in ESI, so ScalarN + explicit sample_size is required
        builder = BESDIndexBuilder(cls.db_path)
        builder.build(WESTRA_BESD, force=True, n_mode='scalar', sample_size=WESTRA_N)

        cls.index = BESDQueryIndex(cls.db_path)

    @classmethod
    def tearDownClass(cls):
        cls.index.close()
        Path(cls.db_path).unlink()

    def test_metadata_loading(self):
        self.assertEqual(self.index.metadata['format_type'], '3F')
        self.assertEqual(self.index.metadata['n_snps'], '506049')
        self.assertEqual(self.index.metadata['n_probes'], '5966')

    def test_single_position_query(self):
        associations = self.index.query_cis_window(
            snp_chr='1', snp_start_kb=1191.87, snp_end_kb=1191.87,
            probe_chr='1', probe_start_kb=1140.818, probe_end_kb=1140.818,
        )
        self.assertEqual(len(associations), 1)
        self.assertEqual(associations[0]['snp_id'], 'rs3818646')
        self.assertEqual(associations[0]['probe_id'], 'ILMN_2349633')

    def test_range_query(self):
        associations = self.index.query_cis_window(
            snp_chr='1', snp_start_kb=100, snp_end_kb=2000,
            probe_chr='1', probe_start_kb=1000, probe_end_kb=2000,
        )
        self.assertEqual(len(associations), 175)

    def test_p_value_calculation(self):
        """Beta/SE reconstructed from float16 z-score within 0.1% of original."""
        associations = self.index.query_cis_window(
            snp_chr='1', snp_start_kb=1191.87, snp_end_kb=1191.87,
            probe_chr='1', probe_start_kb=1140.818, probe_end_kb=1140.818,
        )
        self.assertEqual(len(associations), 1)
        assoc = associations[0]
        # Beta has ~0.1% float16 quantisation error; places=3 ≈ 0.1% tolerance for |beta|~0.4
        self.assertAlmostEqual(assoc['beta'], -0.436080, places=3)
        self.assertAlmostEqual(assoc['se'], 0.040022, places=4)
        self.assertGreaterEqual(assoc['pval'], 0)

    def test_schema_has_zscores_not_betas_ses(self):
        """Confirm probe_data stores zscores/n_scalar, not betas/ses."""
        import sqlite3
        conn = sqlite3.connect(self.db_path)
        cursor = conn.cursor()
        cursor.execute("PRAGMA table_info(probe_data)")
        cols = {row[1] for row in cursor.fetchall()}
        conn.close()
        self.assertIn('zscores', cols)
        self.assertIn('n_scalar', cols)
        self.assertIn('se_vector', cols)
        self.assertNotIn('n_vector', cols)
        self.assertNotIn('betas', cols)

    def test_epi_has_var_y_column(self):
        import sqlite3
        conn = sqlite3.connect(self.db_path)
        cursor = conn.cursor()
        cursor.execute("PRAGMA table_info(epi)")
        cols = {row[1] for row in cursor.fetchall()}
        conn.close()
        self.assertIn('var_y', cols)

    def test_query_by_probe_id(self):
        associations = self.index.query_by_probe_id('ILMN_2349633')
        self.assertEqual(len(associations), 20)
        self.assertTrue(all(a['probe_id'] == 'ILMN_2349633' for a in associations))

    def test_query_by_snp_id(self):
        associations = self.index.query_by_snp_id('rs3818646')
        self.assertEqual(len(associations), 5)
        self.assertTrue(all(a['snp_id'] == 'rs3818646' for a in associations))

    def test_consistency_with_besd_reader(self):
        """Index returns beta/SE within 0.1% of direct BESD reader."""
        engine = BESDQueryEngine(WESTRA_BESD)
        besd_results = engine.query_cis_window(
            snp_chr='1', snp_start_kb=1191.87, snp_end_kb=1191.87,
            probe_chr='1', probe_start_kb=1140.818, probe_end_kb=1140.818,
        )
        index_results = self.index.query_cis_window(
            snp_chr='1', snp_start_kb=1191.87, snp_end_kb=1191.87,
            probe_chr='1', probe_start_kb=1140.818, probe_end_kb=1140.818,
        )
        self.assertEqual(len(besd_results), len(index_results))
        self.assertEqual(besd_results[0]['snp_id'], index_results[0]['snp_id'])
        self.assertEqual(besd_results[0]['probe_id'], index_results[0]['probe_id'])
        # beta: allow up to 0.1% relative error (float16 z-score quantisation)
        beta_ref = besd_results[0]['beta']
        self.assertLess(abs(index_results[0]['beta'] - beta_ref) / abs(beta_ref), 1e-3)
        # se: reconstructed from derived AF, very close to original
        self.assertAlmostEqual(besd_results[0]['se'], index_results[0]['se'], places=5)

    def test_multi_snp_query_index(self):
        assocs1 = self.index.query_by_snp_id('rs3818646')
        assocs2 = self.index.query_by_snp_id('rs7515488')
        self.assertEqual(len(assocs1), 5)
        self.assertEqual(len(assocs2), 6)

    def test_multi_probe_query_index(self):
        assocs1 = self.index.query_by_probe_id('ILMN_2349633')
        assocs2 = self.index.query_by_probe_id('ILMN_2112256')
        self.assertEqual(len(assocs1), 20)
        self.assertEqual(len(assocs2), 20)
        probe_ids = {a['probe_id'] for a in assocs1 + assocs2}
        self.assertEqual(len(probe_ids), 2)

    def test_gene_query_index(self):
        associations = self.index.query_by_gene('TNFRSF18')
        self.assertEqual(len(associations), 20)
        for assoc in associations:
            self.assertEqual(assoc['gene'], 'TNFRSF18')

    def test_consistency_snp_query(self):
        engine = BESDQueryEngine(WESTRA_BESD)
        besd_results = engine.query_by_snp_id('rs3818646')
        index_results = self.index.query_by_snp_id('rs3818646')
        self.assertEqual(len(besd_results), len(index_results))
        besd_ids = {(a['snp_id'], a['probe_id']) for a in besd_results}
        index_ids = {(a['snp_id'], a['probe_id']) for a in index_results}
        self.assertEqual(besd_ids, index_ids)

    def test_consistency_probe_query(self):
        engine = BESDQueryEngine(WESTRA_BESD)
        besd_results = engine.query_by_probe_id('ILMN_2349633')
        index_results = self.index.query_by_probe_id('ILMN_2349633')
        self.assertEqual(len(besd_results), len(index_results))
        besd_ids = {(a['snp_id'], a['probe_id']) for a in besd_results}
        index_ids = {(a['snp_id'], a['probe_id']) for a in index_results}
        self.assertEqual(besd_ids, index_ids)

    def test_consistency_gene_query(self):
        engine = BESDQueryEngine(WESTRA_BESD)
        besd_results = engine.query_by_gene('TNFRSF18')
        index_results = self.index.query_by_gene('TNFRSF18')
        self.assertEqual(len(besd_results), len(index_results))
        besd_ids = {(a['snp_id'], a['probe_id']) for a in besd_results}
        index_ids = {(a['snp_id'], a['probe_id']) for a in index_results}
        self.assertEqual(besd_ids, index_ids)


class TestIndexReaderParsing(unittest.TestCase):
    """Test index file parsing edge cases."""

    def test_read_esi_row_idx_ignores_non_data_lines(self):
        with tempfile.NamedTemporaryFile(mode='w', suffix='.esi', delete=False) as f:
            f.write("# comment\n")
            f.write("\n")
            f.write("1 rsA 0.1 100 A G 0.4\n")
            f.write("bad line\n")
            f.write("1 rsB 0.2 200 C T 0.2\n")
            esi_path = f.name
        try:
            snps = IndexReader.read_esi(esi_path)
            self.assertEqual(len(snps), 2)
            self.assertEqual(snps[0]['row_idx'], 0)
            self.assertEqual(snps[1]['row_idx'], 1)
        finally:
            Path(esi_path).unlink()

    def test_read_epi_row_idx_ignores_non_data_lines(self):
        with tempfile.NamedTemporaryFile(mode='w', suffix='.epi', delete=False) as f:
            f.write("# comment\n")
            f.write("1 probeA 0.0 1000 GENE1 +\n")
            f.write("\n")
            f.write("bad line\n")
            f.write("1 probeB 0.0 2000 GENE2 -\n")
            epi_path = f.name
        try:
            probes = IndexReader.read_epi(epi_path)
            self.assertEqual(len(probes), 2)
            self.assertEqual(probes[0]['row_idx'], 0)
            self.assertEqual(probes[1]['row_idx'], 1)
        finally:
            Path(epi_path).unlink()


class TestPValueEdgeCases(unittest.TestCase):
    """Test numeric stability of p-value computations."""

    def test_p_value_extreme_z_is_bounded_reader(self):
        pval = reader_calculate_p_value(beta=1e9, se=1e-12)
        self.assertGreaterEqual(pval, 0.0)
        self.assertLessEqual(pval, 1.0)

    def test_p_value_extreme_z_is_bounded_sqlite(self):
        pval = sqlite_calculate_p_value(beta=1e9, se=1e-12)
        self.assertGreaterEqual(pval, 0.0)
        self.assertLessEqual(pval, 1.0)

    def test_p_value_zero_or_negative_se_defaults_to_one(self):
        self.assertEqual(reader_calculate_p_value(beta=1.0, se=0.0), 1.0)
        self.assertEqual(reader_calculate_p_value(beta=1.0, se=-1.0), 1.0)
        self.assertEqual(sqlite_calculate_p_value(beta=1.0, se=0.0), 1.0)
        self.assertEqual(sqlite_calculate_p_value(beta=1.0, se=-1.0), 1.0)


class TestCLIValidation(unittest.TestCase):
    """Test CLI parsing and argument validation edge cases."""

    def test_parse_chrpos_reversed_range_raises(self):
        with self.assertRaises(ValueError):
            cli.parse_chrpos("1:200-100")

    def test_parse_chrpos_malformed_position_raises(self):
        with self.assertRaises(ValueError):
            cli.parse_chrpos("1:abc")

    def test_cli_rejects_conflicting_identifier_queries(self):
        test_argv = [
            "besdq",
            "--beqtl-summary", "dummy",
            "--out", "out",
            "--snp", "rs1",
            "--probe", "probe1",
        ]
        with (
            mock.patch("sys.argv", test_argv),
            mock.patch("besdq.cli.BESDQueryEngine"),
            mock.patch("besdq.cli.BESDQueryIndex"),
            contextlib.redirect_stderr(io.StringIO()) as stderr
        ):
            with self.assertRaises(SystemExit) as cm:
                cli.main()
        self.assertEqual(cm.exception.code, 1)
        self.assertIn("mutually exclusive", stderr.getvalue())

    def test_cli_rejects_mixed_identifier_and_region_queries(self):
        test_argv = [
            "besdq",
            "--beqtl-summary", "dummy",
            "--out", "out",
            "--snp", "rs1",
            "--snp-chrpos", "1:100-200",
            "--probe-chrpos", "1:100-200",
        ]
        with (
            mock.patch("sys.argv", test_argv),
            mock.patch("besdq.cli.BESDQueryEngine"),
            mock.patch("besdq.cli.BESDQueryIndex"),
            contextlib.redirect_stderr(io.StringIO()) as stderr
        ):
            with self.assertRaises(SystemExit) as cm:
                cli.main()
        self.assertEqual(cm.exception.code, 1)
        self.assertIn("Cannot combine --snp/--probe/--gene", stderr.getvalue())


if __name__ == '__main__':
    unittest.main(verbosity=2)
