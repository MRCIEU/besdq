"""Two-pass GWAS-SSF builder: parallel Pass 1 (filtering) + serial Pass 2 (index)."""

import json
import sqlite3
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np

from .annotation_reader import TraitConfig
from .gwas_ssf_reader import GwasSsfRow, read_gwas_ssf
from .significance_filter import apply_significance_filter


@dataclass
class _TraitResult:
    """Intermediate result from Pass 1 for one trait."""
    trait: TraitConfig
    rows: List[GwasSsfRow] = field(default_factory=list)
    n_total_read: int = 0
    n_retained: int = 0


def _pass1_worker(args: tuple) -> _TraitResult:
    """Pass 1: stream and filter one GWAS-SSF file."""
    (
        file_path, trait_id, trait_name, trait_chr, trait_bp,
        sample_size, trait_var, gene, context, study_metadata,
        cis_radius, sig_threshold, sug_threshold, plink2_pfile,
        sig_radius, clump_r2, clump_kb,
    ) = args

    trait = TraitConfig(
        file_path=file_path,
        trait_id=trait_id,
        trait_name=trait_name,
        trait_chr=trait_chr,
        trait_bp=trait_bp,
        sample_size=sample_size,
        trait_var=trait_var,
        gene=gene,
        context=context,
        study_metadata=study_metadata,
    )

    result = _TraitResult(trait=trait)
    all_rows = list(read_gwas_ssf(file_path))
    result.n_total_read = len(all_rows)

    filter_result = apply_significance_filter(
        all_rows,
        trait_chr=trait_chr,
        trait_bp=trait_bp,
        cis_radius=cis_radius,
        sig_threshold=sig_threshold,
        sug_threshold=sug_threshold,
    )

    retained = list(filter_result.cis) + list(filter_result.sug_trans)

    # LD clumping for significant trans candidates
    if filter_result.sig_trans_candidates and plink2_pfile:
        try:
            from .ld_clumping import clump_trans_peaks
            trans_retained = clump_trans_peaks(
                filter_result.sig_trans_candidates,
                plink2_pfile=plink2_pfile,
                sig_radius=sig_radius,
                clump_r2=clump_r2,
                clump_kb=clump_kb,
                all_rows=all_rows,
            )
            retained.extend(trans_retained)
        except ImportError:
            # plink2 not available; store all sig candidates as-is
            retained.extend(filter_result.sig_trans_candidates)
    elif filter_result.sig_trans_candidates:
        # No LD reference; store all sig candidates as-is
        retained.extend(filter_result.sig_trans_candidates)

    # Deduplicate by snp_key
    seen: set = set()
    deduped = []
    for row in retained:
        key = row.snp_key
        if key not in seen:
            seen.add(key)
            deduped.append(row)

    result.rows = deduped
    result.n_retained = len(deduped)
    return result


def _snp_sort_key(snp_key: str) -> Tuple:
    """Sort SNP keys by chromosome (numeric then alpha) then position."""
    parts = snp_key.split(':')
    chr_str = parts[0]
    bp = int(parts[1]) if len(parts) > 1 else 0
    try:
        return (0, int(chr_str), bp)
    except ValueError:
        return (1, chr_str, bp)


class GwasSsfIndexBuilder:
    """Build a BESD-compatible SQLite index from GWAS-SSF files."""

    def __init__(self, db_path: str):
        self.db_path = Path(db_path)
        self.db_path.parent.mkdir(parents=True, exist_ok=True)

    def build(
        self,
        traits: List[TraitConfig],
        workers: int = 1,
        cis_radius: int = 1_000_000,
        sig_threshold: float = 5e-8,
        sug_threshold: float = 1e-4,
        plink2_pfile: Optional[str] = None,
        sig_radius: int = 500_000,
        clump_r2: float = 0.01,
        clump_kb: int = 10_000,
    ) -> None:
        """Build the index database from annotation-driven trait configs."""
        if self.db_path.exists():
            self.db_path.unlink()

        # ------- Pass 1: parallel per-file filtering -------
        _log("Pass 1: filtering traits…")
        worker_args = [
            (
                t.file_path, t.trait_id, t.trait_name, t.trait_chr, t.trait_bp,
                t.sample_size, t.trait_var, t.gene, t.context, t.study_metadata,
                cis_radius, sig_threshold, sug_threshold, plink2_pfile,
                sig_radius, clump_r2, clump_kb,
            )
            for t in traits
        ]

        trait_results: Dict[str, _TraitResult] = {}

        if workers == 1:
            for args in worker_args:
                r = _pass1_worker(args)
                trait_results[r.trait.trait_id] = r
                _log(f"  {r.trait.trait_id}: {r.n_retained}/{r.n_total_read} rows retained")
        else:
            with ProcessPoolExecutor(max_workers=workers) as pool:
                futures = {pool.submit(_pass1_worker, args): args[1] for args in worker_args}
                for fut in as_completed(futures):
                    r = fut.result()
                    trait_results[r.trait.trait_id] = r
                    _log(f"  {r.trait.trait_id}: {r.n_retained}/{r.n_total_read} rows retained")

        # ------- Pass 2: serial index construction -------
        _log("Pass 2: building index…")

        # Consolidate SNP universe across all traits, sorted stably
        all_snp_keys: set = set()
        for r in trait_results.values():
            for row in r.rows:
                all_snp_keys.add(row.snp_key)

        sorted_keys = sorted(all_snp_keys, key=_snp_sort_key)
        snp_key_to_idx: Dict[str, int] = {k: i for i, k in enumerate(sorted_keys)}
        _log(f"  ESI size: {len(sorted_keys)} unique SNPs")

        conn = sqlite3.connect(str(self.db_path))
        cursor = conn.cursor()
        self._create_schema(cursor)

        # Write ESI
        snp_data: Dict[str, dict] = {}
        for r in trait_results.values():
            for row in r.rows:
                key = row.snp_key
                if key not in snp_data:
                    snp_data[key] = {
                        'chr': row.chr,
                        'bp': row.bp,
                        'a1': row.a1,
                        'a2': row.a2,
                        'rsid': row.rsid,
                        'eaf': row.eaf,
                    }

        for idx, key in enumerate(sorted_keys):
            d = snp_data[key]
            cursor.execute(
                "INSERT INTO esi (row_idx, chr, snp_id, genetic_dist, bp, a1, a2, freq) "
                "VALUES (?, ?, ?, ?, ?, ?, ?, ?)",
                (idx, d['chr'], d['rsid'], None, d['bp'], d['a1'], d['a2'], d['eaf']),
            )

        _log("  ESI written")

        # Write EPI and probe_data
        study_meta_written = False
        for epi_idx, trait_cfg in enumerate(traits):
            tid = trait_cfg.trait_id
            r = trait_results.get(tid)
            rows = r.rows if r else []

            cursor.execute(
                "INSERT INTO epi (row_idx, trait_id, trait_name, trait_chr, trait_bp, "
                "trait_var, gene, context) VALUES (?, ?, ?, ?, ?, ?, ?, ?)",
                (
                    epi_idx, trait_cfg.trait_id, trait_cfg.trait_name,
                    trait_cfg.trait_chr, trait_cfg.trait_bp,
                    trait_cfg.trait_var if trait_cfg.trait_var != 1.0 else None,
                    trait_cfg.gene, trait_cfg.context,
                ),
            )

            if rows:
                snp_indices = np.array(
                    [snp_key_to_idx[row.snp_key] for row in rows], dtype=np.int32
                )
                zscores = np.array(
                    [row.beta / row.se if row.se > 0 else 0.0 for row in rows],
                    dtype=np.float64,
                ).astype(np.float16)
            else:
                snp_indices = np.array([], dtype=np.int32)
                zscores = np.array([], dtype=np.float16)

            cursor.execute(
                "INSERT INTO probe_data (probe_idx, snp_count, snp_indices, zscores, "
                "n_scalar, se_vector) VALUES (?, ?, ?, ?, ?, ?)",
                (
                    epi_idx, len(rows),
                    snp_indices.tobytes(), zscores.tobytes(),
                    trait_cfg.sample_size, None,
                ),
            )

            # Store study metadata once (first trait with metadata)
            if not study_meta_written and trait_cfg.study_metadata:
                cursor.execute(
                    "INSERT INTO besd_meta (key, value) VALUES (?, ?)",
                    ('study_metadata', json.dumps(trait_cfg.study_metadata)),
                )
                study_meta_written = True

        # Store basic metadata
        cursor.execute(
            "INSERT INTO besd_meta (key, value) VALUES (?, ?)",
            ('n_snps', str(len(sorted_keys))),
        )
        cursor.execute(
            "INSERT INTO besd_meta (key, value) VALUES (?, ?)",
            ('n_traits', str(len(traits))),
        )
        cursor.execute(
            "INSERT INTO besd_meta (key, value) VALUES (?, ?)",
            ('source', 'gwas-ssf'),
        )

        # Indexes
        cursor.execute("CREATE INDEX idx_esi_chr_bp ON esi(chr, bp)")
        cursor.execute("CREATE INDEX idx_esi_snp_id ON esi(snp_id)")
        cursor.execute("CREATE INDEX idx_epi_trait_chr_bp ON epi(trait_chr, trait_bp)")
        cursor.execute("CREATE INDEX idx_epi_trait_id ON epi(trait_id)")

        conn.commit()
        conn.close()
        _log(f"Pass 2 complete. Database written to {self.db_path}")

    def _create_schema(self, cursor: sqlite3.Cursor) -> None:
        cursor.execute("""
            CREATE TABLE besd_meta (
                key TEXT PRIMARY KEY,
                value TEXT NOT NULL
            )
        """)
        cursor.execute("""
            CREATE TABLE esi (
                row_idx INTEGER PRIMARY KEY,
                chr TEXT NOT NULL,
                snp_id TEXT,
                genetic_dist REAL,
                bp INTEGER NOT NULL,
                a1 TEXT,
                a2 TEXT,
                freq REAL
            )
        """)
        cursor.execute("""
            CREATE TABLE epi (
                row_idx INTEGER PRIMARY KEY,
                trait_id TEXT NOT NULL,
                trait_name TEXT NOT NULL,
                trait_chr TEXT,
                trait_bp INTEGER,
                trait_var REAL,
                gene TEXT,
                context TEXT
            )
        """)
        cursor.execute("""
            CREATE TABLE probe_data (
                probe_idx INTEGER PRIMARY KEY,
                snp_count INTEGER NOT NULL,
                snp_indices BLOB NOT NULL,
                zscores BLOB NOT NULL,
                n_scalar INTEGER,
                se_vector BLOB
            )
        """)


def _log(msg: str) -> None:
    import datetime
    ts = datetime.datetime.now().strftime('%H:%M:%S')
    print(f"[{ts}] {msg}", file=sys.stderr)
