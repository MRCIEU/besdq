"""SQLite database builder for BESD data."""

import sqlite3
import sys
from pathlib import Path
from typing import List, Dict, Optional
import numpy as np

from .besd_reader import IndexReader, BESDReader
from .stats import _derive_af, _compute_n_from_data


class BESDIndexBuilder:
    """Build SQLite index from BESD files."""

    def __init__(self, db_path: str):
        self.db_path = Path(db_path)
        self.db_path.parent.mkdir(parents=True, exist_ok=True)

    def build(
        self,
        besd_prefix: str,
        force: bool = False,
        n_mode: str = 'scalar',
        sample_size: Optional[int] = None,
        trait_variance_path: Optional[str] = None,
    ) -> None:
        """Build index database from BESD files.

        Parameters
        ----------
        besd_prefix : path to BESD files without extension
        force : overwrite existing database
        n_mode : 'scalar' (one n per probe) or 'vector' (one n per SNP-probe pair)
        sample_size : explicit scalar sample size; if None for ScalarN+AF, computed from data
        trait_variance_path : optional two-column TSV (probe_id  var_y)
        """
        if n_mode not in ('scalar', 'vector'):
            print(f"Error: --n-mode must be 'scalar' or 'vector', got '{n_mode}'", file=sys.stderr)
            sys.exit(1)

        if self.db_path.exists() and not force:
            raise FileExistsError(
                f"Database {self.db_path} already exists. Use force=True to overwrite."
            )
        if self.db_path.exists():
            self.db_path.unlink()

        esi_path = f"{besd_prefix}.esi"
        epi_path = f"{besd_prefix}.epi"
        besd_path = f"{besd_prefix}.besd"

        print(f"Loading BESD files from {besd_prefix}...")
        snps = IndexReader.read_esi(esi_path)
        probes = IndexReader.read_epi(epi_path)
        besd = BESDReader(besd_path, len(probes))

        print(f"Loaded {len(snps)} SNPs and {len(probes)} probes")
        print(f"BESD format: SPARSE_FILE_TYPE_{besd.format_type}")

        has_af = any(s['freq'] is not None for s in snps)

        # Hard error: ScalarN without AF and without explicit sample size
        if n_mode == 'scalar' and not has_af and sample_size is None:
            print(
                "Error: --n-mode scalar with no allele frequency in ESI requires --sample-size.",
                file=sys.stderr,
            )
            sys.exit(1)

        # Load per-probe trait variance
        var_y_map = _load_var_y(trait_variance_path, {p['probe_id'] for p in probes})

        print(f"Creating database at {self.db_path}...")
        conn = sqlite3.connect(str(self.db_path))
        cursor = conn.cursor()

        self._create_schema(cursor)

        print("Writing metadata...")
        self._write_metadata(cursor, {
            'format_type': besd.format_type,
            'n_snps': str(len(snps)),
            'n_probes': str(len(probes)),
            'besd_path': besd_path,
            'esi_path': esi_path,
            'epi_path': epi_path,
        })

        print("Writing SNP index...")
        self._write_snps(cursor, snps)

        print("Writing probe index...")
        self._write_traits(cursor, probes, var_y_map)

        # AF lookup is only needed for ScalarN reconstruction.
        # VectorN stores SE directly so needs no AF at all.
        if n_mode == 'scalar':
            if has_af:
                af_lookup = {s['row_idx']: s['freq'] for s in snps if s['freq'] is not None}
            else:
                print(
                    "Warning: allele frequencies derived from beta/SE assuming var(y) = 1 "
                    "for all probes. Supply --trait-variance if phenotypes are not standardised.",
                    file=sys.stderr,
                )
                af_lookup = self._derive_snp_afs(besd, len(probes), sample_size, probes, var_y_map)
                for snp_idx, af in af_lookup.items():
                    cursor.execute("UPDATE esi SET freq = ? WHERE row_idx = ?", (af, snp_idx))
        else:
            af_lookup = {}

        print("Writing probe data...")
        self._write_probe_data(cursor, besd, len(probes), n_mode, sample_size, af_lookup, probes, var_y_map)

        print("Creating indices...")
        cursor.execute("CREATE INDEX idx_esi_chr_bp ON esi(chr, bp)")
        cursor.execute("CREATE INDEX idx_esi_snp_id ON esi(snp_id)")
        cursor.execute("CREATE INDEX idx_epi_trait_chr_bp ON epi(trait_chr, trait_bp)")
        cursor.execute("CREATE INDEX idx_epi_trait_id ON epi(trait_id)")

        conn.commit()
        conn.close()
        print(f"Database created successfully at {self.db_path}")

    # ------------------------------------------------------------------
    # Schema
    # ------------------------------------------------------------------

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
                snp_id TEXT NOT NULL,
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

        # Statistics stored as z-scores; n is either scalar (INTEGER) or
        # per-pair SE stored as float16 BLOB aligned to snp_indices (vector mode).
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

    # ------------------------------------------------------------------
    # Metadata / index writers
    # ------------------------------------------------------------------

    def _write_metadata(self, cursor: sqlite3.Cursor, metadata: Dict[str, str]) -> None:
        for key, value in metadata.items():
            cursor.execute(
                "INSERT INTO besd_meta (key, value) VALUES (?, ?)",
                (key, value),
            )

    def _write_snps(self, cursor: sqlite3.Cursor, snps: List[Dict]) -> None:
        for snp in snps:
            cursor.execute("""
                INSERT INTO esi (row_idx, chr, snp_id, genetic_dist, bp, a1, a2, freq)
                VALUES (?, ?, ?, ?, ?, ?, ?, ?)
            """, (
                snp['row_idx'], snp['chr'], snp['snp_id'], snp['genetic_dist'],
                snp['bp'], snp['a1'], snp['a2'], snp['freq'],
            ))

    def _write_traits(
        self,
        cursor: sqlite3.Cursor,
        probes: List[Dict],
        var_y_map: Dict[str, float],
    ) -> None:
        for probe in probes:
            trait_var = var_y_map.get(probe['probe_id'])  # None → NULL
            cursor.execute("""
                INSERT INTO epi (row_idx, trait_id, trait_name, trait_chr, trait_bp, trait_var, gene, context)
                VALUES (?, ?, ?, ?, ?, ?, ?, ?)
            """, (
                probe['row_idx'],
                probe['probe_id'],
                probe['probe_id'],  # BESD has no separate trait name; use probe_id
                probe.get('chr'),
                probe.get('probe_bp'),
                trait_var,
                probe.get('gene'),
                None,
            ))

    # ------------------------------------------------------------------
    # Probe data writer
    # ------------------------------------------------------------------

    def _write_probe_data(
        self,
        cursor: sqlite3.Cursor,
        besd: BESDReader,
        n_probes: int,
        n_mode: str,
        sample_size: Optional[int],
        af_lookup: Dict[int, float],
        probes: List[Dict],
        var_y_map: Dict[str, float],
    ) -> None:
        for probe_idx in range(n_probes):
            assocs = besd.get_probe_associations(probe_idx)
            probe_id = probes[probe_idx]['probe_id']
            probe_var_y = var_y_map.get(probe_id, 1.0)

            if not assocs:
                snp_indices = np.array([], dtype=np.int32)
                zscores_f16 = np.array([], dtype=np.float16)
                n_sc = sample_size
                n_vec_bytes = None
            else:
                raw_snp_idxs, betas_list, ses_list = zip(*assocs)
                snp_indices = np.array(raw_snp_idxs, dtype=np.int32)
                betas = np.array(betas_list, dtype=np.float64)
                ses = np.array(ses_list, dtype=np.float64)

                # Z-scores (float16 storage)
                with np.errstate(divide='ignore', invalid='ignore'):
                    z = np.where(ses > 0, betas / ses, 0.0)
                zscores_f16 = z.astype(np.float16)

                # AF array for this probe's SNPs
                af_array = np.array(
                    [af_lookup.get(int(i), np.nan) for i in snp_indices],
                    dtype=np.float64,
                )
                valid_af = np.isfinite(af_array) & (af_array > 0) & (af_array < 1)

                if n_mode == 'scalar':
                    if sample_size is not None:
                        n_sc = sample_size
                    else:
                        # Compute n from data for this probe
                        if valid_af.any():
                            n_per_pair = _compute_n_from_data(
                                ses[valid_af], af_array[valid_af], var_y=probe_var_y
                            )
                            finite_n = n_per_pair[np.isfinite(n_per_pair) & (n_per_pair > 0)]
                            if len(finite_n):
                                n_sc = int(round(float(np.mean(finite_n))))
                                # Warn if n varies substantially within this probe
                                if len(finite_n) > 1:
                                    cv = float(np.std(finite_n) / np.mean(finite_n))
                                    if cv > 0.1:
                                        print(
                                            f"Warning: probe {probe_id} has high n variance "
                                            f"(CV={cv:.2f}); consider --n-mode vector",
                                            file=sys.stderr,
                                        )
                            else:
                                n_sc = None
                        else:
                            n_sc = None
                    n_vec_bytes = None

                else:  # vector: store SE directly as float16
                    n_sc = None
                    n_vec_bytes = ses.astype(np.float16).tobytes()

            cursor.execute("""
                INSERT INTO probe_data (probe_idx, snp_count, snp_indices, zscores, n_scalar, se_vector)
                VALUES (?, ?, ?, ?, ?, ?)
            """, (
                probe_idx,
                len(assocs) if assocs else 0,
                snp_indices.tobytes(),
                zscores_f16.tobytes(),
                n_sc,
                n_vec_bytes,
            ))

            if (probe_idx + 1) % 1000 == 0:
                print(f"  Wrote {probe_idx + 1} / {n_probes} probes")

    # ------------------------------------------------------------------
    # AF derivation helper (first pass over BESD when ESI lacks freq)
    # ------------------------------------------------------------------

    def _derive_snp_afs(
        self,
        besd: BESDReader,
        n_probes: int,
        n: int,
        probes: List[Dict],
        var_y_map: Dict[str, float],
    ) -> Dict[int, float]:
        """First pass: derive AF per SNP from BESD data.

        AF is derived assuming var_y=1 for all probes (the reconstruction
        formula SE = sqrt(var_y/(n*2*AF*(1-AF))) is inverted with var_y=1).
        The var_y_map values are used only at query-time reconstruction via
        epi.var_y; they do not affect AF derivation here.
        """
        snp_first_se: Dict[int, float] = {}
        for probe_idx in range(n_probes):
            for snp_idx, _beta, se in besd.get_probe_associations(probe_idx):
                if snp_idx not in snp_first_se and se > 0:
                    snp_first_se[snp_idx] = se

        af_lookup: Dict[int, float] = {}
        for snp_idx, se in snp_first_se.items():
            af_val = float(_derive_af(se, n, var_y=1.0))
            if np.isfinite(af_val) and 0.0 < af_val < 1.0:
                af_lookup[snp_idx] = af_val

        return af_lookup


def _load_var_y(path: Optional[str], known_probe_ids: set) -> Dict[str, float]:
    """Load per-probe trait variance from a two-column TSV (probe_id  var_y)."""
    if path is None:
        return {}

    var_y_map: Dict[str, float] = {}
    with open(path, 'r') as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            parts = line.split()
            if len(parts) < 2:
                continue
            probe_id, val = parts[0], parts[1]
            if probe_id.lower() in ('probe_id', 'probe'):
                continue  # header row
            try:
                var_y_map[probe_id] = float(val)
            except ValueError:
                continue

    # Warn for unrecognised trait IDs
    unknown = set(var_y_map) - known_probe_ids
    for pid in sorted(unknown):
        print(
            f"Warning: --trait-variance file contains unrecognised trait ID '{pid}'",
            file=sys.stderr,
        )

    return var_y_map
