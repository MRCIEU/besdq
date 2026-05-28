"""Query module for SQLite-indexed BESD data."""

import sqlite3
import math
from pathlib import Path
from typing import List, Dict, Tuple
import numpy as np

from .stats import _reconstruct_beta_se


def norm_cdf(z: float) -> float:
    """Approximate normal CDF using Abramowitz and Stegun formula.

    Accurate to about 0.00012.
    """
    a1 = 0.254829592
    a2 = -0.284496736
    a3 = 1.421413741
    a4 = -1.453152027
    a5 = 1.061405429
    p = 0.3275911

    sign = 1 if z >= 0 else -1
    z = abs(z) / math.sqrt(2)

    t = 1.0 / (1.0 + p * z)
    y = 1.0 - (((((a5 * t + a4) * t) + a3) * t + a2) * t + a1) * t * math.exp(-z * z)

    return 0.5 * (1.0 + sign * y)


def calculate_p_value(beta: float, se: float) -> float:
    """Calculate two-tailed p-value with bounds checking."""
    if se <= 0 or not math.isfinite(beta) or not math.isfinite(se):
        return 1.0

    z_score = abs(beta / se)
    pval = 2.0 * (1.0 - norm_cdf(z_score))
    return min(1.0, max(0.0, pval))


class BESDQueryIndex:
    """Query BESD data from SQLite index database."""

    def __init__(self, db_path: str):
        self.db_path = Path(db_path)
        if not self.db_path.exists():
            raise FileNotFoundError(f"Database not found: {db_path}")

        self.conn = sqlite3.connect(str(self.db_path))
        self.conn.row_factory = sqlite3.Row
        self._load_metadata()

    def _load_metadata(self) -> None:
        cursor = self.conn.cursor()
        cursor.execute("SELECT key, value FROM besd_meta")
        self.metadata = {row['key']: row['value'] for row in cursor.fetchall()}

    def close(self) -> None:
        self.conn.close()

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()

    # ------------------------------------------------------------------
    # Range / metadata queries (unchanged)
    # ------------------------------------------------------------------

    def query_snp_range(self, chr_val: str, start_kb: float, end_kb: float) -> List[Dict]:
        start_bp = int(start_kb * 1000)
        end_bp = int(end_kb * 1000)

        cursor = self.conn.cursor()
        cursor.execute("""
            SELECT row_idx, chr, snp_id, genetic_dist, bp, a1, a2, freq
            FROM esi
            WHERE chr = ? AND bp >= ? AND bp <= ?
            ORDER BY bp
        """, (chr_val, start_bp, end_bp))

        return [dict(row) for row in cursor.fetchall()]

    def query_probe_range(self, chr_val: str, start_kb: float, end_kb: float) -> List[Dict]:
        start_bp = int(start_kb * 1000)
        end_bp = int(end_kb * 1000)

        cursor = self.conn.cursor()
        cursor.execute("""
            SELECT row_idx, trait_id, trait_name, trait_chr, trait_bp, trait_var, gene, context
            FROM epi
            WHERE trait_chr = ? AND trait_bp >= ? AND trait_bp <= ?
            ORDER BY trait_bp
        """, (chr_val, start_bp, end_bp))

        return [dict(row) for row in cursor.fetchall()]

    # ------------------------------------------------------------------
    # Core statistics reader — reconstructs beta/SE from z-scores
    # ------------------------------------------------------------------

    def get_probe_snps(self, probe_idx: int) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Return (snp_indices, betas, ses) for a probe, reconstructed from z-scores."""
        cursor = self.conn.cursor()
        cursor.execute("""
            SELECT pd.snp_indices, pd.zscores, pd.n_scalar, pd.se_vector, pd.snp_count,
                   e.trait_var
            FROM probe_data pd
            JOIN epi e ON e.row_idx = pd.probe_idx
            WHERE pd.probe_idx = ?
        """, (probe_idx,))

        row = cursor.fetchone()
        empty = (
            np.array([], dtype=np.int32),
            np.array([], dtype=np.float64),
            np.array([], dtype=np.float64),
        )
        if not row or row['snp_count'] == 0:
            return empty

        snp_indices = np.frombuffer(row['snp_indices'], dtype=np.int32)
        zscores = np.frombuffer(row['zscores'], dtype=np.float16).astype(np.float64)

        if row['se_vector'] is not None:
            # Vector mode: SE stored directly; beta = z * se, no AF needed
            ses = np.frombuffer(row['se_vector'], dtype=np.float16).astype(np.float64)
            betas = zscores * ses
            return snp_indices, betas, ses

        # Scalar mode: reconstruct SE from n + AF + trait_var
        if row['n_scalar'] is None:
            raise ValueError(f"No n_scalar or se_vector for probe_idx={probe_idx}")

        var_y = row['trait_var'] if row['trait_var'] is not None else 1.0
        n = float(row['n_scalar'])

        snp_list = [int(i) for i in snp_indices]
        placeholders = ','.join('?' * len(snp_list))
        cursor.execute(
            f"SELECT row_idx, freq FROM esi WHERE row_idx IN ({placeholders})",
            snp_list,
        )
        af_lookup = {r['row_idx']: r['freq'] for r in cursor.fetchall()}
        af_array = np.array(
            [af_lookup.get(i, np.nan) for i in snp_list],
            dtype=np.float64,
        )

        betas, ses = _reconstruct_beta_se(zscores, af_array, n, var_y=var_y)
        return snp_indices, betas, ses

    # ------------------------------------------------------------------
    # High-level query methods (signatures unchanged)
    # ------------------------------------------------------------------

    def query_cis_window(
        self,
        snp_chr: str, snp_start_kb: float, snp_end_kb: float,
        probe_chr: str, probe_start_kb: float, probe_end_kb: float,
    ) -> List[Dict]:
        snps = self.query_snp_range(snp_chr, snp_start_kb, snp_end_kb)
        probes = self.query_probe_range(probe_chr, probe_start_kb, probe_end_kb)

        snp_indices_set = {s['row_idx'] for s in snps}
        snp_by_idx = {s['row_idx']: s for s in snps}

        associations = []
        for probe in probes:
            probe_idx = probe['row_idx']
            snp_indices, betas, ses = self.get_probe_snps(probe_idx)

            for i, snp_idx in enumerate(snp_indices):
                if snp_idx in snp_indices_set:
                    snp = snp_by_idx[snp_idx]
                    beta = float(betas[i])
                    se = float(ses[i])
                    pval = calculate_p_value(beta, se)

                    associations.append({
                        'snp_id': snp['snp_id'],
                        'snp_chr': snp['chr'],
                        'snp_bp': snp['bp'],
                        'a1': snp['a1'],
                        'a2': snp['a2'],
                        'trait_id': probe['trait_id'],
                        'trait_chr': probe['trait_chr'],
                        'trait_bp': probe['trait_bp'],
                        'gene': probe['gene'],
                        'beta': beta,
                        'se': se,
                        'pval': pval,
                    })

        return associations

    def query_by_probe_id(self, probe_id: str) -> List[Dict]:
        cursor = self.conn.cursor()
        cursor.execute("""
            SELECT row_idx, trait_id, trait_chr, trait_bp, gene
            FROM epi WHERE trait_id = ?
        """, (probe_id,))
        probe_row = cursor.fetchone()
        if not probe_row:
            return []

        probe = dict(probe_row)
        snp_indices, betas, ses = self.get_probe_snps(probe['row_idx'])
        if len(snp_indices) == 0:
            return []

        snp_list = [int(i) for i in snp_indices]
        placeholders = ','.join('?' * len(snp_list))
        cursor.execute(
            f"SELECT row_idx, snp_id, chr, bp, a1, a2 FROM esi WHERE row_idx IN ({placeholders})",
            snp_list,
        )
        snp_by_idx = {row['row_idx']: dict(row) for row in cursor.fetchall()}

        associations = []
        for i, snp_idx in enumerate(snp_indices):
            snp_idx_int = int(snp_idx)
            if snp_idx_int in snp_by_idx:
                snp = snp_by_idx[snp_idx_int]
                beta = float(betas[i])
                se = float(ses[i])
                associations.append({
                    'snp_id': snp['snp_id'],
                    'snp_chr': snp['chr'],
                    'snp_bp': snp['bp'],
                    'a1': snp['a1'],
                    'a2': snp['a2'],
                    'trait_id': probe['trait_id'],
                    'trait_chr': probe['trait_chr'],
                    'trait_bp': probe['trait_bp'],
                    'gene': probe['gene'],
                    'beta': beta,
                    'se': se,
                    'pval': calculate_p_value(beta, se),
                })

        return associations

    def query_by_snp_id(self, snp_id: str) -> List[Dict]:
        cursor = self.conn.cursor()
        cursor.execute("SELECT row_idx, chr, snp_id, bp, a1, a2 FROM esi WHERE snp_id = ?", (snp_id,))
        snp_row = cursor.fetchone()
        if not snp_row:
            return []

        snp = dict(snp_row)
        target_snp_idx = snp['row_idx']

        associations = []
        cursor.execute("SELECT row_idx, trait_id, trait_chr, trait_bp, gene FROM epi")
        for probe_row in cursor.fetchall():
            probe_data = dict(probe_row)
            snp_indices, betas, ses = self.get_probe_snps(probe_data['row_idx'])

            match_indices = np.where(snp_indices == target_snp_idx)[0]
            if len(match_indices) > 0:
                match_idx = int(match_indices[0])
                beta = float(betas[match_idx])
                se = float(ses[match_idx])
                associations.append({
                    'snp_id': snp['snp_id'],
                    'snp_chr': snp['chr'],
                    'snp_bp': snp['bp'],
                    'a1': snp['a1'],
                    'a2': snp['a2'],
                    'trait_id': probe_data['trait_id'],
                    'trait_chr': probe_data['trait_chr'],
                    'trait_bp': probe_data['trait_bp'],
                    'gene': probe_data['gene'],
                    'beta': beta,
                    'se': se,
                    'pval': calculate_p_value(beta, se),
                })

        return associations

    def query_by_gene(self, gene_name: str) -> List[Dict]:
        cursor = self.conn.cursor()
        cursor.execute("""
            SELECT row_idx, trait_id, trait_chr, trait_bp, gene
            FROM epi WHERE gene = ?
        """, (gene_name,))
        probes = [dict(row) for row in cursor.fetchall()]
        if not probes:
            return []

        associations = []
        for probe in probes:
            snp_indices, betas, ses = self.get_probe_snps(probe['row_idx'])
            if len(snp_indices) == 0:
                continue

            snp_list = [int(i) for i in snp_indices]
            placeholders = ','.join('?' * len(snp_list))
            cursor.execute(
                f"SELECT row_idx, snp_id, chr, bp, a1, a2 FROM esi WHERE row_idx IN ({placeholders})",
                snp_list,
            )
            snp_by_idx = {row['row_idx']: dict(row) for row in cursor.fetchall()}

            for i, snp_idx in enumerate(snp_indices):
                snp_idx_int = int(snp_idx)
                if snp_idx_int in snp_by_idx:
                    snp = snp_by_idx[snp_idx_int]
                    beta = float(betas[i])
                    se = float(ses[i])
                    associations.append({
                        'snp_id': snp['snp_id'],
                        'snp_chr': snp['chr'],
                        'snp_bp': snp['bp'],
                        'a1': snp['a1'],
                        'a2': snp['a2'],
                        'trait_id': probe['trait_id'],
                        'trait_chr': probe['trait_chr'],
                        'trait_bp': probe['trait_bp'],
                        'gene': probe['gene'],
                        'beta': beta,
                        'se': se,
                        'pval': calculate_p_value(beta, se),
                    })

        return associations
