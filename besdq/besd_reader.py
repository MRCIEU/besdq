"""Core BESD reader functionality."""

import math
import struct
from pathlib import Path
from typing import List, Tuple, Dict


def norm_cdf(z: float) -> float:
    """Approximate normal CDF using Abramowitz and Stegun formula.

    Accurate to about 0.00012.
    """
    # Constants
    a1 =  0.254829592
    a2 = -0.284496736
    a3 =  1.421413741
    a4 = -1.453152027
    a5 =  1.061405429
    p  =  0.3275911

    # Save the sign of z
    sign = 1 if z >= 0 else -1
    z = abs(z) / math.sqrt(2)

    # Calculate t
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


class IndexReader:
    """Read .esi or .epi index files."""

    @staticmethod
    def read_esi(esi_path: str) -> List[Dict]:
        """Read SNP index file."""
        snps = []
        with open(esi_path, 'r') as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith('#'):
                    continue

                parts = line.split()
                if len(parts) < 4:
                    continue

                try:
                    snps.append({
                        'row_idx': len(snps),
                        'chr': parts[0],
                        'snp_id': parts[1],
                        'genetic_dist': float(parts[2]) if parts[2] != 'NA' else None,
                        'bp': int(parts[3]),
                        'a1': parts[4] if len(parts) > 4 else None,
                        'a2': parts[5] if len(parts) > 5 else None,
                        'freq': float(parts[6]) if len(parts) > 6 and parts[6] != 'NA' else None,
                    })
                except (ValueError, IndexError):
                    continue

        return snps

    @staticmethod
    def read_epi(epi_path: str) -> List[Dict]:
        """Read probe index file."""
        probes = []
        with open(epi_path, 'r') as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith('#'):
                    continue

                parts = line.split()
                if len(parts) < 4:
                    continue

                try:
                    probes.append({
                        'row_idx': len(probes),
                        'chr': parts[0],
                        'probe_id': parts[1],
                        'genetic_dist': float(parts[2]) if len(parts) > 2 and parts[2] != 'NA' else None,
                        'probe_bp': int(parts[3]),
                        'gene': parts[4] if len(parts) > 4 else None,
                        'orientation': parts[5] if len(parts) > 5 else None,
                    })
                except (ValueError, IndexError):
                    continue

        return probes


class BESDReader:
    """Read SPARSE_FILE_TYPE_3 and SPARSE_FILE_TYPE_3F BESD format files."""

    MAGIC_SPARSE_3F = 0x40400000  # Sparse format 3F
    MAGIC_SPARSE_3 = 3             # Sparse format 3
    RESERVEDUNITS = 16

    def __init__(self, besd_path: str, n_probes: int):
        """Initialize BESD reader."""
        self.besd_path = besd_path
        self.n_probes = n_probes
        self.format_type = None
        self._cols = None
        self._rowid = None
        self._val = None
        self._val_num = None
        self._open()

    def _open(self):
        """Open and parse BESD file."""
        with open(self.besd_path, 'rb') as f:
            # Read magic
            magic_bytes = f.read(4)
            magic = struct.unpack('<I', magic_bytes)[0]

            if magic == self.MAGIC_SPARSE_3F:
                self.format_type = '3F'
                self._parse_format_3f(f)
            elif magic == self.MAGIC_SPARSE_3:
                self.format_type = '3'
                self._parse_format_3(f)
            else:
                raise ValueError(f"Unsupported BESD magic: 0x{magic:08x}")

    def _parse_format_3f(self, f):
        """Parse SPARSE_FILE_TYPE_3F format."""
        # Read valNum (number of associations)
        val_num_bytes = f.read(8)
        self._val_num = struct.unpack('<Q', val_num_bytes)[0]

        # Calculate structure sizes
        col_num = (self.n_probes << 1) + 1  # probNum * 2 + 1

        # Read columns (offsets)
        col_bytes = f.read(col_num * 8)
        self._cols = list(struct.unpack(f'<{col_num}Q', col_bytes))

        # Read row IDs (SNP indices)
        row_id_bytes = f.read(self._val_num * 4)
        self._rowid = list(struct.unpack(f'<{self._val_num}I', row_id_bytes))

        # Read values (betas and SEs)
        val_bytes = f.read(self._val_num * 4)
        self._val = list(struct.unpack(f'<{self._val_num}f', val_bytes))

    def _parse_format_3(self, f):
        """Parse SPARSE_FILE_TYPE_3 format."""
        # Read reserved units (sample size, n_snps, n_probes, etc.)
        reserved_bytes = f.read((self.RESERVEDUNITS - 1) * 4)
        struct.unpack(f'<{self.RESERVEDUNITS - 1}I', reserved_bytes)

        # Read valNum (number of associations)
        val_num_bytes = f.read(8)
        self._val_num = struct.unpack('<Q', val_num_bytes)[0]

        # Calculate structure sizes
        col_num = (self.n_probes << 1) + 1  # probNum * 2 + 1

        # Read columns (offsets)
        col_bytes = f.read(col_num * 8)
        self._cols = list(struct.unpack(f'<{col_num}Q', col_bytes))

        # Read row IDs (SNP indices)
        row_id_bytes = f.read(self._val_num * 4)
        self._rowid = list(struct.unpack(f'<{self._val_num}I', row_id_bytes))

        # Read values (betas and SEs)
        val_bytes = f.read(self._val_num * 4)
        self._val = list(struct.unpack(f'<{self._val_num}f', val_bytes))

    def get_probe_associations(self, probe_idx: int) -> List[Tuple[int, float, float]]:
        """Get all associations for a probe.

        Returns: list of (snp_idx, beta, se) tuples

        Structure:
        - _cols[2p] = start of betas for probe p
        - _cols[2p+1] = start of SEs for probe p (= end of betas)
        - _cols[2p+2] = end of SEs for probe p
        """
        if probe_idx >= self.n_probes:
            return []

        if self._cols is None or self._rowid is None or self._val is None:
            return []

        # Get column indices for this probe
        beta_start = self._cols[probe_idx << 1]          # Start of betas
        se_start = self._cols[(probe_idx << 1) + 1]      # Start of SEs (= end of betas)

        num_assocs = se_start - beta_start

        if num_assocs <= 0:
            return []

        # Extract associations
        associations = []
        for i in range(int(num_assocs)):
            beta_idx = int(beta_start + i)
            se_idx = int(se_start + i)

            if beta_idx < len(self._rowid) and beta_idx < len(self._val) and se_idx < len(self._val):
                snp_idx = self._rowid[beta_idx]
                beta = self._val[beta_idx]
                se = self._val[se_idx]
                associations.append((snp_idx, beta, se))

        return associations


class BESDQueryEngine:
    """Query BESD data directly without database."""

    def __init__(self, besd_prefix: str):
        """Initialize from BESD file prefix (without extension)."""
        self.prefix = Path(besd_prefix)

        # Load indices
        esi_path = f"{besd_prefix}.esi"
        epi_path = f"{besd_prefix}.epi"
        besd_path = f"{besd_prefix}.besd"

        self.snps = IndexReader.read_esi(esi_path)
        self.probes = IndexReader.read_epi(epi_path)
        self.besd = BESDReader(besd_path, len(self.probes))

        # Build lookup maps
        self.snp_by_idx = {s['row_idx']: s for s in self.snps}
        self.probe_by_idx = {p['row_idx']: p for p in self.probes}

    def query_snp_range(self, chr_val: str, start_kb: float, end_kb: float) -> List[Dict]:
        """Query SNPs in chromosome range (coordinates in kb)."""
        start_bp = int(start_kb * 1000)
        end_bp = int(end_kb * 1000)

        result = []
        for snp in self.snps:
            if snp['chr'] == chr_val and start_bp <= snp['bp'] <= end_bp:
                result.append(snp)
        return sorted(result, key=lambda x: x['bp'])

    def query_probe_range(self, chr_val: str, start_kb: float, end_kb: float) -> List[Dict]:
        """Query probes in chromosome range (coordinates in kb)."""
        start_bp = int(start_kb * 1000)
        end_bp = int(end_kb * 1000)

        result = []
        for probe in self.probes:
            if probe['chr'] == chr_val and start_bp <= probe['probe_bp'] <= end_bp:
                result.append(probe)
        return sorted(result, key=lambda x: x['probe_bp'])

    def query_cis_window(
        self,
        snp_chr: str, snp_start_kb: float, snp_end_kb: float,
        probe_chr: str, probe_start_kb: float, probe_end_kb: float,
    ) -> List[Dict]:
        """Query cis-window: return all SNP-probe associations in both regions."""

        snps = self.query_snp_range(snp_chr, snp_start_kb, snp_end_kb)
        probes = self.query_probe_range(probe_chr, probe_start_kb, probe_end_kb)

        # Build SNP index set for fast lookup
        snp_indices_set = {s['row_idx'] for s in snps}

        # Find all associations
        associations = []
        for probe in probes:
            probe_idx = probe['row_idx']
            assocs = self.besd.get_probe_associations(probe_idx)

            for snp_idx, beta, se in assocs:
                if snp_idx in snp_indices_set:
                    snp = self.snp_by_idx[snp_idx]
                    pval = calculate_p_value(beta, se)

                    associations.append({
                        'snp_id': snp['snp_id'],
                        'snp_chr': snp['chr'],
                        'snp_bp': snp['bp'],
                        'a1': snp['a1'],
                        'a2': snp['a2'],
                        'trait_id': probe['probe_id'],
                        'trait_chr': probe['chr'],
                        'trait_bp': probe['probe_bp'],
                        'gene': probe['gene'],
                        'beta': beta,
                        'se': se,
                        'pval': pval,
                    })

        return associations

    def query_by_snp_id(self, snp_id: str) -> List[Dict]:
        """Query all associations for a specific SNP.

        Args:
            snp_id: SNP ID to query

        Returns:
            List of associations for the SNP
        """
        # Find SNP
        snp_list = [s for s in self.snps if s['snp_id'] == snp_id]
        if not snp_list:
            return []

        snp = snp_list[0]
        target_snp_idx = snp['row_idx']

        # Search all probes for associations with this SNP
        associations = []
        for probe in self.probes:
            probe_idx = probe['row_idx']
            assocs = self.besd.get_probe_associations(probe_idx)

            for snp_idx, beta, se in assocs:
                if snp_idx == target_snp_idx:
                    pval = calculate_p_value(beta, se)

                    associations.append({
                        'snp_id': snp['snp_id'],
                        'snp_chr': snp['chr'],
                        'snp_bp': snp['bp'],
                        'a1': snp['a1'],
                        'a2': snp['a2'],
                        'trait_id': probe['probe_id'],
                        'trait_chr': probe['chr'],
                        'trait_bp': probe['probe_bp'],
                        'gene': probe['gene'],
                        'beta': beta,
                        'se': se,
                        'pval': pval,
                    })

        return associations

    def query_by_probe_id(self, probe_id: str) -> List[Dict]:
        """Query all associations for a specific probe.

        Args:
            probe_id: Probe ID to query

        Returns:
            List of associations for the probe
        """
        # Find probe
        probe_list = [p for p in self.probes if p['probe_id'] == probe_id]
        if not probe_list:
            return []

        probe = probe_list[0]
        probe_idx = probe['row_idx']
        assocs = self.besd.get_probe_associations(probe_idx)

        associations = []
        for snp_idx, beta, se in assocs:
            snp = self.snp_by_idx[snp_idx]
            pval = calculate_p_value(beta, se)

            associations.append({
                'snp_id': snp['snp_id'],
                'snp_chr': snp['chr'],
                'snp_bp': snp['bp'],
                'a1': snp['a1'],
                'a2': snp['a2'],
                'trait_id': probe['probe_id'],
                'trait_chr': probe['chr'],
                'trait_bp': probe['probe_bp'],
                'gene': probe['gene'],
                'beta': beta,
                'se': se,
                'pval': pval,
            })

        return associations

    def query_by_gene(self, gene_name: str) -> List[Dict]:
        """Query all associations for a specific gene.

        Args:
            gene_name: Gene name to query

        Returns:
            List of associations for the gene
        """
        # Find all probes for this gene
        probes_for_gene = [p for p in self.probes if p.get('gene') == gene_name]

        if not probes_for_gene:
            return []

        associations = []
        for probe in probes_for_gene:
            probe_idx = probe['row_idx']
            assocs = self.besd.get_probe_associations(probe_idx)

            for snp_idx, beta, se in assocs:
                snp = self.snp_by_idx[snp_idx]
                pval = calculate_p_value(beta, se)

                associations.append({
                    'snp_id': snp['snp_id'],
                    'snp_chr': snp['chr'],
                    'snp_bp': snp['bp'],
                    'a1': snp['a1'],
                    'a2': snp['a2'],
                    'trait_id': probe['probe_id'],
                    'trait_chr': probe['chr'],
                    'trait_bp': probe['probe_bp'],
                    'gene': probe['gene'],
                    'beta': beta,
                    'se': se,
                    'pval': pval,
                })

        return associations
