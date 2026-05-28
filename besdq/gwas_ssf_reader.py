"""Streaming reader for gzipped GWAS-SSF files with allele normalisation."""

import gzip
from dataclasses import dataclass
from typing import Generator, Optional


@dataclass
class GwasSsfRow:
    chr: str
    bp: int
    a1: str  # alphabetically first allele
    a2: str  # other allele
    rsid: Optional[str]
    beta: float
    se: float
    eaf: float  # effect allele frequency (for a1)
    p: float

    @property
    def snp_key(self) -> str:
        return f"{self.chr}:{self.bp}:{self.a1}:{self.a2}"


def _parse_optional_rsid(val: str) -> Optional[str]:
    if not val or val in ('.', 'NA', 'nan', ''):
        return None
    return val


def read_gwas_ssf(path: str) -> Generator[GwasSsfRow, None, None]:
    """Stream rows from a gzipped GWAS-SSF file with allele normalisation.

    Allele convention: a1 is always the alphabetically first allele.
    When the effect_allele > other_allele, beta is negated and eaf inverted.
    """
    with gzip.open(path, 'rt') as fh:
        header_line = fh.readline()
        cols = header_line.rstrip('\n').split('\t')
        col_idx = {name: i for i, name in enumerate(cols)}

        chr_col = col_idx['chromosome']
        bp_col = col_idx['base_pair_location']
        ea_col = col_idx['effect_allele']
        oa_col = col_idx['other_allele']
        beta_col = col_idx['beta']
        se_col = col_idx['standard_error']
        eaf_col = col_idx['effect_allele_frequency']
        p_col = col_idx['p_value']

        # rsid may be in 'rsid' or 'rs_id' column
        rsid_col = col_idx.get('rsid', col_idx.get('rs_id'))

        for line in fh:
            parts = line.rstrip('\n').split('\t')
            if len(parts) < max(chr_col, bp_col, ea_col, oa_col, beta_col, se_col, eaf_col, p_col) + 1:
                continue

            try:
                effect_allele = parts[ea_col]
                other_allele = parts[oa_col]
                beta = float(parts[beta_col])
                se = float(parts[se_col])
                eaf = float(parts[eaf_col])
                p = float(parts[p_col])
            except (ValueError, IndexError):
                continue

            rsid = _parse_optional_rsid(parts[rsid_col]) if rsid_col is not None else None

            # Normalise: a1 must be alphabetically first
            if effect_allele <= other_allele:
                a1, a2 = effect_allele, other_allele
            else:
                # Swap: negate beta, invert eaf
                a1, a2 = other_allele, effect_allele
                beta = -beta
                eaf = 1.0 - eaf

            try:
                chromosome = parts[chr_col]
                bp = int(parts[bp_col])
            except (ValueError, IndexError):
                continue

            yield GwasSsfRow(
                chr=str(chromosome),
                bp=bp,
                a1=a1,
                a2=a2,
                rsid=rsid,
                beta=beta,
                se=se,
                eaf=eaf,
                p=p,
            )
