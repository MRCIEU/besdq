"""Three-tier significance filter for GWAS-SSF import.

Tiers (applied in a single streaming pass):
  cis           – SNP within cis_radius of trait position (all p-values)
  sig_trans     – different chr (or outside cis), p < sig_threshold
  sug_trans     – different chr (or outside cis), sig_threshold <= p < sug_threshold
  dropped       – p >= sug_threshold outside cis
"""

from dataclasses import dataclass, field
from typing import List, Optional

from .gwas_ssf_reader import GwasSsfRow


@dataclass
class FilterResult:
    cis: List[GwasSsfRow] = field(default_factory=list)
    sig_trans_candidates: List[GwasSsfRow] = field(default_factory=list)
    sug_trans: List[GwasSsfRow] = field(default_factory=list)


def apply_significance_filter(
    rows,
    trait_chr: Optional[str],
    trait_bp: Optional[int],
    cis_radius: int = 1_000_000,
    sig_threshold: float = 5e-8,
    sug_threshold: float = 1e-4,
) -> FilterResult:
    """Classify GWAS rows into cis, sig_trans, and sug_trans buckets.

    Parameters
    ----------
    rows : iterable of GwasSsfRow
    trait_chr : chromosome of the trait (None → no cis tier)
    trait_bp : genomic position of the trait (None → no cis tier)
    cis_radius : bp distance defining the cis window
    sig_threshold : p-value threshold for genome-wide significance
    sug_threshold : p-value threshold for suggestive significance
    """
    result = FilterResult()
    has_cis = trait_chr is not None and trait_bp is not None

    for row in rows:
        in_cis = (
            has_cis
            and row.chr == trait_chr
            and abs(row.bp - trait_bp) <= cis_radius
        )

        if in_cis:
            result.cis.append(row)
        elif row.p < sig_threshold:
            result.sig_trans_candidates.append(row)
        elif row.p < sug_threshold:
            result.sug_trans.append(row)
        # else: dropped

    return result
