# ADR 0003: Three-tier significance mask for GWAS-SSF import

**Status**: Accepted

## Context

GWAS-SSF files from EBI contain genome-wide summary statistics — typically ~7.7 million SNP-trait associations per file. A study may comprise 500 such files (one per molecular phenotype). Storing all associations verbatim is impractical (~3.85 billion rows). A filtering strategy is needed that preserves everything analytically relevant while reducing storage by orders of magnitude.

The primary downstream use is SMR-style cis-QTL analysis and trans-QTL lookup, both of which require different completeness guarantees in different genomic regions.

## Decision

Associations are classified into three tiers at import time, applied independently per trait:

| Tier | Condition | Stored |
|---|---|---|
| Cis | SNP within cis-radius of `trait_bp` | All variants unconditionally |
| Significant trans | p < 5×10⁻⁸ outside cis | All variants within sig-radius of each independent lead SNP |
| Suggestive trans | 5×10⁻⁸ ≤ p < 1×10⁻⁴ outside cis | That variant only |
| Below suggestive | p ≥ 1×10⁻⁴ | Dropped |

Default cis-radius: 1,000,000 bp. Default sig-radius: 500,000 bp. All thresholds and radii are configurable at import time.

Independent significant trans peaks are identified by LD clumping (plink2, r²=0.01, 10,000 kb window) using a user-supplied LD reference panel matched to study ancestry. plink2 must be on PATH.

The cis tier is skipped when `trait_chr`/`trait_bp` are absent from the annotation TSV (trans-only mode). Traits without genomic coordinates are valid inputs.

## Alternatives considered

**P-value threshold only (no cis/trans distinction)**: Simpler, no need for probe coordinates. Rejected because cis regions require complete coverage for SMR — a significance filter inside cis would discard associations needed to model the full LD structure of a locus.

**Shape-based peak detection (no external tool)**: Self-contained, no plink2 dependency. Rejected because it cannot resolve two independent signals within the same genomic neighbourhood without LD information — exactly the case that matters for trans-QTL colocalisation.

**Store everything above suggestive, no window expansion**: Would miss variants in significant trans peaks that fall below the suggestive threshold but are needed for regional LD structure. Rejected for the same reason as P-value-only.

## Consequences

- plink2 must be on PATH at import time; import fails if absent and significant trans associations are present.
- Users must supply an LD reference panel in plink2 `--pfile` format matched to study ancestry.
- Cis completeness is guaranteed; trans coverage is bounded and reproducible given the same thresholds.
- Re-importing with different thresholds requires reprocessing all source files.
