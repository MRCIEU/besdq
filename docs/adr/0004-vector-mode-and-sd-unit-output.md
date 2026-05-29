# ADR 0004: VectorN default for GWAS-SSF and SD-unit query output

**Status**: Accepted

## Context

Initial design of the GWAS-SSF import pipeline assumed ScalarN mode (store z per pair, derive SE from n + AF at query time) would be sufficient for single-cohort studies. In practice, SE reconstruction from `se = sqrt(trait_var / (n × 2 × AF × (1−AF)))` is imprecise when covariates are present: effective n varies per SNP, causing systematic scatter between reconstructed and original SE values. GWAS-SSF files always contain SE directly, making the approximation unnecessary.

Separately: BESDQ indexes are used to compare effects across hundreds of molecular phenotypes in a study, and across multiple studies. Original-unit betas are incomparable across traits (different phenotype scales) and across studies (different measurement instruments). SD-unit betas (beta / sd_y) place all effects on a common scale.

## Decisions

### 1. VectorN is the default storage mode for GWAS-SSF imports

SE is stored directly as float16 per SNP-trait pair (`se_vector BLOB`), alongside z (`zscores BLOB`). n is not stored for GWAS-SSF imports. Reconstruction at query time: `beta = z × se_stored / sd_y`.

ScalarN (store z + n_scalar, reconstruct SE from n and AF) is retained for legacy BESD imports where SE is not available directly.

### 2. SD-unit output is the default

All query methods return beta and SE in standard deviation units of the trait. `--original-scale` flag bypasses the sd_y division. z-scores and p-values are invariant between scales.

### 3. sd_y is auto-estimated from cis SNPs at build time

`trait_var = median( se_i² × n × 2 × eaf_i × (1−eaf_i) )` over cis SNPs. Stored in `epi.trait_var`. User may override via annotation TSV `trait_var` column. If no cis SNPs and no user-supplied value, returns original units with a warning.

## Alternatives considered

**ScalarN for GWAS-SSF**: Rejected. Demonstrated scatter between reconstructed and source SE showed the approximation breaks down under covariate adjustment. SE is available in GWAS-SSF — no reason to derive it.

**Store beta instead of z**: Considered when SE is stored directly (making z vs beta a pure choice). Rejected: z-scores are bounded (eQTL ±40, GWAS rarely >100) and well-scaled for float16 across all GWAS types. Beta spans orders of magnitude for GWAS complex traits (1×10⁻⁴ to 1×10⁰), causing float16 precision loss at small effect sizes. z is the more robust float16 storage unit.

**Original-unit default output**: Rejected. Effects are incomparable across traits in a 500-phenotype study without SD normalisation. SD units make the default output immediately usable for cross-trait and cross-study comparisons.

**Require user-supplied trait_var**: Rejected for usability. Auto-estimation from cis SNPs is robust (median estimator, <1% error under typical covariate adjustment). User override is available for exact control.

## Consequences

- SE reconstruction is exact (modulo float16 rounding) for GWAS-SSF imports
- Default query output requires trait_var to be populated — auto-estimated at build if cis SNPs are available
- Traits with no cis SNPs and no user-supplied trait_var return original-unit output with a warning
- VectorN storage is ~2× the per-pair cost of ScalarN (4 bytes vs 2 bytes per pair), acceptable given cis-window filtering
- ScalarN SD-unit SE simplifies elegantly: `se_sd = 1 / sqrt(n × 2 × AF × (1−AF))` — trait_var cancels
