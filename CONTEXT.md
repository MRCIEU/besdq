# BESDQ Domain Glossary

## Trait

A molecular phenotype being measured — e.g. gene expression of IL10, a CpG methylation level, a protein abundance. The generalisation of "probe" (the BESD binary format's term). Each BESDQ dataset contains one or more traits, each with its own row in the `epi` table and its own statistics block in `probe_data`.

**EPI columns by class:**

| Class | Columns | Behaviour when absent |
|---|---|---|
| Mandatory | `trait_id TEXT`, `trait_name TEXT` | Import fails |
| Optional functional | `trait_var REAL`, `trait_chr TEXT`, `trait_bp INTEGER` | Affects reconstruction / cis-window; warned at build time |
| Optional non-consequential | `gene TEXT`, `context TEXT` | NULL; no effect on queries or reconstruction |

`gene`: human-readable gene symbol (e.g. IL10). `context`: free-text descriptor covering tissue, treatment, time-point, sex, or any other experimental condition that distinguishes traits within a study (e.g. "PBMC_Bbmix_baseline").

`orientation` is dropped — it was carried from the BESD `.epi` convention but is unused in any query or output.

Study-level metadata (publication, year, ancestry, max sample size, study type, source URL) is stored as a single JSON blob in the `besd_meta` table under the key `study_metadata`. This is dataset-wide, not per-trait.

> **Probe** in the BESD binary format is the historical synonym for Trait. The `.epi` file, `probe_data` table, and binary block structure still use "probe" internally.

## SNP Identity and Allele Convention

The canonical SNP key is `chr:pos:A1:A2` where **A1 is always the alphabetically first allele** and A2 is the other. This convention is shared with pleioDB and makes cross-dataset matching on `chr_pos_a1_a2` immediately interpretable without allele harmonisation.

At import time, every SNP-trait association is normalised: if the source file's effect allele is not alphabetically first, the alleles are swapped and the beta (and z-score) is negated. This happens unconditionally for all sources (GWAS-SSF, BESD, any future format).

ESI deduplication uses `chr:pos:A1:A2` as the unique SNP key. The rsid is stored as a lookup field (`snp_id`) but is not the primary key. Novel imputed variants without an rsid are stored with `snp_id = NULL`.

## Statistics Encoding

The representation used in `probe_data` BLOBs. Two modes, selected by data source:

**VectorN mode** (default for GWAS-SSF and any source that provides SE directly): stores `zscores` (float16) + `se_vector` (float16, SE in original study units) per SNP-trait pair. No AF or n needed at query time. SE is stored in original units and divided by `sd_y = sqrt(trait_var)` at query time to produce SD-unit output.

**ScalarN mode** (legacy BESD imports): stores `zscores` (float16) per pair + `n_scalar INTEGER` per trait. SE is reconstructed from n, AF, and trait_var. SD-unit SE simplifies to `1 / sqrt(n × 2 × AF × (1−AF))` — trait_var cancels.

**`probe_data` schema:**
- `zscores BLOB` — float16 numpy array, always present
- `se_vector BLOB` — float16 numpy array, VectorN only (NULL in ScalarN)
- `n_scalar INTEGER` — ScalarN only (NULL in VectorN)

**`epi` schema:** mandatory `trait_id TEXT`, `trait_name TEXT`; optional functional `trait_var REAL`, `trait_chr TEXT`, `trait_bp INTEGER`; optional non-consequential `gene TEXT`, `context TEXT`.

**`esi` schema:** `freq` column populated during build if absent in source (derived from beta, se, user-supplied n in ScalarN mode).

## Z-Score

`z = beta / se`. The primary stored statistic in the Lean Index. Stored as float16 (sufficient precision for eQTL z-scores, which rarely exceed ±40).

## Reconstruction

The process of deriving output beta and se from stored quantities at query time. Default output is **SD units** (see Output Scale). Two paths:

**VectorN path** (GWAS-SSF imports):
```
se_out  = se_stored / sd_y        # sd_y = sqrt(trait_var) from EPI
beta_out = z × se_out
```
Requires: z (stored), se_stored (stored), sd_y from EPI `trait_var`.

**ScalarN path** (legacy BESD imports):
```
se_orig = sqrt(trait_var / (n × 2 × AF × (1 − AF)))
se_out  = se_orig / sd_y = 1 / sqrt(n × 2 × AF × (1 − AF))   # trait_var cancels
beta_out = z × se_out
```
Requires: z (stored), n (stored), AF (stored in ESI). `trait_var` cancels — SD-unit SE is independent of phenotype scale.

With `--original-scale`: `se_out = se_stored` (VectorN) or `se_out = se_orig` (ScalarN); no division by sd_y.

## N (Sample Size)

The number of individuals used to compute each association.

- **ScalarN mode**: n is constant across all SNPs for a trait (typical single-cohort BESD eQTL datasets). Stored once per trait in `probe_data` as `n_scalar INTEGER`. Supplied by the user via `--sample-size N` or read from YAML.
- **VectorN mode**: SE is stored directly per SNP-trait pair (GWAS-SSF imports and any source that provides SE). n is not stored. `se_vector BLOB` (float16) aligned to `snp_indices`.

## Trait Variance (trait_var)

Variance of the phenotype (e.g. gene expression level) for a trait. Per-trait quantity stored in `epi.trait_var REAL`. Used to compute `sd_y = sqrt(trait_var)` for SD-unit output.

**Sources, in precedence order:**
1. User-supplied in annotation TSV `trait_var` column (exact)
2. Auto-estimated at build time from cis SNPs using the median formula:
   `trait_var = median( se_i² × n × 2 × eaf_i × (1 − eaf_i) )` over all cis SNPs
3. 1.0 if no cis SNPs and not user-supplied — output is returned in original units with a warning

Auto-estimation is robust to per-SNP missing data (median estimator) and approximately correct under covariate adjustment (effective-n error ≈ k/n, typically < 1%). Previously named `var_y`.

## Output Scale

All query methods return beta and SE in **SD units** (standard deviation of the trait) by default, making effect sizes directly comparable across studies and traits. Pass `--original-scale` to return beta and SE in original study units.

SD-unit conversion: `beta_sd = beta_orig / sd_y`, `se_sd = se_orig / sd_y`, where `sd_y = sqrt(trait_var)` from EPI. The z-score and p-value are invariant — they are identical in both scales.

## Allele Frequency (AF)

Minor allele frequency of a SNP, used in reconstruction. Sourced from the `.esi` file (column 7, optional). When AF is absent from ESI and ScalarN mode is used, AF is derived per SNP during build via:

```
AF = (1 − sqrt(1 − 2 × var_y / (n × se²))) / 2
```

and written into the ESI table of the Lean Index.

## Significance Mask

The filter applied during import from text files (e.g. GWAS-SSF) to determine which SNP-trait associations are stored. Three tiers:

| Tier | Condition | What is stored |
|---|---|---|
| Cis | SNP within cis-radius of `trait_chr`/`trait_bp` | All variants unconditionally |
| Significant trans | p < 5×10⁻⁸ outside cis | All variants within sig-radius of each independent lead SNP |
| Suggestive trans | 5×10⁻⁸ ≤ p < 1×10⁻⁴ outside cis | That variant only |
| Below suggestive | p ≥ 1×10⁻⁴ outside cis | Dropped |

**Cis-radius**: the distance from `trait_bp` within which all variants are retained unconditionally. Default 1,000,000 bp (±1 Mb). Requires `trait_chr` and `trait_bp`; if absent, cis tier is skipped and only trans tiers apply.

**Sig-radius**: the distance from a significant trans lead SNP within which all variants are retained. Default 500,000 bp (±500 kb).

Independent significant trans peaks are identified by LD clumping with plink2 (must be on PATH) and a user-supplied plink2-format LD reference panel (`--pfile` prefix). Clumping defaults: r²=0.01, window=10,000 kb. Thresholds and radii are configurable at import time.

BESD files imported via the legacy path arrive pre-filtered and are stored verbatim — the Lean Index applies no additional filters to them.

## Build Modes (Lean Index)

| Dataset type | Storage mode | SE source | sd_y source |
|---|---|---|---|
| GWAS-SSF (any source with SE) | VectorN | Direct from file (float16) | Auto-estimated from cis SNPs; user override via `trait_var` |
| Single-cohort BESD, no AF | ScalarN | Reconstructed from n, AF (derived), trait_var | User-supplied `trait_var` or 1.0 |
| Single-cohort BESD, with AF | ScalarN | Reconstructed from n, AF, trait_var | User-supplied `trait_var` or 1.0 |
| Meta-analysis BESD, with AF | ScalarN (per-pair n derived at query time) | Reconstructed | User-supplied `trait_var` or 1.0 |
| Meta-analysis BESD, no AF | — | Unresolvable | **Error at build time** |
