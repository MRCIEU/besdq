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

The representation used in `probe_data` BLOBs. Two modes:

One format only: **z-score encoding**. Stores float16 `zscores` per SNP-trait pair and `n` (scalar per trait in ScalarN mode, float16 BLOB per pair in VectorN mode). Reconstructs beta and se on the fly at query time. No verbatim encoding path.

**Schema changes for z-score encoding:**
- `probe_data` table: `betas`/`ses` (float32 BLOBs) replaced by `zscores` (float16 BLOB); gains `n_scalar INTEGER` (ScalarN mode) and `n_vector BLOB` (VectorN mode, float16 array aligned to `snp_indices`)
- `epi` table: mandatory `trait_id TEXT`, `trait_name TEXT`; optional functional `trait_var REAL` (default 1), `trait_chr TEXT`, `trait_bp INTEGER`; optional non-consequential stored as JSON in `metadata TEXT`
- `esi` table: `freq` column populated during build if absent in source (derived from beta, se, user-supplied n)

## Z-Score

`z = beta / se`. The primary stored statistic in the Lean Index. Stored as float16 (sufficient precision for eQTL z-scores, which rarely exceed ±40).

## Reconstruction

The process of deriving beta and se from stored quantities at query time. Formula:

```
se  = sqrt(var_y / (n × 2 × AF × (1 − AF)))
beta = z × se
```

Requires: z (stored per pair), n (stored per probe or per pair), AF (stored in ESI), var_y (stored per probe in EPI, default 1).

## N (Sample Size)

The number of individuals used to compute each association. May be stored as a **scalar** (one value per probe, ScalarN mode) or a **vector** (one value per SNP-probe pair, VectorN mode).

- **ScalarN mode**: n is constant across all SNPs for a probe (typical single-cohort eQTL datasets). Stored once per probe in `probe_data` as a single INTEGER. Supplied by the user via `--sample-size N`. Selected explicitly via `--n-mode scalar`.
- **VectorN mode**: n varies per SNP-probe pair (meta-analyses such as eQTLgen, GoDMC). Stored as a float16 BLOB per probe in `probe_data`, aligned to `snp_indices`. Selected explicitly by the user via `--n-mode vector`. Requires AF to be present in ESI.

## Trait Variance (trait_var)

Variance of the phenotype (e.g. gene expression level) for a trait. Per-trait quantity. Stored in the `epi` SQLite table as `trait_var REAL`. When absent, assumed to be 1 with a build-time warning. Needed for exact reconstruction when phenotypes are not standardised. Previously named `var_y`.

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

| Dataset type | AF in ESI | n source | Blocked? |
|---|---|---|---|
| Single cohort, no AF | No | User-supplied scalar (`--sample-size`) | No — AF derived and written to ESI |
| Single cohort, with AF | Yes | Derived per trait from (se, AF, trait_var) | No |
| Meta-analysis, with AF | Yes | Derived per pair from (se, AF, trait_var); VectorN mode | No |
| Meta-analysis, no AF | No | Unresolvable | **Yes — error at build time** |
