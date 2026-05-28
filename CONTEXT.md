# BESDQ Domain Glossary

## Statistics Encoding

The representation used in `probe_data` BLOBs. Two modes:

One format only: **z-score encoding**. Stores float16 `zscores` per SNP-probe pair and `n` (scalar per probe in ScalarN mode, float16 BLOB per pair in VectorN mode). Reconstructs beta and se on the fly at query time. No verbatim encoding path.

**Schema changes for z-score encoding:**
- `probe_data` table: `betas`/`ses` (float32 BLOBs) replaced by `zscores` (float16 BLOB); gains `n_scalar INTEGER` (ScalarN mode) and `n_vector BLOB` (VectorN mode, float16 array aligned to `snp_indices`)
- `epi` table: gains optional `var_y REAL` column (trait variance per probe, default 1)
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

## Trait Variance (var_y)

Variance of the phenotype (e.g. gene expression level) for a probe. Per-probe quantity. Supplied at build time as a separate two-column TSV file (`probe_id  var_y`) via `--trait-variance`. Stored in the `epi` SQLite table as `var_y REAL`. When absent, assumed to be 1 with a build-time warning. Needed for exact reconstruction when phenotypes are not standardised.

## Allele Frequency (AF)

Minor allele frequency of a SNP, used in reconstruction. Sourced from the `.esi` file (column 7, optional). When AF is absent from ESI and ScalarN mode is used, AF is derived per SNP during build via:

```
AF = (1 − sqrt(1 − 2 × var_y / (n × se²))) / 2
```

and written into the ESI table of the Lean Index.

## Significance Mask

Out of scope for the current work. BESD files arrive pre-filtered by their own toolchain (e.g. eQTLgen applies p < 1e-5 cis, p < 5e-8 trans). The Lean Index stores everything in the source BESD verbatim — it does not apply additional filters.

A significance mask may become relevant when direct import from VCF or text files is added (future work), at which point cis/trans thresholds will need to be defined explicitly.

## Build Modes (Lean Index)

| Dataset type | AF in ESI | n source | Blocked? |
|---|---|---|---|
| Single cohort, no AF | No | User-supplied scalar (`--sample-size`) | No — AF derived and written to ESI |
| Single cohort, with AF | Yes | Derived per probe from (se, AF, var_y) | No |
| Meta-analysis, with AF | Yes | Derived per pair from (se, AF, var_y); VectorN mode | No |
| Meta-analysis, no AF | No | Unresolvable | **Yes — error at build time** |
