# Vector mode stores SE directly, not derived n_eff

When N is not constant across SNP-probe pairs we need per-pair storage. The alternatives were storing SE (float16) or deriving an effective sample size `n_eff = var_y / (se² × 2 × AF × (1−AF))` and storing that as float16. We store SE directly.

## Considered Options

**n_eff vector**: same 2 bytes/pair, but `n_eff` is not the actual sample size — it equals `n_actual / (1 − r²)`. For strong eQTLs (large explained variance, tiny SE), `n_eff` far exceeds the float16 ceiling of 65504, causing overflow and wrong reconstruction. It also requires AF at query time.

**SE vector**: SE values for typical eQTL/GWAS data fall well within float16 range (~0.001–10). Reconstruction is `beta = z × se` directly — no AF lookup needed at query time, no overflow risk.

## Consequences

Vector mode no longer requires AF in the ESI at either build or query time. The schema column is `se_vector BLOB` (not `n_vector`).
