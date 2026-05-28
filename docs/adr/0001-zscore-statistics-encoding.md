# Z-score encoding for statistics storage

The SQLite index stores float16 z-scores instead of float32 beta and SE per SNP-probe pair, with beta/SE reconstructed at query time. Two modes: **ScalarN** stores z (2 bytes) + n per probe (integer, amortised); **vector** stores z (2 bytes) + SE (2 bytes) per pair. This halves stats storage (12 → 6 bytes/pair in ScalarN mode, 12 → 4 bytes/pair in vector mode) with negligible precision loss for downstream consumers.

## Considered options

**Verbatim float32 beta + SE** — exact, no reconstruction needed, but 12 bytes/pair and no path to further compression.

**Float16 beta + float16 SE (drop z-score)** — 50% savings on stats storage, exact to float16 precision, but requires AF to be in the ESI for no additional benefit over z-score encoding.

**Z-score + n scalar or SE vector (chosen)** — z-scores are dimensionless and stable across datasets. When n is constant per probe (ScalarN mode), only z is stored per pair and n is a single integer per probe. When n varies per pair (vector mode), SE is stored directly per pair — avoids deriving n_eff which overflows float16 for strong effects and still requires AF at query time. See ADR-0002 for the SE-vs-n_eff decision.

## Consequences

- **ScalarN reconstruction** requires AF per SNP (from ESI) and optionally var_y per probe: `se = sqrt(var_y / (n × 2 × AF × (1−AF)))`, `beta = z × se`. When AF is absent from the source ESI and `--n-mode scalar` is used, AF is derived at build time from `(se, n)` assuming `var(y) = 1` and written into the ESI table with a warning.
- **Vector reconstruction** needs only z and SE: `beta = z × se`. No AF lookup at query time.
- Two modes: **ScalarN** (`--n-mode scalar`, n stored as INTEGER per probe); **vector** (`--n-mode vector`, SE stored as float16 BLOB `se_vector` per pair).
- Only z-score encoding is supported. There is no verbatim float32 path; existing databases must be rebuilt.
- Future direct import from VCF/text will need explicit cis/trans thresholds defined at that point; this ADR does not address that case.
