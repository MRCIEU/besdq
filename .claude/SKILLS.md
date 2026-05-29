# BESDQ Project Skills Reference

## Project Overview
`besdq` is a Python library for fast queries of BESD (Binary Efficient Sequential Data) eQTL summary statistics. It provides dual query modes: direct BESD file access and optimized SQLite indexing for repeated queries.

NOTE: The term `probe` may be generalised to `trait`.

## Architecture

### Hybrid Storage Design
- **Metadata & Indices**: SQLite tables (ESI, EPI) with indexed columns for O(log n) range queries
- **Statistics**: Per-probe BLOBs storing binary-serialized numpy arrays to preserve sequential read locality from sparse BESD format

### SQLite Schema
```sql
CREATE TABLE esi (
  row_idx INTEGER PRIMARY KEY,
  chr TEXT, snp_id TEXT, genetic_dist REAL, bp INTEGER, 
  a1 TEXT, a2 TEXT, freq REAL
);
CREATE INDEX idx_esi_chr_bp ON esi(chr, bp);
CREATE INDEX idx_esi_snp_id ON esi(snp_id);

CREATE TABLE epi (
  row_idx INTEGER PRIMARY KEY,
  chr TEXT, probe_id TEXT, genetic_dist REAL, probe_bp INTEGER, 
  gene TEXT, orientation TEXT
);
CREATE INDEX idx_epi_chr_bp ON epi(chr, probe_bp);
CREATE INDEX idx_epi_probe_id ON epi(probe_id);

CREATE TABLE probe_data (
  probe_idx INTEGER PRIMARY KEY,
  snp_count INTEGER,
  snp_indices BLOB,  -- int32 numpy array
  betas BLOB,        -- float32 numpy array
  ses BLOB           -- float32 numpy array
);

CREATE TABLE besd_meta (
  key TEXT PRIMARY KEY,
  value TEXT
);
```

## BESD File Format (Sparse)

### Magic Numbers
- **Sparse**: `0x40400000` or `3` (SPARSE_FILE_TYPE_3F / SPARSE_FILE_TYPE_3)
- **Dense**: Not supported by this parser implementation

Source: [SMR source code](https://yanglab.westlake.edu.cn/software/smr/#Sourcecode) `src/SMR_data_p1.cpp`

### File Structure
- [uint32 magic]
  - The magin number indicating file type, as stated above
- [15 x uint32 reserved]
  - Only present in type 3 but not 3F
- [uint64 val_num]
  - Number of values
- [uint64 cols[2 * n_probes + 1]]
  - Boundaries of list of values of each probe
- [uint32 rowid[val_num]]
  - SNP IDs that the value in the list belongs to
- [float32 val[val_num]]
  - Blocks of values

### Illustrative Example (type 3 and 3F)

|      | probe0<br>beta | probe0<br>se | probe1<br>beta | probe1<br>se |
| ---- | -------------- | ------------ | -------------- | ------------ |
| snp0 | 0.2            | 0.05         |                |              |
| snp1 |                |              | 0.3            | 0.07         |
| snp2 | -0.1           | 0.04         |                |              |

Reformat - read by columns left-right then rows top-bottom:

| index | 0    | 1    | 2    | 3    | 4    | 5    |
| ----- | ---- | ---- | ---- | ---- | ---- | ---- |
| value | 0.2  | -0.1 | 0.05 | 0.04 | 0.3  | 0.07 |
| probe | 0    | 0    | 0    | 0    | 1    | 1    |
| field | beta | beta | se   | se   | beta | se   |
| snp   | 0    | 2    | 0    | 2    | 1    | 1    |

- `val_num = 6`
  - There are six values
- `cols[5] = [0, 2, 4, 5, 6]`
  - Column (field) boundaries are `[0, 2)` `[2, 4)` `[4, 4)` `[5, 5)` and `6`
- `row_id[6] = [0, 2, 0, 2, 1, 1]`
  - SNP IDs
- `val[6] = [0.2, -0.1, 0.05, 0.04, 0.3, 0.07]`

### Reconstruct from BESD
The `i`-th association of the `p`-th probe can be reconstructed from global arrays using `cols` boundaries:
1. `beta_start = cols[2p]`
2. `se_start = cols[2p+1]`
3. `association_count = beta_start - se_start`
4. `snp_idx(i) = rowid[beta_start + i]`
5. `beta(i) = val[beta_start + i]`
6. `se(i) = val[se_start + i]`

**Critical Detail**: `val_num` is the global payload length used to read `rowid` and `val` arrays; per-probe association counts are derived from `cols` boundaries.

### Row Indexing
- ESI/EPI row order must be preserved exactly
- `row_idx` is zero-based in this codebase and assigned by row order in index readers
- Genomic bp values are consumed as provided by the source index files

## Query Modes

### 1. Direct BESD File Queries (`BESDQueryEngine`)
- **Speed**: ~100-200ms per query
- **Memory**: Indices and BESD arrays are loaded during reader initialization
- **Use Case**: One-time queries, exploratory analysis
- **No Setup**: Works immediately on BESD files

### 2. SQLite Index Queries (`BESDQueryIndex`)
- **Speed**: ~10-50ms per query (after index creation)
- **Space**: 50-70% of original BESD file size
- **Setup**: ~2 seconds to create index (one-time)
- **Use Case**: Large-scale analysis, repeated queries

### Supported Query Types
1. **Cis-window**: SNP region + probe region simultaneously
2. **SNP ID**: All associations for specific SNP(s)
3. **Probe ID**: All associations for specific probe(s)
4. **Gene**: All associations for specific gene(s)

### Coordinate Formats
- **Kilobase format**: `--snp-chr 1 --from-snp-kb 100 --to-snp-kb 2000`
- **Base pair format**: `--snp-chr 1 --from-snp-bp 100000 --to-snp-bp 2000000`
- **Chr:pos format**: `--snp-chrpos 1:100000-2000000` or single position `1:1191870`
- **Comma-separated lists**: Support for batch queries with ranges or single positions

## Package Structure

```
besdq/
├── __init__.py              # Public API exports
├── besd_reader.py           # BESDQueryEngine (direct BESD file access)
│   ├── BESDQueryEngine      # Main query interface
│   ├── BESDReader           # Binary BESD file parser
│   ├── IndexReader          # .esi/.epi file reader
│   └── norm_cdf()           # P-value calculation via Abramowitz-Stegun
├── sqlite_query.py          # BESDQueryIndex (SQLite queries)
│   ├── BESDQueryIndex       # Query interface with context manager
│   └── get_probe_snps()     # Deserialize numpy BLOBs
├── builder.py               # BESDIndexBuilder (database creation)
│   ├── build()              # Main build workflow
│   ├── _create_schema()     # SQLite schema setup
│   └── _write_probe_data()  # Serialize numpy arrays to BLOBs
├── queries.py               # Unified interface (re-exports both engines)
└── cli.py                   # Command-line interface
```

## Key Implementation Details

### P-value Calculation
- Formula: `pval = 2 * (1 - norm_cdf(|beta/se|))`
- Uses Abramowitz and Stegun polynomial approximation (~0.00012 accuracy)
- Handles edge cases: SE=0, extreme z-scores (underflow to 0 is expected)

### Numpy Array Serialization
- Use `tobytes()` for serialization, `frombuffer()` for deserialization
- Stored with explicit dtype specification: `int32` for indices, `float32` for betas/SEs
- **Important**: Convert numpy.int32 to Python int for SQLite parameterized queries

### Output Format (SMR-Compatible)
Tab-separated with columns:
```
SNP  SNP_Chr  SNP_bp  A1  A2  Probe  Probe_Chr  Probe_bp  Gene  Beta  SE  P_value
```

### CLI Features
- **Mutual exclusive input**: `--beqtl-summary` (BESD files) OR `--besd-index` (SQLite)
- **Index creation**: `--index` flag with `--beqtl-summary`
- **Data source indication**: Output clearly states which source is used
- **P-value filtering**: `--query` threshold applied before output

## Common Patterns

### Direct Query (No Database)
```python
from besdq import BESDQueryEngine

engine = BESDQueryEngine('data/westra')
assocs = engine.query_cis_window(
    snp_chr='1', snp_start_kb=100, snp_end_kb=2000,
    probe_chr='1', probe_start_kb=1000, probe_end_kb=2000
)
```

### SQLite Query
```python
from besdq import BESDQueryIndex

with BESDQueryIndex('data/westra.db') as index:
    assocs = index.query_by_snp_id('rs3818646')
```

### Create Index
```python
from besdq import BESDIndexBuilder

builder = BESDIndexBuilder('data/westra.db')
builder.build('data/westra', force=True)
```

## Testing

**26 comprehensive unit tests** covering:
- Data loading and format detection
- All query types and coordinate formats
- P-value calculation accuracy
- Consistency between BESD and SQLite modes
- Edge cases (empty queries, boundary conditions)

**Test Data**: Westra eQTL (506,049 SNPs, 5,966 probes)

Run tests:
```bash
python3 -m unittest tests.test_queries -v
```

## Performance Characteristics

| Metric         | Direct BESD | SQLite Index |
| -------------- | ----------- | ------------ |
| Query time     | 100-200ms   | 10-50ms      |
| Index creation | N/A         | ~2s          |
| File size      | ~1.2GB      | ~400-600MB   |
| Memory         | ~500MB      | ~100MB       |

## Error Handling & Edge Cases

- **Zero SE**: norm_cdf skips calculation, returns 0
- **Extreme z-scores**: P-values underflow to 0 (expected for very significant associations)
- **Empty queries**: Return empty list (not an error)
- **Numpy/SQLite compatibility**: Convert numpy.int32 to Python int for parameterized queries
- **Floating point precision**: Use `places=4` or `places=5` in float comparisons

## Dependencies

**Required:**
- Python 3.9+
- numpy>=1.20 (efficient array serialization)

**Development (optional):**
- black>=22.0 (code formatting)
- flake8>=4.0 (linting)

**No external dependencies** for core BESD querying (without SQLite indexing)

## Known Limitations & Future Enhancements

1. Result ordering differs between BESD and SQLite (both contain same data)
2. P-values can underflow to 0 for extremely significant associations
3. Future: Support for compressed BESD formats
4. Future: Parallel batch queries
5. Future: Output format variants (JSON, Parquet, VCF)
6. Future: Web API for remote queries
