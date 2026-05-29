# BESDQ - Fast BESD eQTL Query Tool

Fast queries of BESD (Binary Efficient Sequential Data) eQTL summary statistics files without requiring a database. Supports both SPARSE_FILE_TYPE_3 and SPARSE_FILE_TYPE_3F formats and provides SMR-compatible command-line interface.

## Credits

This work builds on the BESD format designed and implemented in the [SMR software](https://yanglab.westlake.edu.cn/software/smr/#eQTLsummarydata) from [Jian Yang's group](https://yanglab.westlake.edu.cn/).

## Installation

### From GitHub (Recommended)

```bash
git clone https://github.com/explodecomputer/besdq.git
cd besdq
pip install -e .
```

This installs the package in development mode with the `besdq` command-line tool.

### From PyPI (when available)

```bash
pip install besdq
```

### Running without installation

You can also run the CLI directly without installing:

```bash
python3 -m besdq.cli --help
```

## Dependencies

**Required:**
- Python 3.9+
- numpy>=1.20 (for efficient array storage in database indices)
- PyYAML>=5.1 (for reading EBI GWAS Catalog metadata files)

**Development (optional):**
- black>=22.0 (code formatting)
- flake8>=4.0 (linting)

## System Requirements

- ~1GB free disk space per BESD dataset
- ~50-70% additional space for SQLite index databases (optional but recommended)

## Building an Index from GWAS-SSF Files (EBI GWAS Catalog)

If you have summary statistics in [GWAS-SSF format](https://www.ebi.ac.uk/gwas/docs/summary-statistics-format) (e.g. downloaded from the EBI GWAS Catalog), use the `import-gwas-ssf` command to build a queryable index directly — no intermediate conversion required.

### 1. Prepare a trait annotation file

Create a tab-separated file (e.g. `traits.tsv`) with one row per summary-statistics file. Required columns are `file_path`, `trait_id`, and `trait_name`. Provide `trait_chr`/`trait_bp` to enable cis-region storage; omit them for trans-only mode.

```
file_path	trait_id	trait_name	trait_chr	trait_bp	gene	context
data/GCST90275731.h.tsv.gz	GCST90275731	IL10 expression PBMC Bbmix 1e-5	1	206774541	IL10	PBMC_Bbmix_baseline
data/GCST90275739.h.tsv.gz	GCST90275739	IL1Ra expression PBMC Bbmix 1e-4	2	113099315	IL1RN	PBMC_Bbmix_baseline
```

If `sample_size` is omitted from the TSV, it is read automatically from the companion `*-meta.yaml` file that the EBI GWAS Catalog provides alongside each summary-statistics file.

Optional columns:

| Column | Default | Description |
|---|---|---|
| `trait_chr` | — | Chromosome of the trait locus (enables cis storage) |
| `trait_bp` | — | Base-pair position of the trait locus |
| `sample_size` | from YAML | Scalar sample size N |
| `trait_var` | 1.0 | Phenotype variance (for reconstruction) |
| `gene` | — | Gene symbol |
| `context` | — | Free-text experimental context |

### 2. Run the import

```bash
import-gwas-ssf \
  --trait-annotation data/ebi_input/traits.tsv \
  --ld-reference data/ldref/EUR \
  --output data/ebi_import/study.db
```

**Significance filtering applied during import:**

| Tier | Condition | Stored |
|---|---|---|
| Cis | SNP within ±1 Mb of `trait_chr`/`trait_bp` | All variants |
| Significant trans | p < 5×10⁻⁸, different chromosome | ±500 kb window around each independent peak (LD-clumped) |
| Suggestive trans | 5×10⁻⁸ ≤ p < 1×10⁻⁴ | That variant only |
| Below suggestive | p ≥ 1×10⁻⁴ | Dropped |

LD clumping requires [plink2](https://www.cog-genomics.org/plink/2.0/) on `PATH` and a plink2-format reference panel (`--ld-reference` prefix). Install plink2 with:

```bash
conda install -c bioconda plink2
```

**Key options:**

```
--trait-annotation FILE   Trait annotation TSV (required)
--ld-reference PREFIX     plink2 LD reference prefix (required)
--output FILE             Output database (default: <tsv_stem>.db)
--workers N               Parallel workers for Pass 1 (default: 1)
--cis-radius BP           Cis window radius (default: 1000000)
--sig-threshold P         Genome-wide significance (default: 5e-8)
--sug-threshold P         Suggestive threshold (default: 1e-4)
--sig-radius BP           Trans peak window radius (default: 500000)
--clump-r2 R2             LD clumping r² threshold (default: 0.01)
--clump-kb KB             LD clumping window (default: 10000)
```

### 3. Query the resulting index

Once built, the index is queried with the standard `besdq` tool:

```bash
# Query by SNP
besdq --besd-index study.db --snp rs12238997 --out results/out

# Query by trait ID
besdq --besd-index study.db --probe GCST90275731 --out results/out

# Query by genomic region
besdq --besd-index study.db \
  --snp-chrpos 1:206000000-208000000 \
  --probe-chrpos 1:205000000-208000000 \
  --out results/out
```

---

## Basic Usage

### Two Query Options

BESDQ supports querying data from two sources:

1. **Direct BESD files** - `--beqtl-summary`
   - Reads from binary `.besd`, `.esi`, `.epi` files directly
   - No setup required, works immediately
   - Suitable for one-time queries or exploratory analysis

2. **SQLite index** - `--besd-index`
   - Queries from an indexed SQLite database
   - Requires one-time database creation (fast)
   - Optimized for repeated queries and large-scale analysis
   - 50-70% smaller than original BESD files

### Creating an SQLite Index Database

For improved performance on repeated queries, create a SQLite index database (one-time setup):

```bash
# Create index from BESD files
besdq --beqtl-summary data/westra_eqtl_hg19 --index data/westra_eqtl_hg19.db
```

This creates a database with:
- SQLite tables for SNP and probe metadata with indexed columns for fast range queries
- Per-probe BLOBs storing numpy arrays of associations (snp_indices, betas, SEs)
- Full metadata preservation from the original BESD files

### Command-line Interface

#### Option 1: Query BESD Files Directly

```bash
# Using chr:pos format (recommended)
time besdq --beqtl-summary data/westra_eqtl_hg19 \
  --snp-chrpos 1:100000-2000000 \
  --probe-chrpos 1:1000000-2000000 \
  --out results/output \
  --query 1e-4
```

#### Option 2: Query SQLite Index (Faster for Repeated Queries)

```bash
# Using SQLite index
time besdq --besd-index data/westra_eqtl_hg19.db \
  --snp-chrpos 1:100000-2000000 \
  --probe-chrpos 1:1000000-2000000 \
  --out results/output \
  --query 1e-4
```

#### Coordinate Format Options

Both query modes support the same coordinate format options:

**Kilobase format:**
```bash
besdq --beqtl-summary data/westra_eqtl_hg19 \
  --snp-chr 1 --from-snp-kb 100 --to-snp-kb 2000 \
  --probe-chr 1 --from-probe-kb 1000 --to-probe-kb 2000 \
  --out results/output
```

**Base pair format:**
```bash
besdq --beqtl-summary data/westra_eqtl_hg19 \
  --snp-chr 1 --from-snp-bp 100000 --to-snp-bp 2000000 \
  --probe-chr 1 --from-probe-bp 1000000 --to-probe-bp 2000000 \
  --out results/output
```

**Chr:pos format (range or single position):**
```bash
# Range query
besdq --beqtl-summary data/westra_eqtl_hg19 \
  --snp-chrpos 1:100000-2000000 \
  --probe-chrpos 1:1000000-2000000 \
  --out results/output

# Single position query
besdq --beqtl-summary data/westra_eqtl_hg19 \
  --snp-chrpos 1:1191870 \
  --probe-chrpos 1:1140818 \
  --out results/output
```

Validation notes:
- For `chr:start-end`, `start` must be less than or equal to `end`.
- `--snp`, `--probe`, and `--gene` are mutually exclusive query modes.
- Identifier query modes (`--snp`, `--probe`, `--gene`) cannot be combined with region options (`--snp-chr*`, `--probe-chr*`, `--from-*`, `--to-*`).

#### P-value Filtering

Filter results by p-value threshold (applies to both query modes):

```bash
# Query with p-value filter
besdq --beqtl-summary data/westra_eqtl_hg19 \
  --snp-chrpos 1:100000-2000000 \
  --probe-chrpos 1:1000000-2000000 \
  --query 1e-4 \
  --out results/output
```

### Python API

BESDQ provides two query engines with identical interfaces:

#### Option 1: Query BESD Files Directly

```python
from besdq import BESDQueryEngine

# Initialize with BESD file prefix
engine = BESDQueryEngine('data/westra_eqtl_hg19')

# Query associations
associations = engine.query_cis_window(
    snp_chr='1', snp_start_kb=100, snp_end_kb=2000,
    probe_chr='1', probe_start_kb=1000, probe_end_kb=2000,
)

# Results include SNP/trait metadata and statistics
for assoc in associations:
    print(f"{assoc['snp_id']} - {assoc['trait_id']}: "
          f"beta={assoc['beta']:.4f}, p={assoc['pval']:.2e}")
```

#### Option 2: Query SQLite Index

```python
from besdq import BESDQueryIndex

# Open index database (context manager handles connection)
with BESDQueryIndex('data/westra_eqtl_hg19.db') as index:
    # Query by cis-window (same interface as BESDQueryEngine)
    associations = index.query_cis_window(
        snp_chr='1', snp_start_kb=100, snp_end_kb=2000,
        probe_chr='1', probe_start_kb=1000, probe_end_kb=2000,
    )
    
    # Additional query methods for SQLite index
    # Query all associations for a probe
    assocs = index.query_by_probe_id('ILMN_2349633')
    
    # Query all associations for an SNP
    assocs = index.query_by_snp_id('rs3818646')
```

## Output Format

Results are written in tab-separated format compatible with SMR:

```
SNP    SNP_Chr  SNP_bp   A1  A2   Probe           Probe_Chr  Probe_bp  Gene     Beta       SE         P_value
rs123  1        1191870  T   C    ILMN_2349633    1          1140818   TNFRSF18 -0.436080  0.040022   1.23e-25
...
```

## Unit Testing

Run the test suite:

```bash
python3 -m unittest tests.test_queries -v
```

Run a specific test:

```bash
python3 -m unittest tests.test_queries.TestBESDQueryIndex.test_query_by_probe_id -v
```

### Test Coverage

The test suite includes:

**BESDQueryEngine (original BESD reader):**
- Data loading verification (SNP/probe counts, format detection)
- Single position and range queries
- P-value calculation accuracy
- Beta and SE storage
- Metadata indexing
- Chromosome and position filtering
- Empty query handling

**BESDQueryIndex (SQLite index):**
- Metadata loading from database
- Single position and range queries
- P-value calculation
- Query by probe ID
- Query by SNP ID
- **Consistency verification**: Confirms SQLite index produces identical results to BESD reader

All tests use the `westra_eqtl_hg19` dataset as reference data.

## Performance

The SQLite index provides significant performance advantages:

- **Metadata queries**: Indexed (chr, bp) columns enable O(log n) range queries on millions of SNPs/probes
- **Association storage**: Per-probe BLOBs with numpy arrays maintain sequential read locality
- **Deserialization**: Efficient numpy array deserialization directly from binary format

The index database is typically 50-70% the size of the original BESD files and enables much faster repeated queries.

## File Format

BESDQ expects three files with a common prefix:

- `.besd` - Binary BESD file with association statistics
- `.esi` - SNP index file (chr, rsid, genetic_distance, bp, allele1, allele2, frequency)
- `.epi` - Probe index file (chr, probe_id, genetic_distance, probe_bp, gene, orientation)

## License

MIT
