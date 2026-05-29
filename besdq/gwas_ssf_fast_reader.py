"""Fast GWAS-SSF reader: pytabix + awk fast path with Python streaming fallback.

Fast path (requires the `pytabix` package and a .tbi index alongside the BGZF file):
  1. pytabix: fetch the cis window unconditionally (random-access, no CLI dependency).
  2. awk:     stream the full file and keep rows where p_value < p_threshold.
  3. Merge and deduplicate.

Fallback (always safe, no extra dependencies):
  Stream the full file via gzip.open, apply the same cis/trans logic in Python.

Use read_gwas_ssf_candidates() — it picks the fast path automatically when
pytabix is importable and a .tbi file exists next to the data file.
"""

import gzip
import os
import shlex
import subprocess
from typing import Generator, List, Optional

from .gwas_ssf_reader import GwasSsfRow, _parse_optional_rsid


# ---------------------------------------------------------------------------
# Utilities
# ---------------------------------------------------------------------------

def _pytabix_available() -> bool:
    """Return True if the pytabix package is importable."""
    try:
        import tabix  # noqa: F401
        return True
    except ImportError:
        return False


def _get_header(path: str) -> List[str]:
    """Return the list of column names from the header of a gzipped TSV."""
    with gzip.open(path, 'rt') as fh:
        return fh.readline().rstrip('\n').split('\t')


def _get_pval_col(path: str) -> int:
    """Return the 1-based awk column index of the p_value column."""
    return _get_header(path).index('p_value') + 1


def _count_data_lines(path: str) -> int:
    """Count data rows (excluding the header) in a gzipped file via subprocess."""
    result = subprocess.run(
        f'gzip -dc {shlex.quote(path)} | wc -l',
        shell=True,
        capture_output=True,
        text=True,
        check=True,
    )
    total = int(result.stdout.strip())
    return total - 1  # subtract header line


# ---------------------------------------------------------------------------
# Row parsing (shared between fallback and fast-path output)
# ---------------------------------------------------------------------------

def _make_col_idx(cols: List[str]) -> dict:
    return {name: i for i, name in enumerate(cols)}


def _parse_row(parts: List[str], col_idx: dict) -> Optional[GwasSsfRow]:
    """Parse a split TSV row into a GwasSsfRow with allele normalisation.

    Accepts both plain text lines (split on '\\t') and pytabix field lists.
    Returns None for rows that cannot be parsed.
    """
    chr_col = col_idx['chromosome']
    bp_col = col_idx['base_pair_location']
    ea_col = col_idx['effect_allele']
    oa_col = col_idx['other_allele']
    beta_col = col_idx['beta']
    se_col = col_idx['standard_error']
    eaf_col = col_idx['effect_allele_frequency']
    p_col = col_idx['p_value']
    rsid_col = col_idx.get('rsid', col_idx.get('rs_id'))

    required = max(chr_col, bp_col, ea_col, oa_col, beta_col, se_col, eaf_col, p_col)
    if len(parts) <= required:
        return None

    try:
        effect_allele = parts[ea_col].strip()
        other_allele = parts[oa_col].strip()
        beta = float(parts[beta_col])
        se = float(parts[se_col])
        eaf = float(parts[eaf_col])
        p = float(parts[p_col])
        chromosome = str(parts[chr_col]).strip()
        bp = int(parts[bp_col])
    except (ValueError, IndexError):
        return None

    rsid = _parse_optional_rsid(parts[rsid_col].strip()) if rsid_col is not None else None

    # Normalise: a1 is always alphabetically first
    if effect_allele <= other_allele:
        a1, a2 = effect_allele, other_allele
    else:
        a1, a2 = other_allele, effect_allele
        beta = -beta
        eaf = 1.0 - eaf

    return GwasSsfRow(
        chr=chromosome,
        bp=bp,
        a1=a1,
        a2=a2,
        rsid=rsid,
        beta=beta,
        se=se,
        eaf=eaf,
        p=p,
    )


# ---------------------------------------------------------------------------
# Python fallback path
# ---------------------------------------------------------------------------

def _fallback_candidates(
    path: str,
    cis_chr: Optional[str],
    cis_start: Optional[int],
    cis_end: Optional[int],
    p_threshold: float,
) -> Generator[GwasSsfRow, None, None]:
    """Stream the whole file; yield cis rows unconditionally + trans rows below threshold."""
    has_cis = cis_chr is not None and cis_start is not None and cis_end is not None

    with gzip.open(path, 'rt') as fh:
        cols = fh.readline().rstrip('\n').split('\t')
        col_idx = _make_col_idx(cols)

        for line in fh:
            parts = line.rstrip('\n').split('\t')
            row = _parse_row(parts, col_idx)
            if row is None:
                continue

            if has_cis and row.chr == cis_chr and cis_start <= row.bp <= cis_end:
                yield row
            elif row.p < p_threshold:
                yield row


def _force_fallback_candidates(
    path: str,
    cis_chr: Optional[str] = None,
    cis_start: Optional[int] = None,
    cis_end: Optional[int] = None,
    p_threshold: float = 1e-4,
) -> Generator[GwasSsfRow, None, None]:
    """Always use the Python fallback (for testing / comparison against the fast path)."""
    yield from _fallback_candidates(path, cis_chr, cis_start, cis_end, p_threshold)


# ---------------------------------------------------------------------------
# Fast path: pytabix (cis lookup) + awk (trans filter)
# ---------------------------------------------------------------------------

def _pytabix_cis_rows(
    path: str,
    cis_chr: str,
    cis_start: int,
    cis_end: int,
    col_idx: dict,
) -> Generator[GwasSsfRow, None, None]:
    """Use pytabix to fetch all rows in the cis window; parse and yield GwasSsfRow.

    pytabix returns each record as a list of field strings (no trailing newline).
    The last field may carry a '\\r' from Windows line endings — _parse_row
    strips whitespace from each field before conversion.
    """
    import tabix  # imported lazily so the module loads without pytabix installed
    tb = tabix.open(path)
    try:
        for parts in tb.query(cis_chr, cis_start, cis_end):
            row = _parse_row(parts, col_idx)
            if row is not None:
                yield row
    except tabix.TabixError:
        # Region not present in index — no cis rows
        return


def _awk_trans_candidates(
    path: str,
    p_threshold: float,
    pval_col: int,
    col_idx: dict,
) -> Generator[GwasSsfRow, None, None]:
    """Use awk to filter rows where p_value < threshold; parse and yield GwasSsfRow."""
    awk_prog = f'NR>1 && ${pval_col}+0 < {p_threshold}'
    cmd = f'gzip -dc {shlex.quote(path)} | awk -F\'\\t\' \'{awk_prog}\''
    proc = subprocess.Popen(
        cmd, shell=True, stdout=subprocess.PIPE, text=True,
    )
    try:
        for line in proc.stdout:
            parts = line.rstrip('\n').split('\t')
            row = _parse_row(parts, col_idx)
            if row is not None:
                yield row
    finally:
        proc.stdout.close()
        proc.wait()


def _fast_path_candidates(
    path: str,
    cis_chr: Optional[str],
    cis_start: Optional[int],
    cis_end: Optional[int],
    p_threshold: float,
) -> Generator[GwasSsfRow, None, None]:
    """Fast path: pytabix for cis + awk for trans, merged and deduplicated."""
    pval_col = _get_pval_col(path)
    cols = _get_header(path)
    col_idx = _make_col_idx(cols)

    seen: set = set()

    # pytabix cis rows first (unconditional random-access lookup)
    has_cis = cis_chr is not None and cis_start is not None and cis_end is not None
    if has_cis:
        for row in _pytabix_cis_rows(path, cis_chr, cis_start, cis_end, col_idx):
            key = row.snp_key
            if key not in seen:
                seen.add(key)
                yield row

    # awk trans rows (p < threshold; cis rows that also pass are dedup'd away)
    for row in _awk_trans_candidates(path, p_threshold, pval_col, col_idx):
        key = row.snp_key
        if key not in seen:
            seen.add(key)
            yield row


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def read_gwas_ssf_candidates(
    path: str,
    cis_chr: Optional[str] = None,
    cis_start: Optional[int] = None,
    cis_end: Optional[int] = None,
    p_threshold: float = 1e-4,
) -> Generator[GwasSsfRow, None, None]:
    """Yield candidate rows from a GWAS-SSF file.

    Candidate rows are:
      - All rows inside the cis window (cis_chr:cis_start-cis_end), unconditionally.
      - All rows with p_value < p_threshold, regardless of location.

    Uses the pytabix + awk fast path when:
      - the `pytabix` package is installed, and
      - a .tbi index exists alongside the file.

    Otherwise falls back to Python streaming (always correct, no extra dependencies).

    Args:
        path:        Path to a BGZF-compressed (.tsv.gz) GWAS-SSF file.
        cis_chr:     Chromosome of the cis window (or None for no cis window).
        cis_start:   Start bp of the cis window (inclusive).
        cis_end:     End bp of the cis window (inclusive).
        p_threshold: P-value threshold for trans candidates.

    Yields:
        GwasSsfRow instances with alleles normalised (a1 <= a2 alphabetically).
    """
    tbi_path = path + '.tbi'
    if _pytabix_available() and os.path.exists(tbi_path):
        yield from _fast_path_candidates(path, cis_chr, cis_start, cis_end, p_threshold)
    else:
        yield from _fallback_candidates(path, cis_chr, cis_start, cis_end, p_threshold)
