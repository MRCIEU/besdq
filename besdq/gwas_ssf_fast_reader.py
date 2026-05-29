"""Fast GWAS-SSF reader: pytabix cis lookup + optimised Python gzip scan.

Two paths, both use the same partial-split gzip scan for the p-value filter:

  Fast path (pytabix available + .tbi present):
    1. pytabix  — random-access cis window (milliseconds, no full-file scan).
    2. Partial-split gzip scan — trans rows with p < threshold, skipping the
       cis window (already collected in step 1).

  Fallback (always works, no extra dependencies):
    Partial-split gzip scan — cis rows unconditionally + trans rows with
    p < threshold.

The "partial-split" trick: for every line, only split up to the p_value column
(and the chromosome/bp columns needed for the cis check).  The full tab-split
is deferred to the ~1 % of rows that actually pass — cutting gzip scan time
by ~3× vs a naive full split on every line.

Install pytabix for the fast path:  pip install besdq[fast]
"""

import gzip
import logging
import os
from pathlib import Path
import shlex
import subprocess
from typing import Generator, List, Optional

from .gwas_ssf_reader import GwasSsfRow, _parse_optional_rsid

logger = logging.getLogger(__name__)


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
    """Return the 1-based awk column index of the p_value column (kept for API compat)."""
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
# Row parsing
# ---------------------------------------------------------------------------

def _make_col_idx(cols: List[str]) -> dict:
    return {name: i for i, name in enumerate(cols)}


def _parse_row(parts: List[str], col_idx: dict) -> Optional[GwasSsfRow]:
    """Parse a pre-split TSV row into a GwasSsfRow with allele normalisation.

    Accepts both plain text line parts and pytabix field lists.
    Strips whitespace from each field to handle Windows line endings.
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
# Core: partial-split gzip scan
# ---------------------------------------------------------------------------

def _gzip_scan(
    path: str,
    cis_chr: Optional[str],
    cis_start: Optional[int],
    cis_end: Optional[int],
    p_threshold: float,
    skip_cis: bool = False,
) -> Generator[GwasSsfRow, None, None]:
    """Stream the file with a fast partial-split early exit for non-candidate rows.

    For each line only the columns needed for the cis-check and p-value check
    are split; the remaining columns are split only for the ~1% of rows that
    pass.  This cuts scan time by ~3× versus a naive full split on every line.

    Args:
        skip_cis: When True, cis rows are identified but NOT yielded (they will
                  be provided by a pytabix lookup instead).  Non-cis rows that
                  pass p_threshold are still yielded.
    """
    has_cis = cis_chr is not None and cis_start is not None and cis_end is not None

    with gzip.open(path, 'rt') as fh:
        cols = fh.readline().rstrip('\n').split('\t')
        col_idx = _make_col_idx(cols)
        chr_col = col_idx['chromosome']
        bp_col  = col_idx['base_pair_location']
        p_col   = col_idx['p_value']
        # Split just enough columns to check cis membership and p-value
        split_n = max(chr_col, bp_col, p_col) + 1

        for line in fh:
            quick = line.split('\t', split_n)
            if len(quick) <= p_col:
                continue

            try:
                p = float(quick[p_col])
            except ValueError:
                continue

            in_cis = False
            if has_cis:
                try:
                    bp = int(quick[bp_col])
                    in_cis = quick[chr_col] == cis_chr and cis_start <= bp <= cis_end
                except (ValueError, IndexError):
                    pass

            if in_cis:
                if skip_cis:
                    continue  # pytabix already collected this row
                # fallback path: yield cis rows unconditionally
            elif p >= p_threshold:
                continue  # trans row that doesn't pass threshold

            # Full parse only for rows we're going to yield
            parts = line.rstrip('\n').split('\t')
            row = _parse_row(parts, col_idx)
            if row is not None:
                yield row


# ---------------------------------------------------------------------------
# Fallback path (Python only)
# ---------------------------------------------------------------------------

def _fallback_candidates(
    path: str,
    cis_chr: Optional[str],
    cis_start: Optional[int],
    cis_end: Optional[int],
    p_threshold: float,
) -> Generator[GwasSsfRow, None, None]:
    """Partial-split gzip scan: cis rows unconditionally + trans rows below threshold."""
    yield from _gzip_scan(path, cis_chr, cis_start, cis_end, p_threshold, skip_cis=False)


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
# Fast path: pytabix cis lookup + partial-split gzip scan for trans
# ---------------------------------------------------------------------------

def _pytabix_cis_rows(
    path: str,
    cis_chr: str,
    cis_start: int,
    cis_end: int,
    col_idx: dict,
) -> Generator[GwasSsfRow, None, None]:
    """Use pytabix to fetch all rows in the cis window.

    pytabix returns each record as a pre-split list of field strings.
    The last field may carry a Windows '\\r'; _parse_row strips it.
    """
    import tabix  # lazy import — module loads fine without pytabix installed
    tb = tabix.open(path)
    try:
        for parts in tb.query(cis_chr, cis_start, cis_end):
            row = _parse_row(parts, col_idx)
            if row is not None:
                yield row
    except tabix.TabixError:
        return  # region not in index — no cis rows


def _fast_path_candidates(
    path: str,
    cis_chr: Optional[str],
    cis_start: Optional[int],
    cis_end: Optional[int],
    p_threshold: float,
) -> Generator[GwasSsfRow, None, None]:
    """Fast path: pytabix for instant cis lookup + partial-split gzip scan for trans."""
    has_cis = cis_chr is not None and cis_start is not None and cis_end is not None
    seen: set = set()

    # 1. pytabix cis rows (random-access, no full-file scan needed)
    if has_cis:
        col_idx = _make_col_idx(_get_header(path))
        for row in _pytabix_cis_rows(path, cis_chr, cis_start, cis_end, col_idx):
            key = row.snp_key
            if key not in seen:
                seen.add(key)
                yield row

    # 2. Partial-split gzip scan for trans rows (skip cis — already collected above)
    for row in _gzip_scan(path, cis_chr, cis_start, cis_end, p_threshold, skip_cis=has_cis):
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
    force_fallback: bool = False,
) -> Generator[GwasSsfRow, None, None]:
    """Yield candidate rows from a GWAS-SSF file.

    Candidate rows are:
      - All rows inside the cis window (cis_chr:cis_start-cis_end), unconditionally.
      - All rows with p_value < p_threshold, regardless of location.

    Both paths use a partial-split optimisation that only fully parses lines
    that pass the threshold — approximately 3× faster than a naive scan.

    Uses the pytabix fast path when all of the following hold:
      - force_fallback is False
      - the pytabix package is installed  (pip install besdq[fast])
      - a .tbi index exists alongside the file

    The fast path additionally uses pytabix for instant cis-window random access
    (useful for cis-only queries), then sweeps the file for trans rows.

    Args:
        path:           BGZF-compressed (.tsv.gz) GWAS-SSF file.
        cis_chr:        Chromosome of the cis window (None → no cis window).
        cis_start:      Start bp of the cis window (inclusive).
        cis_end:        End bp of the cis window (inclusive).
        p_threshold:    P-value threshold for trans candidates.
        force_fallback: Always use Python streaming, ignore pytabix.

    Yields:
        GwasSsfRow instances with alleles normalised (a1 <= a2 alphabetically).
    """
    filename = Path(path).name
    tbi_path = path + '.tbi'

    if force_fallback:
        logger.info("%s: Python streaming (forced)", filename)
        yield from _fallback_candidates(path, cis_chr, cis_start, cis_end, p_threshold)
    elif not _pytabix_available():
        logger.info("%s: Python streaming (pytabix not installed)", filename)
        yield from _fallback_candidates(path, cis_chr, cis_start, cis_end, p_threshold)
    elif not os.path.exists(tbi_path):
        logger.info("%s: Python streaming (no .tbi index found)", filename)
        yield from _fallback_candidates(path, cis_chr, cis_start, cis_end, p_threshold)
    else:
        logger.info("%s: pytabix fast path", filename)
        yield from _fast_path_candidates(path, cis_chr, cis_start, cis_end, p_threshold)
