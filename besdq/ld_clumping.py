"""LD clumping via plink2 to identify independent significant trans peaks."""

import shutil
import subprocess
import tempfile
from pathlib import Path
from typing import List, Optional

from .gwas_ssf_reader import GwasSsfRow


def _check_plink2() -> str:
    """Return path to plink2 or raise ImportError if not found."""
    path = shutil.which('plink2')
    if path is None:
        raise ImportError(
            "plink2 not found on PATH. Install with: conda install -c bioconda plink2  "
            "or: mamba install -c bioconda plink2"
        )
    return path


def clump_trans_peaks(
    candidates: List[GwasSsfRow],
    plink2_pfile: Optional[str],
    sig_radius: int = 500_000,
    clump_r2: float = 0.01,
    clump_kb: int = 10_000,
    all_rows: Optional[List[GwasSsfRow]] = None,
) -> List[GwasSsfRow]:
    """Run LD clumping on significant trans candidates and expand to windows.

    Parameters
    ----------
    candidates : genome-wide significant associations to clump
    plink2_pfile : prefix for plink2 --pfile LD reference (required unless candidates is empty)
    sig_radius : bp half-window around each lead SNP to retain
    clump_r2 : r² threshold for LD clumping
    clump_kb : kb window for LD clumping
    all_rows : full set of rows to expand windows from (if None, expands from candidates only)

    Returns
    -------
    List of GwasSsfRow within sig_radius of each independent lead SNP
    """
    if not candidates:
        return []

    plink2 = _check_plink2()

    if plink2_pfile is None:
        raise ValueError("plink2_pfile is required for LD clumping")

    pfile_path = Path(plink2_pfile)
    # Check at least one of the plink2 format files exists
    if not (pfile_path.with_suffix('.pgen').exists() or pfile_path.with_suffix('.bed').exists()):
        raise FileNotFoundError(f"plink2 reference files not found at: {plink2_pfile}")

    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir = Path(tmpdir)

        # Write candidates as plink2-compatible GWAS summary stats
        ssf_file = tmpdir / 'candidates.ssf'
        with open(ssf_file, 'w') as fh:
            fh.write('ID\tP\n')
            for row in candidates:
                fh.write(f"{row.snp_key}\t{row.p}\n")

        out_prefix = str(tmpdir / 'clump')
        cmd = [
            plink2,
            '--pfile', str(plink2_pfile),
            '--clump', str(ssf_file),
            '--clump-p1', str(5e-8),
            '--clump-r2', str(clump_r2),
            '--clump-kb', str(clump_kb),
            '--out', out_prefix,
        ]

        result = subprocess.run(cmd, capture_output=True, text=True)
        if result.returncode != 0:
            raise RuntimeError(
                f"plink2 clumping failed (exit {result.returncode}):\n{result.stderr}"
            )

        # Parse lead SNP positions from clump output
        clump_out = Path(out_prefix + '.clumps')
        lead_positions: list[tuple[str, int]] = []

        if clump_out.exists():
            with open(clump_out) as fh:
                header = fh.readline().strip().split()
                chr_col = header.index('#CHROM') if '#CHROM' in header else 0
                pos_col = header.index('POS') if 'POS' in header else 1
                for line in fh:
                    parts = line.strip().split()
                    if len(parts) > max(chr_col, pos_col):
                        try:
                            lead_positions.append((parts[chr_col], int(parts[pos_col])))
                        except (ValueError, IndexError):
                            continue

    if not lead_positions:
        # No clumped leads found; fall back to returning all candidates
        return candidates

    # Expand each lead SNP to a ±sig_radius window
    source = all_rows if all_rows is not None else candidates
    retained = []
    seen = set()

    for row in source:
        for lead_chr, lead_bp in lead_positions:
            if row.chr == lead_chr and abs(row.bp - lead_bp) <= sig_radius:
                key = row.snp_key
                if key not in seen:
                    seen.add(key)
                    retained.append(row)
                break

    return retained
