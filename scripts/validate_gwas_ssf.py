#!/usr/bin/env python3
"""Validate that betas stored in a BESDQ index match the original GWAS-SSF source files.

Queries the index in original-scale units (bypassing SD normalisation) and compares
directly to normalised betas from the source files. Produces a scatter plot and
prints summary statistics.

Usage:
    python scripts/validate_gwas_ssf.py \\
        --index data/ebi_input/study.db \\
        --annotation data/ebi_input/traits.tsv \\
        --out validation

Output:
    <out>_betas.png   — scatter: original beta vs index beta, per trait
    <out>_report.txt  — correlation and max absolute difference per trait
"""

import argparse
import sys
import warnings
from pathlib import Path

import numpy as np


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Validate BESDQ index betas against original GWAS-SSF files"
    )
    parser.add_argument("--index", required=True, help="Path to BESDQ SQLite index")
    parser.add_argument(
        "--annotation", required=True, help="Trait annotation TSV used to build the index"
    )
    parser.add_argument("--out", default="validation", help="Output file prefix (default: validation)")
    args = parser.parse_args()

    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print("ERROR: matplotlib is required. Install with: pip install matplotlib", file=sys.stderr)
        sys.exit(1)

    # Add repo root to path so besdq is importable when run from any directory
    repo_root = Path(__file__).resolve().parent.parent
    if str(repo_root) not in sys.path:
        sys.path.insert(0, str(repo_root))

    from besdq.annotation_reader import read_trait_annotation
    from besdq.gwas_ssf_reader import read_gwas_ssf
    from besdq.sqlite_query import BESDQueryIndex

    annotation_path = Path(args.annotation)
    if not annotation_path.exists():
        print(f"ERROR: annotation file not found: {args.annotation}", file=sys.stderr)
        sys.exit(1)

    index_path = Path(args.index)
    if not index_path.exists():
        print(f"ERROR: index not found: {args.index}", file=sys.stderr)
        sys.exit(1)

    traits = read_trait_annotation(str(annotation_path))
    print(f"Loaded {len(traits)} traits from annotation TSV")

    # Collect per-trait data
    trait_data = {}  # trait_id -> {'orig': [], 'idx': []}

    with BESDQueryIndex(str(index_path), original_scale=True) as idx:
        for trait in traits:
            tid = trait.trait_id
            print(f"  Querying {tid}...", end=" ", flush=True)

            assocs = idx.query_by_probe_id(tid)
            if not assocs:
                print("0 associations in index, skipping")
                continue

            # Build lookup: snp_key -> beta from index
            idx_lookup: dict = {}
            for a in assocs:
                key = f"{a['snp_chr']}:{a['snp_bp']}:{a['a1']}:{a['a2']}"
                idx_lookup[key] = a["beta"]

            # Stream original file and match
            orig_betas = []
            idx_betas = []
            file_path = trait.file_path
            # Resolve relative to annotation file directory if needed
            if not Path(file_path).is_absolute():
                file_path = str(annotation_path.parent.parent / file_path)

            for row in read_gwas_ssf(file_path):
                if row.snp_key in idx_lookup:
                    orig_betas.append(row.beta)
                    idx_betas.append(idx_lookup[row.snp_key])

            n_matched = len(orig_betas)
            n_index = len(assocs)
            print(f"{n_index} in index, {n_matched} matched in source")

            if n_matched > 0:
                trait_data[tid] = {
                    "orig": np.array(orig_betas),
                    "idx": np.array(idx_betas),
                    "gene": trait.gene or tid,
                }

    if not trait_data:
        print("No matched associations found. Check that file paths in annotation TSV are correct.")
        sys.exit(1)

    # --- Plot ---
    n_traits = len(trait_data)
    ncols = min(n_traits, 3)
    nrows = (n_traits + ncols - 1) // ncols
    fig, axes = plt.subplots(nrows, ncols, figsize=(5 * ncols, 5 * nrows), squeeze=False)

    report_lines = ["trait_id\tn_matched\tcorrelation\tmax_abs_diff\tmean_abs_diff"]

    for ax_idx, (tid, d) in enumerate(trait_data.items()):
        row_i, col_i = divmod(ax_idx, ncols)
        ax = axes[row_i][col_i]

        orig = d["orig"]
        idx_b = d["idx"]
        label = d["gene"]

        ax.scatter(orig, idx_b, alpha=0.4, s=4, color="steelblue", rasterized=True)

        # y = x reference line
        lim = np.abs(np.concatenate([orig, idx_b])).max() * 1.05
        ax.plot([-lim, lim], [-lim, lim], "r-", linewidth=1, label="y = x")
        ax.set_xlim(-lim, lim)
        ax.set_ylim(-lim, lim)
        ax.set_xlabel("Original beta (GWAS-SSF)")
        ax.set_ylabel("Index beta (original scale)")
        ax.set_title(f"{label}\n(n={len(orig):,})")
        ax.legend(fontsize=8)
        ax.set_aspect("equal")

        # Stats
        corr = float(np.corrcoef(orig, idx_b)[0, 1]) if len(orig) > 1 else float("nan")
        residuals = idx_b - orig
        max_diff = float(np.abs(residuals).max())
        mean_diff = float(np.abs(residuals).mean())

        ax.text(
            0.05, 0.95,
            f"r = {corr:.6f}\nmax |Δ| = {max_diff:.2e}",
            transform=ax.transAxes,
            verticalalignment="top",
            fontsize=8,
            bbox=dict(boxstyle="round", facecolor="wheat", alpha=0.5),
        )

        report_lines.append(f"{tid}\t{len(orig)}\t{corr:.8f}\t{max_diff:.6e}\t{mean_diff:.6e}")

    # Hide unused axes
    for ax_idx in range(n_traits, nrows * ncols):
        row_i, col_i = divmod(ax_idx, ncols)
        axes[row_i][col_i].set_visible(False)

    fig.suptitle("BESDQ index validation: stored vs original betas", fontsize=13)
    plt.tight_layout()

    plot_path = f"{args.out}_betas.png"
    plt.savefig(plot_path, dpi=150)
    print(f"\nPlot saved: {plot_path}")

    report_path = f"{args.out}_report.txt"
    with open(report_path, "w") as f:
        f.write("\n".join(report_lines) + "\n")
    print(f"Report saved: {report_path}")

    # Print summary to stdout
    print("\n--- Validation summary ---")
    for line in report_lines:
        print(line)


if __name__ == "__main__":
    main()
