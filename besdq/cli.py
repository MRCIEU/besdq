"""Command-line interface for BESD query tool."""

import argparse
import sys
from pathlib import Path
from typing import List, Dict

from .besd_reader import BESDQueryEngine
from .builder import BESDIndexBuilder
from .sqlite_query import BESDQueryIndex



def parse_chrpos(chrpos_str: str) -> tuple[str, int, int]:
    """Parse chr:pos or chr:start-end format and return (chr, start_bp, end_bp).

    Supports:
    - Single position: 1:1191870
    - Range: 1:100000-2000000
    """
    if ':' not in chrpos_str:
        raise ValueError(f"Invalid chrpos format: {chrpos_str}. Expected format: chr:pos or chr:start-end")

    chr_part, pos_part = chrpos_str.split(':', 1)

    # Check if it's a range or single position
    if '-' in pos_part:
        # Range format: start-end
        try:
            start_bp, end_bp = map(int, pos_part.split('-'))
        except ValueError:
            raise ValueError(f"Invalid position format in {chrpos_str}. Expected: start-end (integers)")
    else:
        # Single position format: pos
        try:
            pos = int(pos_part)
            start_bp = pos
            end_bp = pos
        except ValueError:
            raise ValueError(f"Invalid position format in {chrpos_str}. Expected: integer position")

    if start_bp > end_bp:
        raise ValueError(f"Invalid range in {chrpos_str}. Start position must be <= end position")

    return chr_part, start_bp, end_bp


def write_output(associations: List[Dict], output_path: str, pval_threshold: float = 0.05):
    """Write results in SMR-compatible format."""
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    # Filter by p-value threshold
    filtered = [a for a in associations if a['pval'] < pval_threshold]

    with open(f"{output_path}.txt", 'w') as f:
        # Header
        f.write('\t'.join([
            'SNP', 'SNP_Chr', 'SNP_bp', 'A1', 'A2',
            'Probe', 'Probe_Chr', 'Probe_bp', 'Gene',
            'Beta', 'SE', 'P_value'
        ]) + '\n')

        # Data
        for assoc in filtered:
            f.write('\t'.join([
                assoc['snp_id'],
                str(assoc['snp_chr']),
                str(assoc['snp_bp']),
                assoc['a1'] or 'NA',
                assoc['a2'] or 'NA',
                assoc['trait_id'],
                str(assoc['trait_chr']) if assoc['trait_chr'] is not None else 'NA',
                str(assoc['trait_bp']) if assoc['trait_bp'] is not None else 'NA',
                assoc['gene'] or 'NA',
                f"{assoc['beta']:.6f}",
                f"{assoc['se']:.6f}",
                f"{assoc['pval']:.6e}",
            ]) + '\n')

    print(f"Results written to {output_path}.txt")
    print(f"  Associations (before filter): {len(associations)}")
    print(f"  Associations (after p-value filter < {pval_threshold}): {len(filtered)}")


def main():
    parser = argparse.ArgumentParser(
        description='Direct BESD query tool (SMR-compatible interface)'
    )

    # Input source: either BESD files or SQLite index
    input_group = parser.add_mutually_exclusive_group(required=True)
    input_group.add_argument('--beqtl-summary',
                        help='Path to BESD files (without extension) - reads from binary files')
    input_group.add_argument('--besd-index',
                        help='Path to SQLite index database file (created with --index) - reads from database')

    parser.add_argument('--index',
                        help='Create SQLite index database at specified path (e.g. data/westra.db). Requires --beqtl-summary.')
    parser.add_argument('--n-mode', choices=['scalar', 'vector'], default='scalar',
                        help='Sample-size storage mode: scalar (one n per probe) or vector (one n per SNP-probe pair). Default: scalar.')
    parser.add_argument('--sample-size', type=int, default=None,
                        help='Scalar sample size N used for --n-mode scalar. If omitted and ESI has allele frequencies, n is computed from the data.')
    parser.add_argument('--trait-variance', default=None, dest='trait_variance',
                        help='Two-column TSV (probe_id  var_y) supplying per-probe trait variance. Probes absent from the file default to var_y=1.0.')
    parser.add_argument('--query', type=float, default=1.0,
                        help='P-value threshold for filtering results (default: 1.0, i.e. no filter)')
    parser.add_argument('--original-scale', action='store_true', dest='original_scale',
                        help='Return beta and SE in original study units instead of SD units (default: SD units)')
    parser.add_argument('--snp-chr',
                        help='SNP chromosome')
    parser.add_argument('--from-snp-kb', type=float,
                        help='SNP region start (kb)')
    parser.add_argument('--to-snp-kb', type=float,
                        help='SNP region end (kb)')
    parser.add_argument('--from-snp-bp', type=int,
                        help='SNP region start (bp)')
    parser.add_argument('--to-snp-bp', type=int,
                        help='SNP region end (bp)')
    parser.add_argument('--snp-chrpos',
                        help='SNP region as chr:pos or chr:start-end (e.g. 1:1191870 or 1:100000-2000000). Comma-separated for multiple ranges.')
    parser.add_argument('--probe-chr',
                        help='Probe chromosome')
    parser.add_argument('--from-probe-kb', type=float,
                        help='Probe region start (kb)')
    parser.add_argument('--to-probe-kb', type=float,
                        help='Probe region end (kb)')
    parser.add_argument('--from-probe-bp', type=int,
                        help='Probe region start (bp)')
    parser.add_argument('--to-probe-bp', type=int,
                        help='Probe region end (bp)')
    parser.add_argument('--probe-chrpos',
                        help='Probe region as chr:pos or chr:start-end (e.g. 1:1140818 or 1:1000000-2000000). Comma-separated for multiple ranges.')
    parser.add_argument('--snp',
                        help='SNP ID(s) to query (comma-separated for multiple, e.g. rs123,rs456)')
    parser.add_argument('--probe',
                        help='Probe ID(s) to query (comma-separated for multiple)')
    parser.add_argument('--gene',
                        help='Gene name(s) to query (comma-separated for multiple)')
    parser.add_argument('--out',
                        help='Output file prefix')

    args = parser.parse_args()

    # Handle indexing mode
    if args.index:
        if not args.beqtl_summary:
            print("Error: --beqtl-summary is required when using --index", file=sys.stderr)
            sys.exit(1)
        try:
            print(f"Building SQLite index from BESD files: {args.beqtl_summary}")
            builder = BESDIndexBuilder(args.index)
            builder.build(
                args.beqtl_summary,
                force=False,
                n_mode=args.n_mode,
                sample_size=args.sample_size,
                trait_variance_path=args.trait_variance,
            )
            print(f"Index created: {args.index}")
            return
        except Exception as e:
            print(f"Error: {e}", file=sys.stderr)
            import traceback
            traceback.print_exc()
            sys.exit(1)

    # Query mode requires --out
    if not args.out:
        print("Error: --out is required for query mode", file=sys.stderr)
        sys.exit(1)

    # Determine data source
    if args.besd_index:
        print(f"Using SQLite index: {args.besd_index}")
        query_engine = BESDQueryIndex(args.besd_index, original_scale=args.original_scale)
        is_index = True
    else:
        print(f"Using BESD files: {args.beqtl_summary}")
        query_engine = BESDQueryEngine(args.beqtl_summary)
        is_index = False

    # Determine query type and validate arguments
    has_snp_query = bool(
        args.snp or args.snp_chrpos or args.snp_chr
        or args.from_snp_kb is not None or args.to_snp_kb is not None
        or args.from_snp_bp is not None or args.to_snp_bp is not None
    )
    has_probe_query = bool(
        args.probe or args.probe_chrpos or args.probe_chr or args.gene
        or args.from_probe_kb is not None or args.to_probe_kb is not None
        or args.from_probe_bp is not None or args.to_probe_bp is not None
    )
    id_query_count = sum(bool(x) for x in (args.snp, args.probe, args.gene))
    has_region_query_args = bool(
        args.snp_chrpos or args.probe_chrpos
        or args.snp_chr or args.probe_chr
        or args.from_snp_kb is not None or args.to_snp_kb is not None
        or args.from_snp_bp is not None or args.to_snp_bp is not None
        or args.from_probe_kb is not None or args.to_probe_kb is not None
        or args.from_probe_bp is not None or args.to_probe_bp is not None
    )

    if id_query_count > 1:
        print("Error: --snp, --probe, and --gene are mutually exclusive", file=sys.stderr)
        sys.exit(1)
    if id_query_count == 1 and has_region_query_args:
        print("Error: Cannot combine --snp/--probe/--gene with region query options", file=sys.stderr)
        sys.exit(1)
    if id_query_count == 0 and (not has_snp_query or not has_probe_query):
        print("Error: Cis-window queries require both SNP and probe query arguments", file=sys.stderr)
        sys.exit(1)

    if args.snp:
        # Query by SNP ID(s)
        snp_ids = [s.strip() for s in args.snp.split(',')]
        print(f"Querying {len(snp_ids)} SNP(s)...")
        associations = []
        for snp_id in snp_ids:
            assocs = query_engine.query_by_snp_id(snp_id)
            associations.extend(assocs)
            print(f"  {snp_id}: {len(assocs)} associations")

    elif args.probe:
        # Query by probe ID(s)
        probe_ids = [p.strip() for p in args.probe.split(',')]
        print(f"Querying {len(probe_ids)} probe(s)...")
        associations = []
        for probe_id in probe_ids:
            assocs = query_engine.query_by_probe_id(probe_id)
            associations.extend(assocs)
            print(f"  {probe_id}: {len(assocs)} associations")

    elif args.gene:
        # Query by gene name(s)
        genes = [g.strip() for g in args.gene.split(',')]
        print(f"Querying {len(genes)} gene(s)...")
        associations = []
        for gene_name in genes:
            assocs = query_engine.query_by_gene(gene_name)
            associations.extend(assocs)
            print(f"  {gene_name}: {len(assocs)} associations")

    else:
        # Cis-window query (original logic)
        # Handle comma-separated chrpos lists
        snp_ranges = []
        probe_ranges = []

        if args.snp_chrpos:
            if args.snp_chr or args.from_snp_kb or args.to_snp_kb or args.from_snp_bp or args.to_snp_bp:
                print("Error: Cannot use --snp-chrpos with --snp-chr, --from-snp-kb, --to-snp-kb, --from-snp-bp, or --to-snp-bp", file=sys.stderr)
                sys.exit(1)
            # Split by comma and parse each range
            try:
                for chrpos in args.snp_chrpos.split(','):
                    chrpos = chrpos.strip()
                    snp_chr, start_bp, end_bp = parse_chrpos(chrpos)
                    snp_ranges.append((snp_chr, start_bp / 1000.0, end_bp / 1000.0))
            except ValueError as e:
                print(f"Error parsing --snp-chrpos: {e}", file=sys.stderr)
                sys.exit(1)

        if args.probe_chrpos:
            if args.probe_chr or args.from_probe_kb or args.to_probe_kb or args.from_probe_bp or args.to_probe_bp:
                print("Error: Cannot use --probe-chrpos with --probe-chr, --from-probe-kb, --to-probe-kb, --from-probe-bp, or --to-probe-bp", file=sys.stderr)
                sys.exit(1)
            # Split by comma and parse each range
            try:
                for chrpos in args.probe_chrpos.split(','):
                    chrpos = chrpos.strip()
                    probe_chr, start_bp, end_bp = parse_chrpos(chrpos)
                    probe_ranges.append((probe_chr, start_bp / 1000.0, end_bp / 1000.0))
            except ValueError as e:
                print(f"Error parsing --probe-chrpos: {e}", file=sys.stderr)
                sys.exit(1)

        # If chrpos was used, combine results from all ranges
        if snp_ranges or probe_ranges:
            if not snp_ranges:
                print("Error: Must specify both --snp-chrpos and --probe-chrpos together", file=sys.stderr)
                sys.exit(1)
            if not probe_ranges:
                print("Error: Must specify both --snp-chrpos and --probe-chrpos together", file=sys.stderr)
                sys.exit(1)

            associations = []
            for snp_chr, snp_start_kb, snp_end_kb in snp_ranges:
                for probe_chr, probe_start_kb, probe_end_kb in probe_ranges:
                    print(f"Querying SNP {snp_chr}:{snp_start_kb*1000:.0f}-{snp_end_kb*1000:.0f} bp, Probe {probe_chr}:{probe_start_kb*1000:.0f}-{probe_end_kb*1000:.0f} bp")
                    range_assocs = query_engine.query_cis_window(
                        snp_chr=snp_chr,
                        snp_start_kb=snp_start_kb,
                        snp_end_kb=snp_end_kb,
                        probe_chr=probe_chr,
                        probe_start_kb=probe_start_kb,
                        probe_end_kb=probe_end_kb,
                    )
                    associations.extend(range_assocs)
        else:
            # Standard range query (not using chrpos)
            snp_chr = args.snp_chr
            snp_start_kb = None
            snp_end_kb = None
            probe_chr = args.probe_chr
            probe_start_kb = None
            probe_end_kb = None

            if args.from_snp_kb is not None and args.from_snp_bp is not None:
                print("Error: Cannot specify both --from-snp-kb and --from-snp-bp", file=sys.stderr)
                sys.exit(1)
            if args.to_snp_kb is not None and args.to_snp_bp is not None:
                print("Error: Cannot specify both --to-snp-kb and --to-snp-bp", file=sys.stderr)
                sys.exit(1)

            if args.from_snp_bp is not None:
                snp_start_kb = args.from_snp_bp / 1000.0
            elif args.from_snp_kb is not None:
                snp_start_kb = args.from_snp_kb
            else:
                print("Error: Must specify either --snp-chrpos, --from-snp-kb, or --from-snp-bp", file=sys.stderr)
                sys.exit(1)

            if args.to_snp_bp is not None:
                snp_end_kb = args.to_snp_bp / 1000.0
            elif args.to_snp_kb is not None:
                snp_end_kb = args.to_snp_kb
            else:
                print("Error: Must specify either --snp-chrpos, --to-snp-kb, or --to-snp-bp", file=sys.stderr)
                sys.exit(1)

            # Validate and convert probe position arguments
            if args.from_probe_kb is not None and args.from_probe_bp is not None:
                print("Error: Cannot specify both --from-probe-kb and --from-probe-bp", file=sys.stderr)
                sys.exit(1)
            if args.to_probe_kb is not None and args.to_probe_bp is not None:
                print("Error: Cannot specify both --to-probe-kb and --to-probe-bp", file=sys.stderr)
                sys.exit(1)

            if args.from_probe_bp is not None:
                probe_start_kb = args.from_probe_bp / 1000.0
            elif args.from_probe_kb is not None:
                probe_start_kb = args.from_probe_kb
            else:
                print("Error: Must specify either --probe-chrpos, --from-probe-kb, or --from-probe-bp", file=sys.stderr)
                sys.exit(1)

            if args.to_probe_bp is not None:
                probe_end_kb = args.to_probe_bp / 1000.0
            elif args.to_probe_kb is not None:
                probe_end_kb = args.to_probe_kb
            else:
                print("Error: Must specify either --probe-chrpos, --to-probe-kb, or --to-probe-bp", file=sys.stderr)
                sys.exit(1)

            print(f"Querying cis-window:")
            print(f"  SNP:   {snp_chr}:{snp_start_kb*1000:.0f}-{snp_end_kb*1000:.0f} bp")
            print(f"  Probe: {probe_chr}:{probe_start_kb*1000:.0f}-{probe_end_kb*1000:.0f} bp")

            associations = query_engine.query_cis_window(
                snp_chr=snp_chr,
                snp_start_kb=snp_start_kb,
                snp_end_kb=snp_end_kb,
                probe_chr=probe_chr,
                probe_start_kb=probe_start_kb,
                probe_end_kb=probe_end_kb,
            )

    try:
        print(f"\nFound {len(associations)} associations")

        if associations:
            write_output(associations, args.out, pval_threshold=args.query)

    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        import traceback
        traceback.print_exc()
        sys.exit(1)
    finally:
        # Clean up database connection if using index
        if is_index:
            query_engine.close()



def import_gwas_ssf_main() -> None:
    """Entry point for import-gwas-ssf command."""
    parser = argparse.ArgumentParser(
        prog='import-gwas-ssf',
        description='Import GWAS-SSF summary statistics into a queryable SQLite index',
    )
    parser.add_argument(
        '--trait-annotation', required=True, dest='trait_annotation',
        help='Path to trait annotation TSV file (required columns: file_path, trait_id, trait_name)',
    )
    parser.add_argument(
        '--ld-reference', required=True, dest='ld_reference',
        help='Prefix for plink2-format LD reference panel (--pfile prefix)',
    )
    parser.add_argument(
        '--output', default=None,
        help='Output SQLite database path (default: <annotation_tsv_stem>.db)',
    )
    parser.add_argument(
        '--workers', type=int, default=1,
        help='Number of parallel workers for Pass 1 (default: 1)',
    )
    parser.add_argument(
        '--sig-threshold', type=float, default=5e-8, dest='sig_threshold',
        help='Genome-wide significance threshold (default: 5e-8)',
    )
    parser.add_argument(
        '--sug-threshold', type=float, default=1e-4, dest='sug_threshold',
        help='Suggestive significance threshold (default: 1e-4)',
    )
    parser.add_argument(
        '--cis-radius', type=int, default=1_000_000, dest='cis_radius',
        help='Cis window radius in bp (default: 1000000)',
    )
    parser.add_argument(
        '--sig-radius', type=int, default=500_000, dest='sig_radius',
        help='Significant trans window radius in bp (default: 500000)',
    )
    parser.add_argument(
        '--clump-r2', type=float, default=0.01, dest='clump_r2',
        help='LD clumping r² threshold (default: 0.01)',
    )
    parser.add_argument(
        '--clump-kb', type=int, default=10_000, dest='clump_kb',
        help='LD clumping window in kb (default: 10000)',
    )
    parser.add_argument(
        '--force-fallback', action='store_true', dest='force_fallback',
        help='Always use Python streaming reader instead of the pytabix+awk fast path',
    )

    args = parser.parse_args()

    # Route besdq.gwas_ssf_fast_reader INFO messages to stderr so the user can
    # see which read path is selected for each file.
    import logging as _logging
    _logging.basicConfig(
        level=_logging.INFO,
        format='[%(asctime)s] %(message)s',
        datefmt='%H:%M:%S',
        stream=sys.stderr,
    )

    # Validate annotation file exists before starting processing
    annotation_path = Path(args.trait_annotation)
    if not annotation_path.exists():
        print(f"Error: annotation file not found: {args.trait_annotation}", file=sys.stderr)
        sys.exit(1)

    # Determine output path
    output_path = args.output or str(annotation_path.stem + '.db')

    try:
        from .annotation_reader import read_trait_annotation
        from .gwas_ssf_builder import GwasSsfIndexBuilder

        print(f"Reading trait annotations from: {args.trait_annotation}", file=sys.stderr)
        traits = read_trait_annotation(args.trait_annotation)
        print(f"  {len(traits)} traits found", file=sys.stderr)

        builder = GwasSsfIndexBuilder(output_path)
        builder.build(
            traits,
            workers=args.workers,
            cis_radius=args.cis_radius,
            sig_threshold=args.sig_threshold,
            sug_threshold=args.sug_threshold,
            plink2_pfile=args.ld_reference,
            sig_radius=args.sig_radius,
            clump_r2=args.clump_r2,
            clump_kb=args.clump_kb,
            force_fallback=args.force_fallback,
        )
        print(f"Index written to: {output_path}", file=sys.stderr)

    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        import traceback
        traceback.print_exc()
        sys.exit(1)


if __name__ == '__main__':
    main()
