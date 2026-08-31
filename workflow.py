#!/usr/bin/env python3
"""
IMGT antibody (IG) germline gene cross-species comparison workflow.
this tool only covers IGH/IGK/IGL.

USAGE (structured):
    python3 workflow.py --species_a "macaque" --species_b "mouse" --chain heavy --gene V

USAGE (freeform):
    python3 workflow.py --query "macaque vs mouse, heavy chain, V genes"

Chain options:  heavy, kappa, lambda, light (light = runs kappa AND lambda)
Gene options:   V, D, J, C   (C = constant region; light-chain C-REGION is used
                directly, heavy-chain C is approximated by the CH1 exon since
                IGH constant genes span multiple exons)

Species names accept common names (human, mouse, macaque, rabbit, ...) or
exact/partial IMGT scientific names, with fuzzy-match fallback.
Run `python3 list_species.py` to see everything available in your bulk file.
"""
import argparse
import sys
from imgt_parser import parse_bulk_fasta
from query_resolver import resolve_query, parse_freeform_query, QueryError
from run_comparison import compare_two_species

DEFAULT_BULK = '/mnt/user-data/uploads/IMGTGENEDB-ReferenceSequences.fasta-AA-WithGaps-F_ORF_inframeP'
DEFAULT_OUT = '/mnt/user-data/outputs'

def build_argparser():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--query', help='Freeform query, e.g. "macaque vs mouse, heavy chain, V genes"')
    p.add_argument('--species_a', help='First species (common or scientific name)')
    p.add_argument('--species_b', help='Second species (common or scientific name)')
    p.add_argument('--chain', default='heavy', help='heavy | kappa | lambda | light (default: heavy)')
    p.add_argument('--gene', default='V', help='V | D | J | C (default: V)')
    p.add_argument('--bulk_path', default=DEFAULT_BULK, help='Path to IMGT/GENE-DB bulk FASTA file')
    p.add_argument('--out_dir', default=DEFAULT_OUT, help='Output directory')
    return p

def run(species_a_raw, species_b_raw, chain_raw, gene_raw, bulk_path, out_dir):
    print(f'Loading {bulk_path} ...')
    records = parse_bulk_fasta(bulk_path)  # ig_only=True by default -- TCR excluded
    print(f'Loaded {len(records)} antibody (IG) gene records (V/D/J/C, all functionality classes).')

    available_species = sorted(set(r['species'] for r in records))

    try:
        resolved = resolve_query(species_a_raw, species_b_raw, chain_raw, gene_raw, available_species)
    except QueryError as e:
        print(f'\nERROR: {e}')
        print('\nAvailable species in this file:')
        for s in available_species:
            print(f'  - {s}')
        sys.exit(1)

    print(f"\nResolved species A: '{species_a_raw}' -> {resolved['species_a']} ({resolved['species_a_match']} match)")
    print(f"Resolved species B: '{species_b_raw}' -> {resolved['species_b']} ({resolved['species_b_match']} match)")
    print(f"Chain: {chain_raw} -> locus(es): {resolved['loci']}")
    print(f"Gene type: {resolved['gene_type_label']} -> region(s): {resolved['region_set']}")

    results = []
    for locus in resolved['loci']:
        print(f'\n=== {locus} ===')
        result = compare_two_species(
            records, resolved['species_a'], resolved['species_b'], locus,
            out_dir=out_dir, region_set=resolved['region_set'],
            gene_type_label=resolved['gene_type_label'],
        )
        results.append(result)
    return results

def main():
    args = build_argparser().parse_args()
    if args.query:
        try:
            species_a_raw, species_b_raw, chain_raw, gene_raw = parse_freeform_query(args.query)
        except QueryError as e:
            print(f'ERROR: {e}')
            sys.exit(1)
        print(f'Parsed freeform query -> species_a="{species_a_raw}", species_b="{species_b_raw}", '
              f'chain="{chain_raw}", gene="{gene_raw}"')
    else:
        if not args.species_a or not args.species_b:
            print('ERROR: provide --query OR both --species_a and --species_b')
            sys.exit(1)
        species_a_raw, species_b_raw, chain_raw, gene_raw = args.species_a, args.species_b, args.chain, args.gene

    run(species_a_raw, species_b_raw, chain_raw, gene_raw, args.bulk_path, args.out_dir)

if __name__ == '__main__':
    main()
