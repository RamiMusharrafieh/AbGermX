#!/usr/bin/env python3
"""List all species available in the bulk IMGT file, with antibody (IG) gene record
counts and loci present. TCR loci are excluded and this tool only covers IGH/IGK/IGL."""
import sys
from collections import defaultdict
from imgt_parser import parse_bulk_fasta

DEFAULT_BULK = '/mnt/user-data/uploads/IMGTGENEDB-ReferenceSequences.fasta-AA-WithGaps-F_ORF_inframeP'

def main():
    bulk_path = sys.argv[1] if len(sys.argv) > 1 else DEFAULT_BULK
    records = parse_bulk_fasta(bulk_path)
    by_species = defaultdict(lambda: defaultdict(int))
    for r in records:
        by_species[r['species']][r['locus']] += 1
    for sp in sorted(by_species, key=lambda s: -sum(by_species[s].values())):
        loci = ', '.join(f'{l}:{c}' for l, c in sorted(by_species[sp].items(), key=lambda x:-x[1]))
        print(f'{sp:35s} {loci}')

if __name__ == '__main__':
    main()
