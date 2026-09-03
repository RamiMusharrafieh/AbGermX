#!/usr/bin/env python3
"""
Lists every species available in an IMGT/GENE-DB bulk file, with antibody (IG)
gene record counts per locus. TCR loci are excluded: this toolkit covers
IGH/IGK/IGL only.
"""
import argparse
from collections import defaultdict

import _bootstrap  # noqa: F401  (puts ../src on sys.path)
from config import DEFAULT_AA_BULK
from imgt_parser import parse_bulk_fasta


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('bulk_path', nargs='?', default=DEFAULT_AA_BULK,
                   help='Path to IMGT/GENE-DB bulk FASTA (default: data/ copy)')
    args = p.parse_args()

    records = parse_bulk_fasta(args.bulk_path)
    by_species = defaultdict(lambda: defaultdict(int))
    for r in records:
        by_species[r['species']][r['locus']] += 1
    for sp in sorted(by_species, key=lambda s: -sum(by_species[s].values())):
        loci = ', '.join(f'{locus}:{count}'
                         for locus, count in sorted(by_species[sp].items(), key=lambda x: -x[1]))
        print(f'{sp:35s} {loci}')


if __name__ == '__main__':
    main()
