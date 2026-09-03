"""
General-purpose parser for IMGT/GENE-DB bulk FASTA files.
Works on the standard 15-field pipe-delimited IMGT header format.
"""
import re
from collections import defaultdict

IG_LOCI = {'IGH', 'IGK', 'IGL'}  # antibody (immunoglobulin) loci only; TR* and IGI excluded

def parse_bulk_fasta(path, region_filter=None, ig_only=True):
    """
    Parse an IMGT/GENE-DB bulk FASTA file.
    region_filter: if set, only keep entries whose field-5 (exon/region name) is in
                   this set/list (e.g. {'V-REGION'}). None keeps all region types.
    ig_only: if True (default), keep only antibody loci (IGH/IGK/IGL) and drop
             T-cell receptor (TRA/TRB/TRG/TRD) and other non-Ig entries (e.g. fish IGI).
    Returns list of dict records.
    """
    if isinstance(region_filter, str):
        region_filter = {region_filter}
    records = []
    with open(path, encoding='utf-8', errors='replace') as f:
        header = None
        seq_lines = []
        def flush():
            if header is None:
                return
            fields = header.split('|')
            acc = fields[0] if len(fields) > 0 else ''
            gene_allele = fields[1] if len(fields) > 1 else ''
            species_full = fields[2] if len(fields) > 2 else ''
            functionality = fields[3] if len(fields) > 3 else ''
            region = fields[4] if len(fields) > 4 else ''
            if region_filter is not None and region not in region_filter:
                return
            gene = gene_allele.split('*')[0] if gene_allele else ''
            allele = gene_allele.split('*')[1] if '*' in gene_allele else ''
            locus_m = re.match(r'(IG[HKL]|TR[ABGD]|IGI)', gene)
            locus = locus_m.group(1) if locus_m else gene[:3]
            if ig_only and locus not in IG_LOCI:
                return
            species = species_full.split('_')[0]
            strain = species_full[len(species):].lstrip('_') if '_' in species_full else ''
            subgroup_m = re.match(r'([A-Z]+\d+)', gene)
            subgroup = subgroup_m.group(1) if subgroup_m else gene
            seq = ''.join(seq_lines)
            records.append(dict(
                acc=acc, gene_allele=gene_allele, gene=gene, allele=allele,
                species=species, strain=strain, species_full=species_full,
                functionality=functionality, region=region, locus=locus,
                subgroup=subgroup, seq=seq, raw_len=len(seq)
            ))
        for line in f:
            line = line.rstrip('\n')
            if line.startswith('>'):
                flush()
                header = line[1:]
                seq_lines = []
            else:
                if line.strip():
                    seq_lines.append(line.strip())
        flush()
    return records

def norm_functionality(f):
    f = f.strip()
    if f in ('F', '(F)', '[F]'):
        return 'Functional'
    if f == 'ORF':
        return 'ORF'
    if f.startswith('P') or '(P)' in f or '[P]' in f:
        return 'Pseudogene'
    return 'Other'

def list_species(records):
    counts = defaultdict(int)
    for r in records:
        counts[r['species']] += 1
    return dict(sorted(counts.items(), key=lambda x: -x[1]))

if __name__ == '__main__':
    import sys
    from config import DEFAULT_AA_BULK
    path = sys.argv[1] if len(sys.argv) > 1 else DEFAULT_AA_BULK
    recs = parse_bulk_fasta(path)
    print(f'Parsed {len(recs)} IG gene records')
    print('Species:', list_species(recs))
