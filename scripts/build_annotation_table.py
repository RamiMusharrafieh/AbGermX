"""
Builds a CSV of FR1/CDR1/FR2/CDR2/FR3/(partial-CDR3) annotations for antibody
heavy-chain (or any locus) V genes, across one or more species.
"""
import pandas as pd
from imgt_parser import parse_bulk_fasta, norm_functionality
from annotate_regions import annotate_sequence
from aa_properties import mean_hydrophobicity, count_aromatic, count_polar

def build_table(records, locus='IGH', species_list=None, functional_only=True, allele_scope='all'):
    """
    allele_scope: 'all' (every allele) or 'representative' (one *01-preferred allele per gene)
    """
    recs = [r for r in records if r['locus'] == locus and r['region'] == 'V-REGION']
    if species_list is not None:
        recs = [r for r in recs if r['species'] in species_list]
    if functional_only:
        recs = [r for r in recs if norm_functionality(r['functionality']) == 'Functional']

    if allele_scope == 'representative':
        from collections import defaultdict
        by_gene = defaultdict(list)
        for r in recs:
            by_gene[(r['species'], r['gene'])].append(r)
        chosen = []
        for (sp, gene), rs in by_gene.items():
            rs01 = [r for r in rs if r['allele'] == '01']
            chosen.append(rs01[0] if rs01 else rs[0])
        recs = chosen

    rows = []
    for r in recs:
        ann = annotate_sequence(r['seq'])
        cdr1_seq = ann.cdr1_ungapped
        cdr2_seq = ann.cdr2_ungapped
        cdr1_cdr2_seq = cdr1_seq + cdr2_seq
        rows.append(dict(
            species=r['species'], locus=r['locus'], gene=r['gene'], allele=r['allele'],
            gene_allele=r['gene_allele'], accession=r['acc'], functionality=r['functionality'].strip(),
            gapped_len=r['raw_len'],
            FR1=ann.fr1_ungapped, CDR1=cdr1_seq, FR2=ann.fr2_ungapped,
            CDR2=cdr2_seq, FR3=ann.fr3_ungapped, CDR3_partial=ann.cdr3_partial_ungapped,
            CDR1_length=len(cdr1_seq), CDR2_length=len(cdr2_seq),
            CDR1_hydrophobicity=mean_hydrophobicity(cdr1_seq),
            CDR2_hydrophobicity=mean_hydrophobicity(cdr2_seq),
            CDR1_CDR2_hydrophobicity=mean_hydrophobicity(cdr1_cdr2_seq),
            CDR1_aromatic_count=count_aromatic(cdr1_seq),
            CDR2_aromatic_count=count_aromatic(cdr2_seq),
            CDR1_CDR2_aromatic_count=count_aromatic(cdr1_cdr2_seq),
            CDR1_polar_count=count_polar(cdr1_seq),
            CDR2_polar_count=count_polar(cdr2_seq),
            CDR1_CDR2_polar_count=count_polar(cdr1_cdr2_seq),
            anchor_cys23_ok=ann.anchor_cys23_ok, anchor_trp41_ok=ann.anchor_trp41_ok,
            anchor_cys104_ok=ann.anchor_cys104_ok, qc_pass=ann.qc_pass,
            cdr1_insertion=ann.cdr1_insertion, cdr2_insertion=ann.cdr2_insertion,
        ))
    return pd.DataFrame(rows)

if __name__ == '__main__':
    import sys
    bulk_path = sys.argv[1] if len(sys.argv) > 1 else '/mnt/user-data/uploads/IMGTGENEDB-ReferenceSequences.fasta-AA-WithGaps-F_ORF_inframeP'
    records = parse_bulk_fasta(bulk_path)
    df = build_table(records, locus='IGH', species_list=['Homo sapiens'], functional_only=True, allele_scope='all')
    print(df.shape)
    print(df.head(10).to_string())
    print()
    print('QC pass rate:', df.qc_pass.mean())
    print(df[~df.qc_pass][['gene_allele','functionality','gapped_len']].head(10))
