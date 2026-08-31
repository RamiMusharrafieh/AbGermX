"""
Builds a numeric feature matrix (one row per V gene) for dimensionality
reduction (PCA / t-SNE), combining:
  - per-region biochemical summary stats (length, hydrophobicity, aromatic,
    polar, acidic, basic, net charge) for CDR1, CDR2, FR1+FR2+FR3 combined
  - full 20-amino-acid composition vectors for CDR1+CDR2 combined and for
    FR1+FR2+FR3 combined
Requires an annotation DataFrame from build_annotation_table.build_table().
"""
import pandas as pd
from aa_properties import (
    mean_hydrophobicity, count_aromatic, count_polar,
    count_acidic, count_basic, net_charge, aa_composition, ALL_20_AA
)

def _region_stats(seq, prefix, out):
    out[f'{prefix}_length'] = len(seq)
    out[f'{prefix}_hydrophobicity'] = mean_hydrophobicity(seq) or 0.0
    out[f'{prefix}_aromatic_count'] = count_aromatic(seq)
    out[f'{prefix}_polar_count'] = count_polar(seq)
    out[f'{prefix}_acidic_count'] = count_acidic(seq)
    out[f'{prefix}_basic_count'] = count_basic(seq)
    out[f'{prefix}_net_charge'] = net_charge(seq)

def build_feature_matrix(annotation_df):
    """
    annotation_df: output of build_annotation_table.build_table() (must have
    FR1, CDR1, FR2, CDR2, FR3 columns).
    Returns (feature_df, meta_df): feature_df is numeric-only (for the reducer),
    meta_df carries species/gene_allele/etc. for labeling points afterward,
    aligned by row.
    """
    rows = []
    for _, r in annotation_df.iterrows():
        cdr1, cdr2 = r['CDR1'], r['CDR2']
        cdr_combined = cdr1 + cdr2
        fr_combined = r['FR1'] + r['FR2'] + r['FR3']

        out = {}
        _region_stats(cdr1, 'CDR1', out)
        _region_stats(cdr2, 'CDR2', out)
        _region_stats(fr_combined, 'FR', out)

        comp_cdr = aa_composition(cdr_combined)
        for aa in ALL_20_AA:
            out[f'CDR_comp_{aa}'] = comp_cdr[aa]
        comp_fr = aa_composition(fr_combined)
        for aa in ALL_20_AA:
            out[f'FR_comp_{aa}'] = comp_fr[aa]

        rows.append(out)

    feature_df = pd.DataFrame(rows)
    meta_cols = ['species', 'gene', 'allele', 'gene_allele', 'locus']
    meta_df = annotation_df[[c for c in meta_cols if c in annotation_df.columns]].reset_index(drop=True)
    return feature_df, meta_df

if __name__ == '__main__':
    from imgt_parser import parse_bulk_fasta
    from build_annotation_table import build_table
    records = parse_bulk_fasta('/mnt/user-data/uploads/IMGTGENEDB-ReferenceSequences.fasta-AA-WithGaps-F_ORF_inframeP')
    ann = build_table(records, locus='IGH', species_list=None, functional_only=True, allele_scope='representative')
    ann = ann[ann.qc_pass]  # drop QC-failed (partial) sequences for cleaner features
    feat, meta = build_feature_matrix(ann)
    print('Feature matrix shape:', feat.shape)
    print('Features:', list(feat.columns))
    print(feat.head())
