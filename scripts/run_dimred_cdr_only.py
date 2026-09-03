#!/usr/bin/env python3
"""
Same pipeline as run_dimred.py, restricted to CDR1/CDR2-derived features only
(dropping the FR1+FR2+FR3 biochemical stats and composition vector).

The point of the comparison is to ask whether CDR properties alone separate
species, or whether the framework composition is doing most of that work.
"""
import argparse

import _bootstrap  # noqa: F401  (puts ../src on sys.path)
import matplotlib
matplotlib.use('Agg')
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE

from config import DEFAULT_AA_BULK, ensure_out_dir
from imgt_parser import parse_bulk_fasta
from build_annotation_table import build_table
from feature_matrix import build_feature_matrix
from run_dimred import plot_embedding

TITLE = '{locus}V germline genes in CDR1/CDR2-only feature space'
SUBTITLE = ('{n} QC-passed genes across {n_species} species, {n_feat} CDR-derived features '
            '(framework composition excluded)')


def run(locus='IGH', out_dir=None, perplexity=30, random_state=42, bulk_path=DEFAULT_AA_BULK):
    out_dir = ensure_out_dir(out_dir)
    records = parse_bulk_fasta(bulk_path)
    ann = build_table(records, locus=locus, species_list=None, functional_only=True,
                      allele_scope='representative')
    ann = ann[ann.qc_pass].reset_index(drop=True)
    print(f'{len(ann)} QC-passed genes across {ann.species.nunique()} species')

    feat_full, meta = build_feature_matrix(ann)
    cdr_cols = [c for c in feat_full.columns
                if c.startswith(('CDR1_', 'CDR2_', 'CDR_comp_'))]
    feat = feat_full[cdr_cols]
    print(f'CDR-only feature count: {feat.shape[1]} (vs {feat_full.shape[1]} full)')

    X = StandardScaler().fit_transform(feat.values)

    pca = PCA(n_components=2, random_state=random_state)
    pca_coords = pca.fit_transform(X)
    print(f'PCA explained variance: PC1={pca.explained_variance_ratio_[0]:.1%}, '
          f'PC2={pca.explained_variance_ratio_[1]:.1%}')

    perp = min(perplexity, max(5, len(ann) // 4))
    tsne = TSNE(n_components=2, perplexity=perp, random_state=random_state,
                init='pca', learning_rate='auto')
    tsne_coords = tsne.fit_transform(X)

    meta = meta.copy()
    meta['pca_x'], meta['pca_y'] = pca_coords[:, 0], pca_coords[:, 1]
    meta['tsne_x'], meta['tsne_y'] = tsne_coords[:, 0], tsne_coords[:, 1]
    meta.to_csv(f'{out_dir}/{locus}_Vgene_dimred_coords_CDRonly.csv', index=False)

    plot_embedding(
        meta, pca, perp,
        title=TITLE.format(locus=locus),
        subtitle=SUBTITLE.format(n=len(ann), n_species=meta.species.nunique(), n_feat=feat.shape[1]),
        out_path=f'{out_dir}/{locus}_Vgene_dimred_by_species_CDRonly.png',
    )
    print('Saved plot and coordinates.')
    return meta, feat


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--locus', default='IGH', help='IGH | IGK | IGL (default: IGH)')
    p.add_argument('--bulk_path', default=DEFAULT_AA_BULK, help='Path to IMGT/GENE-DB AA bulk FASTA')
    p.add_argument('--out_dir', default=None, help='Output directory')
    p.add_argument('--perplexity', type=int, default=30, help='t-SNE perplexity (default: 30)')
    args = p.parse_args()
    run(locus=args.locus, out_dir=args.out_dir, perplexity=args.perplexity, bulk_path=args.bulk_path)


if __name__ == '__main__':
    main()
