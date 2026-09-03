#!/usr/bin/env python3
"""
Runs PCA and t-SNE on the V-gene feature matrix and plots the embedding,
coloured by species.
"""
import argparse

import _bootstrap  # noqa: F401  (puts ../src on sys.path)
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE

import plot_style
from config import DEFAULT_AA_BULK, ensure_out_dir
from imgt_parser import parse_bulk_fasta
from build_annotation_table import build_table
from feature_matrix import build_feature_matrix

TITLE = '{locus}V germline genes in biochemical feature space'
SUBTITLE = ('{n} QC-passed genes across {n_species} species, {n_feat} features '
            '(region length, hydrophobicity, charge and amino-acid composition)')


def plot_embedding(meta, pca, perp, title, subtitle, out_path, n_species_cols=6):
    """Two panels (PCA, t-SNE) sharing one legend.

    Species identity is carried by colour AND marker shape together: there are
    more species here than any palette can separate by hue alone, so the two
    channels are combined rather than generating extra hues.
    """
    plot_style.use_style()
    species_list = sorted(meta.species.unique())
    styles = plot_style.species_styles(species_list)

    fig, axes = plt.subplots(1, 2, figsize=(15, 8.4))
    for sp in species_list:
        sub = meta[meta.species == sp]
        plot_style.scatter_species(axes[0], sub.pca_x, sub.pca_y, styles[sp])
        plot_style.scatter_species(axes[1], sub.tsne_x, sub.tsne_y, styles[sp], label=sp)

    axes[0].set_title(f'PCA  ·  PC1 {pca.explained_variance_ratio_[0]:.1%}, '
                      f'PC2 {pca.explained_variance_ratio_[1]:.1%} of variance', loc='left')
    axes[0].set_xlabel('PC1')
    axes[0].set_ylabel('PC2')
    axes[1].set_title(f't-SNE  ·  perplexity {perp}', loc='left')
    axes[1].set_xlabel('t-SNE 1')
    axes[1].set_ylabel('t-SNE 2')
    for ax in axes:
        plot_style.style_axes(ax)

    handles, labels = axes[1].get_legend_handles_labels()
    fig.subplots_adjust(left=0.055, right=0.985, top=0.845, bottom=0.255, wspace=0.16)
    plot_style.figure_title(fig, title, subtitle, x=0.055, y=0.945)
    plot_style.species_legend(fig, handles, labels, ncol=n_species_cols, y=0.185)
    fig.savefig(out_path)
    plt.close(fig)


def run(locus='IGH', out_dir=None, perplexity=30, random_state=42, bulk_path=DEFAULT_AA_BULK):
    out_dir = ensure_out_dir(out_dir)
    records = parse_bulk_fasta(bulk_path)
    ann = build_table(records, locus=locus, species_list=None, functional_only=True,
                      allele_scope='representative')
    ann = ann[ann.qc_pass].reset_index(drop=True)
    print(f'{len(ann)} QC-passed genes across {ann.species.nunique()} species')

    feat, meta = build_feature_matrix(ann)
    X = StandardScaler().fit_transform(feat.values)

    pca = PCA(n_components=2, random_state=random_state)
    pca_coords = pca.fit_transform(X)
    print(f'PCA explained variance: PC1={pca.explained_variance_ratio_[0]:.1%}, '
          f'PC2={pca.explained_variance_ratio_[1]:.1%}')

    # perplexity capped below n_samples/3 as a safe default
    perp = min(perplexity, max(5, len(ann) // 4))
    tsne = TSNE(n_components=2, perplexity=perp, random_state=random_state,
                init='pca', learning_rate='auto')
    tsne_coords = tsne.fit_transform(X)

    meta = meta.copy()
    meta['pca_x'], meta['pca_y'] = pca_coords[:, 0], pca_coords[:, 1]
    meta['tsne_x'], meta['tsne_y'] = tsne_coords[:, 0], tsne_coords[:, 1]
    meta.to_csv(f'{out_dir}/{locus}_Vgene_dimred_coords.csv', index=False)

    plot_embedding(
        meta, pca, perp,
        title=TITLE.format(locus=locus),
        subtitle=SUBTITLE.format(n=len(ann), n_species=meta.species.nunique(), n_feat=feat.shape[1]),
        out_path=f'{out_dir}/{locus}_Vgene_dimred_by_species.png',
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
