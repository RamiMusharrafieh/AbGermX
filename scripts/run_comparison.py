#!/usr/bin/env python3
"""
Orchestrates one species-pair / locus comparison: repertoire summary, pairwise
identity matrix, best-ortholog table, and the clustered identity heatmap.
"""
import time

import _bootstrap  # noqa: F401  (puts ../src on sys.path)
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import seaborn as sns

import plot_style
from config import DEFAULT_AA_BULK, ensure_out_dir
from imgt_parser import parse_bulk_fasta
from compare_species import repertoire_summary, pairwise_identity_matrix, best_ortholog_table

MAX_FIG_DIM_IN = 30


def _label_size(n, scale=560, lo=3.5, hi=8.0):
    """Tick label size that shrinks with the number of ticks, within sane bounds."""
    return float(np.clip(scale / max(n, 1), lo, hi))


def plot_identity_clustermap(mat, species_a, species_b, locus, gene_type_label, out_path,
                             cmap=None, center=None):
    """Hierarchically clustered identity heatmap.

    Percent identity is a magnitude, so it gets a single-hue light-to-dark ramp
    rather than a multi-hue one: with one hue, "darker" reads as "more similar"
    without the reader having to consult the colour bar to order two cells.
    """
    plot_style.use_style()
    cmap = cmap if cmap is not None else plot_style.SEQUENTIAL_CMAP

    fig_w = min(MAX_FIG_DIM_IN, max(11, mat.shape[1] * 0.20))
    fig_h = min(MAX_FIG_DIM_IN, max(8, mat.shape[0] * 0.17))
    show_xticks = mat.shape[1] <= 150
    show_yticks = mat.shape[0] <= 150

    if center is None:
        norm = None
    else:
        norm = mcolors.TwoSlopeNorm(
            vmin=float(np.floor(mat.values.min())), vcenter=center,
            vmax=max(float(np.ceil(mat.values.max())), center + 1),
        )

    g = sns.clustermap(
        mat, cmap=cmap, figsize=(fig_w, fig_h), norm=norm,
        xticklabels=show_xticks, yticklabels=show_yticks,
        cbar_kws={'label': '% amino acid identity'},
        dendrogram_ratio=(0.13, 0.09), method='average',
        linewidths=0, rasterized=True,
    )
    g.fig.set_facecolor(plot_style.SURFACE)
    # Reserve a headroom band for the title/subtitle block above the dendrogram.
    # This re-lays out every axes in the clustermap gridspec, the colour bar
    # included, so the bar is positioned afterwards rather than via cbar_pos.
    g.gs.update(top=0.90)
    g.ax_cbar.set_position((0.012, 0.30, 0.22 / fig_w, 0.22))

    heat = g.ax_heatmap
    heat.set_xlabel(f'{species_b}  ·  {mat.shape[1]} {locus} {gene_type_label} genes',
                    fontsize=11, color=plot_style.TEXT_SECONDARY, labelpad=10)
    heat.set_ylabel(f'{species_a}  ·  {mat.shape[0]} {locus} {gene_type_label} genes',
                    fontsize=11, color=plot_style.TEXT_SECONDARY, labelpad=10)
    heat.tick_params(axis='x', labelsize=_label_size(mat.shape[1]), rotation=90,
                     length=0, pad=2, colors=plot_style.TEXT_SECONDARY)
    heat.tick_params(axis='y', labelsize=_label_size(mat.shape[0]),
                     length=0, pad=2, colors=plot_style.TEXT_SECONDARY)

    # Dendrograms are structure, not data: keep them as hairlines that recede.
    for ax in (g.ax_row_dendrogram, g.ax_col_dendrogram):
        ax.set_facecolor(plot_style.SURFACE)
        for coll in ax.collections:
            coll.set_linewidth(0.6)
            coll.set_color(plot_style.TEXT_MUTED)

    cbar = g.ax_cbar
    cbar.set_facecolor(plot_style.SURFACE)
    cbar.tick_params(labelsize=9, length=0, pad=3, colors=plot_style.TEXT_SECONDARY)
    cbar.yaxis.label.set_size(9.5)
    cbar.yaxis.label.set_color(plot_style.TEXT_SECONDARY)
    for spine in cbar.spines.values():
        spine.set_visible(False)

    plot_style.figure_title(
        g.fig,
        f'{species_a} vs {species_b}: {locus} {gene_type_label}-gene germline identity',
        'Needleman-Wunsch global alignment, one representative functional allele per gene; '
        'rows and columns ordered by average-linkage clustering',
        x=0.015, y=0.955,
    )
    g.savefig(out_path)
    plt.close('all')
    return g


def compare_two_species(records, species_a, species_b, locus, out_dir='.', label_a=None, label_b=None,
                        region_set=None, gene_type_label='V', cmap=None,
                        center=None):
    label_a = label_a or species_a.split()[0].lower()
    label_b = label_b or species_b.split()[0].lower()
    gene_slug = gene_type_label.strip().replace(' ', '').replace('-', '')
    prefix = f'{out_dir}/{label_a}_vs_{label_b}_{locus}_{gene_slug}'

    sum_a = repertoire_summary(records, species_a, locus, region_set)
    sum_b = repertoire_summary(records, species_b, locus, region_set)

    print(f'--- {species_a} {locus} ---')
    print(sum_a['functionality_counts'], 'unique functional genes:', sum_a['n_unique_functional_genes'])
    print(f'--- {species_b} {locus} ---')
    print(sum_b['functionality_counts'], 'unique functional genes:', sum_b['n_unique_functional_genes'])

    t0 = time.time()
    mat = pairwise_identity_matrix(records, species_a, species_b, locus, region_set)
    print(f'Identity matrix {mat.shape} computed in {time.time()-t0:.1f}s')

    if mat.shape[0] == 0 or mat.shape[1] == 0:
        print(f'No functional {gene_type_label} genes found for this species/locus combination '
              'on one or both sides, skipping matrix/plot outputs.')
        return dict(summary_a=sum_a, summary_b=sum_b, identity_matrix=mat,
                    best_orthologs=None, prefix=prefix)

    best = best_ortholog_table(mat, label_a, label_b)
    best.to_csv(f'{prefix}_best_orthologs.csv', index=False)
    mat.round(2).to_csv(f'{prefix}_identity_matrix.csv')

    print('Median best-match identity:', best.pct_identity.median())
    print('Mean best-match identity:', round(best.pct_identity.mean(), 2))

    plot_identity_clustermap(mat, species_a, species_b, locus, gene_type_label,
                             f'{prefix}_clustermap.png', cmap=cmap, center=center)

    return dict(summary_a=sum_a, summary_b=sum_b, identity_matrix=mat, best_orthologs=best,
                prefix=prefix)


def main():
    import argparse
    p = argparse.ArgumentParser(
        description='Compare two species at one locus, using exact IMGT species names. '
                    'For common names and freeform queries, use workflow.py instead.')
    p.add_argument('--species_a', default='Homo sapiens', help='Exact IMGT species name')
    p.add_argument('--species_b', default='Mus musculus', help='Exact IMGT species name')
    p.add_argument('--locus', default='IGH', help='IGH | IGK | IGL (default: IGH)')
    p.add_argument('--region', default='V-REGION',
                   help='IMGT region label to compare (default: V-REGION)')
    p.add_argument('--gene_label', default='V', help='Gene-type label used in titles and filenames')
    p.add_argument('--bulk_path', default=DEFAULT_AA_BULK, help='IMGT/GENE-DB AA bulk FASTA')
    p.add_argument('--out_dir', default=None, help='Output directory')
    args = p.parse_args()

    records = parse_bulk_fasta(args.bulk_path)
    compare_two_species(records, args.species_a, args.species_b, args.locus,
                        out_dir=ensure_out_dir(args.out_dir),
                        region_set={args.region}, gene_type_label=args.gene_label)


if __name__ == '__main__':
    main()
