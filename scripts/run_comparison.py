import sys, time
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import seaborn as sns

from imgt_parser import parse_bulk_fasta
from compare_species import repertoire_summary, pairwise_identity_matrix, best_ortholog_table

BLUE_YELLOW_RED = mcolors.LinearSegmentedColormap.from_list(
    'blue_yellow_red', ['#2138C4', '#FFF200', '#D42A2A']
)

def compare_two_species(records, species_a, species_b, locus, out_dir='.', label_a=None, label_b=None,
                         region_set=None, gene_type_label='V', cmap='viridis', cbar_pos=(0.0, 0.3, 0.015, 0.4),
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
        print(f'No functional {gene_type_label} genes found for this species/locus combination on one or both sides -- skipping matrix/plot outputs.')
        return dict(summary_a=sum_a, summary_b=sum_b, identity_matrix=mat, best_orthologs=None, prefix=prefix)

    best = best_ortholog_table(mat, label_a, label_b)
    best.to_csv(f'{prefix}_best_orthologs.csv', index=False)
    mat.round(2).to_csv(f'{prefix}_identity_matrix.csv')

    print('Median best-match identity:', best.pct_identity.median())
    print('Mean best-match identity:', round(best.pct_identity.mean(),2))

    max_dim_in = 34  
    fig_w = min(max_dim_in, max(10, mat.shape[1]*0.22))
    fig_h = min(max_dim_in, max(8, mat.shape[0]*0.18))
    show_xticks = mat.shape[1] <= 150
    show_yticks = mat.shape[0] <= 150

    if center is None:
        norm = None
    else:
        data_min = float(np.floor(mat.values.min()))
        data_max = float(np.ceil(mat.values.max()))
        vmax = max(data_max, center + 1)
        norm = mcolors.TwoSlopeNorm(vmin=data_min, vcenter=center, vmax=vmax)

    g = sns.clustermap(
        mat, cmap=cmap, figsize=(fig_w, fig_h),
        xticklabels=show_xticks, yticklabels=show_yticks,
        cbar_kws={'label': '% amino acid identity'},
        cbar_pos=cbar_pos, norm=norm,
        dendrogram_ratio=(0.16, 0.1), method='average',
    )
    g.ax_heatmap.set_xlabel(f'{species_b} {locus} {gene_type_label} genes (n={mat.shape[1]})', fontsize=11)
    g.ax_heatmap.set_ylabel(f'{species_a} {locus} {gene_type_label} genes (n={mat.shape[0]})', fontsize=11)
    if show_xticks:
        g.ax_heatmap.tick_params(axis='x', labelsize=max(3, 500/mat.shape[1]), rotation=90)
    if show_yticks:
        g.ax_heatmap.tick_params(axis='y', labelsize=max(3, 500/mat.shape[0]))
    g.fig.suptitle(f'{species_a} vs {species_b}: {locus} {gene_type_label}-gene germline identity', y=1.02, fontsize=13)
    g.savefig(f'{prefix}_clustermap.png', dpi=150, bbox_inches='tight')
    plt.close('all')

    return dict(summary_a=sum_a, summary_b=sum_b, identity_matrix=mat, best_orthologs=best,
                prefix=prefix)

if __name__ == '__main__':
    records = parse_bulk_fasta('/mnt/user-data/uploads/IMGTGENEDB-ReferenceSequences.fasta-AA-WithGaps-F_ORF_inframeP')
    result = compare_two_species(records, 'Homo sapiens', 'Mus musculus', 'IGH', out_dir='/home/claude/imgt_tool')
