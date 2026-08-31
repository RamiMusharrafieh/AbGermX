"""
Runs PCA and t-SNE on the V-gene feature matrix and plots the embedding,
colored by species.
"""
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE

from imgt_parser import parse_bulk_fasta
from build_annotation_table import build_table
from feature_matrix import build_feature_matrix

def run(locus='IGH', out_dir='/mnt/user-data/outputs', perplexity=30, random_state=42,
        bulk_path='/mnt/user-data/uploads/IMGTGENEDB-ReferenceSequences.fasta-AA-WithGaps-F_ORF_inframeP'):
    records = parse_bulk_fasta(bulk_path)
    ann = build_table(records, locus=locus, species_list=None, functional_only=True, allele_scope='representative')
    ann = ann[ann.qc_pass].reset_index(drop=True)
    print(f'{len(ann)} QC-passed genes across {ann.species.nunique()} species')

    feat, meta = build_feature_matrix(ann)
    X = StandardScaler().fit_transform(feat.values)

    # PCA
    pca = PCA(n_components=2, random_state=random_state)
    pca_coords = pca.fit_transform(X)
    print(f'PCA explained variance: PC1={pca.explained_variance_ratio_[0]:.1%}, PC2={pca.explained_variance_ratio_[1]:.1%}')

    # t-SNE (perplexity capped below n_samples/3 as a safe default)
    perp = min(perplexity, max(5, len(ann)//4))
    tsne = TSNE(n_components=2, perplexity=perp, random_state=random_state, init='pca', learning_rate='auto')
    tsne_coords = tsne.fit_transform(X)

    meta = meta.copy()
    meta['pca_x'], meta['pca_y'] = pca_coords[:,0], pca_coords[:,1]
    meta['tsne_x'], meta['tsne_y'] = tsne_coords[:,0], tsne_coords[:,1]
    meta.to_csv(f'{out_dir}/{locus}_Vgene_dimred_coords.csv', index=False)

    # Plot: species as color, two panels (PCA left, t-SNE right)
    species_list = sorted(meta.species.unique())
    cmap = plt.get_cmap('tab20' if len(species_list) <= 20 else 'gist_ncar')
    color_map = {sp: cmap(i / max(1,len(species_list)-1)) for i, sp in enumerate(species_list)}

    fig, axes = plt.subplots(1, 2, figsize=(20, 9))
    for sp in species_list:
        sub = meta[meta.species == sp]
        axes[0].scatter(sub.pca_x, sub.pca_y, s=14, alpha=0.75, color=color_map[sp], label=sp)
        axes[1].scatter(sub.tsne_x, sub.tsne_y, s=14, alpha=0.75, color=color_map[sp], label=sp)

    axes[0].set_title(f'PCA (PC1 {pca.explained_variance_ratio_[0]:.1%}, PC2 {pca.explained_variance_ratio_[1]:.1%})')
    axes[0].set_xlabel('PC1'); axes[0].set_ylabel('PC2')
    axes[1].set_title(f't-SNE (perplexity={perp})')
    axes[1].set_xlabel('t-SNE 1'); axes[1].set_ylabel('t-SNE 2')
    fig.suptitle(f'{locus}V gene biochemical/sequence feature space, colored by species (n={len(ann)} genes, {feat.shape[1]} features)', fontsize=13)
    axes[1].legend(bbox_to_anchor=(1.02, 1), loc='upper left', fontsize=7, ncol=1, title='Species')
    fig.tight_layout()
    fig.savefig(f'{out_dir}/{locus}_Vgene_dimred_by_species.png', dpi=150, bbox_inches='tight')
    plt.close('all')
    print('Saved plot and coordinates.')
    return meta, feat

if __name__ == '__main__':
    run()
