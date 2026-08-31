"""
Projects a query CDR3 sequence (e.g. a designed HCDR3, or the full CDR3 loop of
a designed VHH/antibody) onto the simulated per-species germline HCDR3 feature
space, and ranks species by biochemical distance + sequence identity to the query.
"""
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE

from run_vdj_simulation import hcdr3_features
from compare_species import nw_identity

def compare_cdr3_query(query_name, query_cdr3, sim_csv='/mnt/user-data/outputs/simulated_HCDR3_all_species.csv',
                        out_dir='/mnt/user-data/outputs'):
    sim = pd.read_csv(sim_csv)

    feat_rows = [hcdr3_features(s) for s in sim.cdr3_aa]
    feat_df = pd.DataFrame(feat_rows)
    scaler = StandardScaler().fit(feat_df.values)
    X = scaler.transform(feat_df.values)

    query_feat = pd.DataFrame([hcdr3_features(query_cdr3)])[feat_df.columns]
    Xq = scaler.transform(query_feat.values)

    dists = np.linalg.norm(X - Xq, axis=1)
    sim = sim.copy()
    sim['feature_distance_to_query'] = dists
    sim['seq_identity_to_query'] = [nw_identity(query_cdr3, c) for c in sim.cdr3_aa]

    species_summary = sim.groupby('species').agg(
        n=('cdr3_aa', 'size'),
        mean_feature_distance=('feature_distance_to_query', 'mean'),
        min_feature_distance=('feature_distance_to_query', 'min'),
        mean_seq_identity=('seq_identity_to_query', 'mean'),
        max_seq_identity=('seq_identity_to_query', 'max'),
    ).reset_index().sort_values('mean_feature_distance')
    species_summary.to_csv(f'{out_dir}/{query_name}_CDR3_species_ranking.csv', index=False)

    # PCA projection
    pca = PCA(n_components=2, random_state=42)
    pca_coords = pca.fit_transform(X)
    query_pca = pca.transform(Xq)

    species_list = sorted(sim.species.unique())
    cmap = plt.get_cmap('tab20' if len(species_list) <= 20 else 'gist_ncar')
    color_map = {sp: cmap(i / max(1, len(species_list)-1)) for i, sp in enumerate(species_list)}

    fig, ax = plt.subplots(figsize=(11, 9))
    for sp in species_list:
        mask = (sim.species == sp).values
        ax.scatter(pca_coords[mask, 0], pca_coords[mask, 1], s=10, alpha=0.45, color=color_map[sp], label=sp)
    ax.scatter(query_pca[:, 0], query_pca[:, 1], s=280, marker='*', color='black',
               edgecolor='white', linewidth=1.2, zorder=10, label=f'{query_name} (query)')
    ax.annotate(f'{query_name}\n"{query_cdr3}"', (query_pca[0, 0], query_pca[0, 1]),
                textcoords='offset points', xytext=(10, 8), fontsize=9, fontweight='bold')
    ax.set_xlabel(f'PC1 ({pca.explained_variance_ratio_[0]:.1%})')
    ax.set_ylabel(f'PC2 ({pca.explained_variance_ratio_[1]:.1%})')
    ax.set_title(f'{query_name} CDR3 projected onto simulated germline HCDR3 feature space')
    ax.legend(bbox_to_anchor=(1.02, 1), loc='upper left', fontsize=6.5, title='Species')
    fig.tight_layout()
    fig.savefig(f'{out_dir}/{query_name}_CDR3_PCA_projection.png', dpi=150, bbox_inches='tight')
    plt.close('all')

    # t-SNE with query merged in
    X_with_query = np.vstack([X, Xq])
    perp = min(30, max(5, len(sim)//4))
    tsne = TSNE(n_components=2, perplexity=perp, random_state=42, init='pca', learning_rate='auto')
    tsne_coords = tsne.fit_transform(X_with_query)
    tsne_natural, tsne_query = tsne_coords[:-1], tsne_coords[-1:]

    fig, ax = plt.subplots(figsize=(11, 9))
    for sp in species_list:
        mask = (sim.species == sp).values
        ax.scatter(tsne_natural[mask, 0], tsne_natural[mask, 1], s=10, alpha=0.45, color=color_map[sp], label=sp)
    ax.scatter(tsne_query[:, 0], tsne_query[:, 1], s=280, marker='*', color='black',
               edgecolor='white', linewidth=1.2, zorder=10, label=f'{query_name} (query)')
    ax.annotate(f'{query_name}\n"{query_cdr3}"', (tsne_query[0, 0], tsne_query[0, 1]),
                textcoords='offset points', xytext=(10, 8), fontsize=9, fontweight='bold')
    ax.set_xlabel('t-SNE 1'); ax.set_ylabel('t-SNE 2')
    ax.set_title(f'{query_name} CDR3 in simulated germline HCDR3 feature space (t-SNE, query merged into fit)')
    ax.legend(bbox_to_anchor=(1.02, 1), loc='upper left', fontsize=6.5, title='Species')
    fig.tight_layout()
    fig.savefig(f'{out_dir}/{query_name}_CDR3_tSNE_projection.png', dpi=150, bbox_inches='tight')
    plt.close('all')

    return species_summary

if __name__ == '__main__':
    from example_designed_loops import EXAMPLE_LOOPS
    from compare_query_to_species import extract_query_cdrs
    vhh = EXAMPLE_LOOPS['VHH_flu_01']
    extracted = extract_query_cdrs(vhh['sequence'])
    print('Query CDR3:', extracted['CDR3_full'])
    summary = compare_cdr3_query('VHH_flu_01', extracted['CDR3_full'])
    print(summary.head(10).to_string())
