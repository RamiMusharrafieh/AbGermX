#!/usr/bin/env python3
"""
Projects a query CDR3 sequence (e.g. a designed HCDR3, or the full CDR3 loop of
a designed VHH/antibody) onto the simulated per-species germline HCDR3 feature
space, and ranks species by biochemical distance and sequence identity.

Requires run_vdj_simulation.py to have been run first: it reads that script's
simulated_HCDR3_all_species.csv output.
"""
import argparse

import _bootstrap  # noqa: F401  (puts ../src on sys.path)
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE

from config import OUT_DIR, ensure_out_dir
from compare_query_to_species import _projection_figure, N_HIGHLIGHT
from run_vdj_simulation import hcdr3_features
from compare_species import nw_identity

DEFAULT_SIM_CSV = str(OUT_DIR / 'simulated_HCDR3_all_species.csv')


def compare_cdr3_query(query_name, query_cdr3, sim_csv=DEFAULT_SIM_CSV, out_dir=None):
    out_dir = ensure_out_dir(out_dir)
    sim = pd.read_csv(sim_csv)

    feat_df = pd.DataFrame([hcdr3_features(s) for s in sim.cdr3_aa])
    scaler = StandardScaler().fit(feat_df.values)
    X = scaler.transform(feat_df.values)

    query_feat = pd.DataFrame([hcdr3_features(query_cdr3)])[feat_df.columns]
    Xq = scaler.transform(query_feat.values)

    sim = sim.copy()
    sim['feature_distance_to_query'] = np.linalg.norm(X - Xq, axis=1)
    sim['seq_identity_to_query'] = [nw_identity(query_cdr3, c) for c in sim.cdr3_aa]

    species_summary = sim.groupby('species').agg(
        n=('cdr3_aa', 'size'),
        mean_feature_distance=('feature_distance_to_query', 'mean'),
        min_feature_distance=('feature_distance_to_query', 'min'),
        mean_seq_identity=('seq_identity_to_query', 'mean'),
        max_seq_identity=('seq_identity_to_query', 'max'),
    ).reset_index().sort_values('mean_feature_distance')
    species_summary.to_csv(f'{out_dir}/{query_name}_CDR3_species_ranking.csv', index=False)

    highlight = species_summary.species.head(N_HIGHLIGHT).tolist()
    subtitle = (f'{len(sim)} simulated HCDR3s across {sim.species.nunique()} species; '
                f'the {N_HIGHLIGHT} species closest to the query are picked out')

    pca = PCA(n_components=2, random_state=42)
    pca_coords = pca.fit_transform(X)
    _projection_figure(
        pca_coords, pca.transform(Xq), sim.species, highlight, query_name,
        title=f'{query_name} CDR3 against simulated germline HCDR3 space',
        subtitle=subtitle,
        xlabel=f'PC1 ({pca.explained_variance_ratio_[0]:.1%} of variance)',
        ylabel=f'PC2 ({pca.explained_variance_ratio_[1]:.1%} of variance)',
        out_path=f'{out_dir}/{query_name}_CDR3_PCA_projection.png',
        point_size=12,
    )

    perp = min(30, max(5, len(sim) // 4))
    tsne = TSNE(n_components=2, perplexity=perp, random_state=42, init='pca',
                learning_rate='auto')
    tsne_coords = tsne.fit_transform(np.vstack([X, Xq]))
    _projection_figure(
        tsne_coords[:-1], tsne_coords[-1:], sim.species, highlight, query_name,
        title=f'{query_name} CDR3 against simulated germline HCDR3 space',
        subtitle=f'{subtitle}. t-SNE (perplexity {perp}) with the query included in the fit',
        xlabel='t-SNE 1', ylabel='t-SNE 2',
        out_path=f'{out_dir}/{query_name}_CDR3_tSNE_projection.png',
        point_size=12,
    )

    return species_summary


def main():
    from example_designed_loops import EXAMPLE_LOOPS
    from compare_query_to_species import extract_query_cdrs

    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--name', default='VHH_flu_01', help='Label for the query, used in filenames')
    p.add_argument('--cdr3', default=None,
                   help='Query CDR3 (default: extracted from the bundled VHH_flu_01 example)')
    p.add_argument('--sim_csv', default=DEFAULT_SIM_CSV,
                   help='simulated_HCDR3_all_species.csv from run_vdj_simulation.py')
    p.add_argument('--out_dir', default=None, help='Output directory')
    args = p.parse_args()

    cdr3 = args.cdr3 or extract_query_cdrs(EXAMPLE_LOOPS[args.name]['sequence'])['CDR3_full']
    print('Query CDR3:', cdr3)
    summary = compare_cdr3_query(args.name, cdr3, sim_csv=args.sim_csv, out_dir=args.out_dir)
    print(summary.head(10).to_string(index=False))


if __name__ == '__main__':
    main()
