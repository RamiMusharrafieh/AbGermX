#!/usr/bin/env python3
"""
Runs the stochastic V(D)J simulation across all species with sufficient V/D/J
gene data, builds the resulting synthetic HCDR3 dataset, computes
biochemical/sequence features (same schema as the CDR1/CDR2 analysis), and
produces PCA + t-SNE embeddings plus a per-species length distribution.
"""
import argparse
import random

import _bootstrap  # noqa: F401  (puts ../src on sys.path)
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE

import plot_style
from config import DEFAULT_AA_BULK, DEFAULT_NT_BULK, ensure_out_dir
from imgt_parser import parse_bulk_fasta
from vdj_recombination import build_species_vdj_pools, generate_hcdr3s
from aa_properties import (mean_hydrophobicity, count_aromatic, count_polar,
                           count_acidic, count_basic, net_charge, aa_composition, ALL_20_AA)
from run_dimred import plot_embedding


def simulate_all_species(n_per_species=300, seed=42, locus='IGH',
                         aa_path=DEFAULT_AA_BULK, nt_path=DEFAULT_NT_BULK):
    aa_records = parse_bulk_fasta(aa_path)
    nt_records = parse_bulk_fasta(nt_path)
    species_list = sorted(set(r['species'] for r in aa_records if r['locus'] == locus))

    all_rows = []
    productivity_log = []
    rng = random.Random(seed)
    for species in species_list:
        pools = build_species_vdj_pools(aa_records, nt_records, species, locus=locus)
        n_v, n_d, n_j = len(pools['v_tails']), len(pools['d_segments']), len(pools['j_prefixes'])
        if not (n_v and n_d and n_j):
            productivity_log.append((species, n_v, n_d, n_j, 0, 0))
            continue
        results, attempts = generate_hcdr3s(pools, n_target=n_per_species, rng=rng)
        productivity_log.append((species, n_v, n_d, n_j, len(results), attempts))
        for r in results:
            r['species'] = species
            all_rows.append(r)

    df = pd.DataFrame(all_rows)
    log_df = pd.DataFrame(productivity_log,
                          columns=['species', 'n_V', 'n_D', 'n_J', 'n_productive', 'n_attempts'])
    return df, log_df


def hcdr3_features(seq):
    out = {
        'CDR3_length': len(seq),
        'CDR3_hydrophobicity': mean_hydrophobicity(seq) or 0.0,
        'CDR3_aromatic_count': count_aromatic(seq),
        'CDR3_polar_count': count_polar(seq),
        'CDR3_acidic_count': count_acidic(seq),
        'CDR3_basic_count': count_basic(seq),
        'CDR3_net_charge': net_charge(seq),
    }
    comp = aa_composition(seq)
    for aa in ALL_20_AA:
        out[f'CDR3_comp_{aa}'] = comp[aa]
    return out


def plot_length_distribution(df, out_path, n_per_species):
    """Per-species HCDR3 length distributions, ordered by median.

    All 24 species share one hue and the single species with the longest tail is
    picked out in a second: the chart has one finding to deliver (a right-tailed
    outlier), so colour is spent on emphasis rather than on 24-way identity that
    the axis labels already carry.
    """
    plot_style.use_style()
    order = (df.groupby('species').cdr3_aa.apply(lambda s: s.str.len().median())
             .sort_values().index.tolist())
    lengths = {sp: df.loc[df.species == sp, 'cdr3_aa'].str.len().to_numpy() for sp in order}
    # Emphasise whichever species reaches furthest into the long tail.
    standout = max(order, key=lambda sp: np.percentile(lengths[sp], 99))

    fig, ax = plt.subplots(figsize=(11.5, 9.6))
    positions = np.arange(len(order))
    parts = ax.violinplot([lengths[sp] for sp in order], positions=positions,
                          vert=False, widths=0.86, showextrema=False, showmedians=False)
    for sp, body in zip(order, parts['bodies']):
        emphasised = sp == standout
        body.set_facecolor(plot_style.CATEGORICAL[1] if emphasised else plot_style.CATEGORICAL[0])
        body.set_alpha(0.9 if emphasised else 0.32)
        body.set_edgecolor(plot_style.SURFACE)
        body.set_linewidth(1.0)
        body.set_zorder(3)

    # Median tick per species: a single value read is what the eye wants first.
    for pos, sp in zip(positions, order):
        median = float(np.median(lengths[sp]))
        ax.plot([median, median], [pos - 0.19, pos + 0.19], solid_capstyle='butt',
                color=plot_style.TEXT_PRIMARY if sp == standout else plot_style.TEXT_SECONDARY,
                linewidth=1.8, zorder=5)

    ax.set_yticks(positions)
    ax.set_yticklabels(order, fontstyle='italic')
    for tick, sp in zip(ax.get_yticklabels(), order):
        if sp == standout:
            tick.set_color(plot_style.TEXT_PRIMARY)
            tick.set_fontweight('bold')
    ax.set_ylim(-0.8, len(order) - 0.2)
    ax.set_xlabel('Simulated HCDR3 length (amino acids)')
    ax.set_ylabel('')
    plot_style.style_axes(ax, grid_axis='x')
    ax.grid(False, axis='y')

    standout_pos = positions[order.index(standout)]
    tail_x = float(np.percentile(lengths[standout], 99.5))
    ax.annotate(
        f'{standout} reaches {int(lengths[standout].max())} aa,\n'
        'consistent with its ultralong CDR-H3 biology',
        xy=(tail_x, standout_pos), xycoords='data',
        xytext=(tail_x * 0.58, standout_pos - 3.0), textcoords='data',
        fontsize=9.5, color=plot_style.TEXT_SECONDARY, ha='left', va='center',
        arrowprops=dict(arrowstyle='-', color=plot_style.AXIS_LINE, linewidth=0.9,
                        shrinkA=2, shrinkB=2, connectionstyle='arc3,rad=0.25'),
    )

    fig.subplots_adjust(left=0.235, right=0.975, top=0.885, bottom=0.075)
    plot_style.figure_title(
        fig, 'Simulated germline HCDR3 length by species',
        f'{n_per_species} productive recombinations per species, ordered by median length; '
        'vertical rule marks the median',
        x=0.02, y=0.945,
    )
    fig.savefig(out_path)
    plt.close(fig)
    return order, standout


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--n_per_species', type=int, default=300,
                   help='Productive HCDR3s to simulate per species (default: 300)')
    p.add_argument('--locus', default='IGH', help='Locus to simulate (default: IGH)')
    p.add_argument('--seed', type=int, default=42, help='Random seed (default: 42)')
    p.add_argument('--aa_path', default=DEFAULT_AA_BULK, help='IMGT/GENE-DB amino-acid bulk FASTA')
    p.add_argument('--nt_path', default=DEFAULT_NT_BULK, help='IMGT/GENE-DB nucleotide bulk FASTA')
    p.add_argument('--out_dir', default=None, help='Output directory')
    args = p.parse_args()

    out_dir = ensure_out_dir(args.out_dir)
    df, log_df = simulate_all_species(n_per_species=args.n_per_species, seed=args.seed,
                                      locus=args.locus, aa_path=args.aa_path, nt_path=args.nt_path)
    print(log_df.to_string(index=False))
    df.to_csv(f'{out_dir}/simulated_HCDR3_all_species.csv', index=False)
    print(f'\nTotal simulated productive HCDR3s: {len(df)} across {df.species.nunique()} species')

    feat_df = pd.DataFrame([hcdr3_features(s) for s in df.cdr3_aa])
    feat_df.to_csv(f'{out_dir}/simulated_HCDR3_features.csv', index=False)

    X = StandardScaler().fit_transform(feat_df.values)
    pca = PCA(n_components=2, random_state=args.seed)
    pca_coords = pca.fit_transform(X)
    perp = min(30, max(5, len(df) // 4))
    tsne = TSNE(n_components=2, perplexity=perp, random_state=args.seed,
                init='pca', learning_rate='auto')
    tsne_coords = tsne.fit_transform(X)

    df['pca_x'], df['pca_y'] = pca_coords[:, 0], pca_coords[:, 1]
    df['tsne_x'], df['tsne_y'] = tsne_coords[:, 0], tsne_coords[:, 1]
    df.to_csv(f'{out_dir}/simulated_HCDR3_dimred_coords.csv', index=False)

    plot_embedding(
        df, pca, perp,
        title='Simulated germline HCDR3s in biochemical feature space',
        subtitle=(f'{len(df)} productive recombinations across {df.species.nunique()} species, '
                  f'{feat_df.shape[1]} features'),
        out_path=f'{out_dir}/simulated_HCDR3_dimred_by_species.png',
    )
    plot_length_distribution(
        df, f'{out_dir}/simulated_HCDR3_length_distribution_by_species.png',
        n_per_species=args.n_per_species,
    )
    print('Saved plots and coordinate files.')


if __name__ == '__main__':
    main()
