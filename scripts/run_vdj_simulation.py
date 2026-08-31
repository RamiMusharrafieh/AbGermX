"""
Runs the stochastic VDJ simulation across all species with sufficient V/D/J gene
data, builds the resulting synthetic HCDR3 dataset, computes biochemical/sequence
features (same style as the CDR1/CDR2 analysis), and runs PCA + t-SNE colored by
species.
"""
import random
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE

from imgt_parser import parse_bulk_fasta
from vdj_recombination import build_species_vdj_pools, generate_hcdr3s
from aa_properties import mean_hydrophobicity, count_aromatic, count_polar, count_acidic, count_basic, net_charge, aa_composition, ALL_20_AA

def simulate_all_species(n_per_species=300, seed=42, locus='IGH',
                          aa_path='/mnt/user-data/uploads/IMGTGENEDB-ReferenceSequences.fasta-AA-WithGaps-F_ORF_inframeP',
                          nt_path='/mnt/user-data/uploads/IMGTGENEDB-ReferenceSequences.fasta-nt-WithGaps-F_ORF_inframeP'):
    aa_records = parse_bulk_fasta(aa_path)
    nt_records = parse_bulk_fasta(nt_path)
    species_list = sorted(set(r['species'] for r in aa_records if r['locus']==locus))

    all_rows = []
    productivity_log = []
    rng = random.Random(seed)
    for species in species_list:
        pools = build_species_vdj_pools(aa_records, nt_records, species, locus=locus)
        if not (pools['v_tails'] and pools['d_segments'] and pools['j_prefixes']):
            productivity_log.append((species, len(pools['v_tails']), len(pools['d_segments']), len(pools['j_prefixes']), 0, 0))
            continue
        results, attempts = generate_hcdr3s(pools, n_target=n_per_species, rng=rng)
        productivity_log.append((species, len(pools['v_tails']), len(pools['d_segments']), len(pools['j_prefixes']), len(results), attempts))
        for r in results:
            r['species'] = species
            all_rows.append(r)

    df = pd.DataFrame(all_rows)
    log_df = pd.DataFrame(productivity_log, columns=['species','n_V','n_D','n_J','n_productive','n_attempts'])
    return df, log_df

def hcdr3_features(seq):
    out = {}
    out['CDR3_length'] = len(seq)
    out['CDR3_hydrophobicity'] = mean_hydrophobicity(seq) or 0.0
    out['CDR3_aromatic_count'] = count_aromatic(seq)
    out['CDR3_polar_count'] = count_polar(seq)
    out['CDR3_acidic_count'] = count_acidic(seq)
    out['CDR3_basic_count'] = count_basic(seq)
    out['CDR3_net_charge'] = net_charge(seq)
    comp = aa_composition(seq)
    for aa in ALL_20_AA:
        out[f'CDR3_comp_{aa}'] = comp[aa]
    return out

if __name__ == '__main__':
    out_dir = '/mnt/user-data/outputs'
    df, log_df = simulate_all_species(n_per_species=300)
    print(log_df.to_string(index=False))
    df.to_csv(f'{out_dir}/simulated_HCDR3_all_species.csv', index=False)
    print(f'\nTotal simulated productive HCDR3s: {len(df)} across {df.species.nunique()} species')

    feat_rows = [hcdr3_features(s) for s in df.cdr3_aa]
    feat_df = pd.DataFrame(feat_rows)
    feat_df.to_csv(f'{out_dir}/simulated_HCDR3_features.csv', index=False)

    X = StandardScaler().fit_transform(feat_df.values)
    pca = PCA(n_components=2, random_state=42)
    pca_coords = pca.fit_transform(X)
    perp = min(30, max(5, len(df)//4))
    tsne = TSNE(n_components=2, perplexity=perp, random_state=42, init='pca', learning_rate='auto')
    tsne_coords = tsne.fit_transform(X)

    df['pca_x'], df['pca_y'] = pca_coords[:,0], pca_coords[:,1]
    df['tsne_x'], df['tsne_y'] = tsne_coords[:,0], tsne_coords[:,1]
    df.to_csv(f'{out_dir}/simulated_HCDR3_dimred_coords.csv', index=False)

    species_list = sorted(df.species.unique())
    cmap = plt.get_cmap('tab20' if len(species_list) <= 20 else 'gist_ncar')
    color_map = {sp: cmap(i / max(1,len(species_list)-1)) for i, sp in enumerate(species_list)}

    fig, axes = plt.subplots(1, 2, figsize=(20, 9))
    for sp in species_list:
        sub = df[df.species == sp]
        axes[0].scatter(sub.pca_x, sub.pca_y, s=10, alpha=0.5, color=color_map[sp], label=sp)
        axes[1].scatter(sub.tsne_x, sub.tsne_y, s=10, alpha=0.5, color=color_map[sp], label=sp)
    axes[0].set_title(f'PCA (PC1 {pca.explained_variance_ratio_[0]:.1%}, PC2 {pca.explained_variance_ratio_[1]:.1%})')
    axes[0].set_xlabel('PC1'); axes[0].set_ylabel('PC2')
    axes[1].set_title(f't-SNE (perplexity={perp})')
    axes[1].set_xlabel('t-SNE 1'); axes[1].set_ylabel('t-SNE 2')
    fig.suptitle(f'Simulated germline HCDR3 biochemical feature space, by species (n={len(df)})', fontsize=13)
    axes[1].legend(bbox_to_anchor=(1.02, 1), loc='upper left', fontsize=7, title='Species')
    fig.tight_layout()
    fig.savefig(f'{out_dir}/simulated_HCDR3_dimred_by_species.png', dpi=150, bbox_inches='tight')
    plt.close('all')
    print('Saved plot and coordinate files.')
