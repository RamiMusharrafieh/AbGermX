"""
Compares a query antibody loop sequence (e.g. an RFdiffusion+ProteinMPNN designed
VHH or CDR sequence) against the natural cross-species IGHV CDR1/CDR2 repertoire,
both biochemically (same feature space as run_dimred_cdr_only.py) and by direct
sequence alignment.

only CDR1 and CDR2 are compared. 
"""
import re
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
from annotate_regions import annotate_sequence
from aa_properties import mean_hydrophobicity, count_aromatic, count_polar, count_acidic, count_basic, net_charge, aa_composition, ALL_20_AA
from compare_species import nw_identity

FR4_MOTIF = re.compile(r'WG.G')  # conserved Trp-Gly-x-Gly at the CDR3->FR4 boundary

def extract_query_cdrs(full_sequence):
    """Extract CDR1/CDR2 (and, if present, full CDR3) from a full VH/VHH sequence
    using the same anchor-residue logic used for germline genes. Works on plain
    (ungapped) sequences directly. the anchor search tolerates the natural
    absence of IMGT gap characters in a real expressed sequence."""
    ann = annotate_sequence(full_sequence)
    cdr3_full = ''
    if ann.fr3: 
        
        fr3_end_idx = full_sequence.find(ann.fr3) + len(ann.fr3) if ann.fr3 in full_sequence else None
        if fr3_end_idx is not None:
            m = FR4_MOTIF.search(full_sequence, fr3_end_idx)
            if m:
                cdr3_full = full_sequence[fr3_end_idx:m.start()]
    return dict(CDR1=ann.cdr1_ungapped, CDR2=ann.cdr2_ungapped, CDR3_full=cdr3_full,
                qc_pass=ann.qc_pass, cdr1_insertion=ann.cdr1_insertion, cdr2_insertion=ann.cdr2_insertion)

def _feature_vector(cdr1, cdr2):
    """Same 34-feature CDR1/CDR2 schema as run_dimred_cdr_only.py, for one sequence."""
    out = {}
    for prefix, seq in [('CDR1', cdr1), ('CDR2', cdr2)]:
        out[f'{prefix}_length'] = len(seq)
        out[f'{prefix}_hydrophobicity'] = mean_hydrophobicity(seq) or 0.0
        out[f'{prefix}_aromatic_count'] = count_aromatic(seq)
        out[f'{prefix}_polar_count'] = count_polar(seq)
        out[f'{prefix}_acidic_count'] = count_acidic(seq)
        out[f'{prefix}_basic_count'] = count_basic(seq)
        out[f'{prefix}_net_charge'] = net_charge(seq)
    comp = aa_composition(cdr1 + cdr2)
    for aa in ALL_20_AA:
        out[f'CDR_comp_{aa}'] = comp[aa]
    return out

def compare_query(query_name, query_cdr1, query_cdr2, locus='IGH', out_dir='/mnt/user-data/outputs',
                   bulk_path='/mnt/user-data/uploads/IMGTGENEDB-ReferenceSequences.fasta-AA-WithGaps-F_ORF_inframeP'):
    records = parse_bulk_fasta(bulk_path)
    ann = build_table(records, locus=locus, species_list=None, functional_only=True, allele_scope='representative')
    ann = ann[ann.qc_pass].reset_index(drop=True)

    # --- biochemical feature space (natural genes) ---
    feat_rows = [_feature_vector(r.CDR1, r.CDR2) for r in ann.itertuples()]
    feat_df = pd.DataFrame(feat_rows)
    scaler = StandardScaler().fit(feat_df.values)
    X = scaler.transform(feat_df.values)

    query_feat = pd.DataFrame([_feature_vector(query_cdr1, query_cdr2)])[feat_df.columns]
    Xq = scaler.transform(query_feat.values)

    # per-gene Euclidean distance from the query, in standardized feature space
    dists = np.linalg.norm(X - Xq, axis=1)
    ann = ann.copy()
    ann['feature_distance_to_query'] = dists

    # --- sequence identity (NW alignment) ---
    ann['CDR1_identity_to_query'] = [nw_identity(query_cdr1, c) if query_cdr1 and c else 0.0 for c in ann.CDR1]
    ann['CDR2_identity_to_query'] = [nw_identity(query_cdr2, c) if query_cdr2 and c else 0.0 for c in ann.CDR2]
    ann['CDR1_CDR2_identity_to_query'] = [
        nw_identity(query_cdr1 + query_cdr2, c1 + c2) if (query_cdr1 or query_cdr2) else 0.0
        for c1, c2 in zip(ann.CDR1, ann.CDR2)
    ]

    # --- per-species summary ranking ---
    species_summary = ann.groupby('species').agg(
        n_genes=('gene', 'size'),
        mean_feature_distance=('feature_distance_to_query', 'mean'),
        min_feature_distance=('feature_distance_to_query', 'min'),
        mean_seq_identity=('CDR1_CDR2_identity_to_query', 'mean'),
        max_seq_identity=('CDR1_CDR2_identity_to_query', 'max'),
    ).reset_index()
    species_summary['closest_gene'] = species_summary['species'].map(
        lambda sp: ann.loc[ann[ann.species == sp].feature_distance_to_query.idxmin(), 'gene_allele']
    )
    species_summary = species_summary.sort_values('mean_feature_distance')

    per_gene_out = ann[['species', 'gene_allele', 'CDR1', 'CDR2', 'feature_distance_to_query',
                         'CDR1_identity_to_query', 'CDR2_identity_to_query', 'CDR1_CDR2_identity_to_query']]
    per_gene_out = per_gene_out.sort_values('feature_distance_to_query')

    species_summary.to_csv(f'{out_dir}/{query_name}_species_ranking.csv', index=False)
    per_gene_out.to_csv(f'{out_dir}/{query_name}_per_gene_comparison.csv', index=False)

    # --- PCA projection (query overlaid on natural CDR feature space) ---
    pca = PCA(n_components=2, random_state=42)
    pca_coords = pca.fit_transform(X)
    query_pca = pca.transform(Xq)

    fig, ax = plt.subplots(figsize=(11, 9))
    species_list = sorted(ann.species.unique())
    cmap = plt.get_cmap('tab20' if len(species_list) <= 20 else 'gist_ncar')
    color_map = {sp: cmap(i / max(1, len(species_list)-1)) for i, sp in enumerate(species_list)}
    for sp in species_list:
        mask = (ann.species == sp).values
        ax.scatter(pca_coords[mask, 0], pca_coords[mask, 1], s=14, alpha=0.55, color=color_map[sp], label=sp)
    ax.scatter(query_pca[:, 0], query_pca[:, 1], s=260, marker='*', color='black',
               edgecolor='white', linewidth=1.2, zorder=10, label=f'{query_name} (query)')
    ax.annotate(query_name, (query_pca[0, 0], query_pca[0, 1]), textcoords='offset points',
                xytext=(10, 8), fontsize=10, fontweight='bold')
    ax.set_xlabel(f'PC1 ({pca.explained_variance_ratio_[0]:.1%})')
    ax.set_ylabel(f'PC2 ({pca.explained_variance_ratio_[1]:.1%})')
    ax.set_title(f'{query_name} projected onto natural {locus}V CDR1/CDR2 biochemical feature space')
    ax.legend(bbox_to_anchor=(1.02, 1), loc='upper left', fontsize=6.5, title='Species')
    fig.tight_layout()
    fig.savefig(f'{out_dir}/{query_name}_PCA_projection.png', dpi=150, bbox_inches='tight')
    plt.close('all')

    # --- t-SNE with query merged in (no out-of-sample projection for t-SNE, so refit including it) ---
    X_with_query = np.vstack([X, Xq])
    perp = min(30, max(5, len(ann)//4))
    tsne = TSNE(n_components=2, perplexity=perp, random_state=42, init='pca', learning_rate='auto')
    tsne_coords = tsne.fit_transform(X_with_query)
    tsne_natural, tsne_query = tsne_coords[:-1], tsne_coords[-1:]

    fig, ax = plt.subplots(figsize=(11, 9))
    for sp in species_list:
        mask = (ann.species == sp).values
        ax.scatter(tsne_natural[mask, 0], tsne_natural[mask, 1], s=14, alpha=0.55, color=color_map[sp], label=sp)
    ax.scatter(tsne_query[:, 0], tsne_query[:, 1], s=260, marker='*', color='black',
               edgecolor='white', linewidth=1.2, zorder=10, label=f'{query_name} (query)')
    ax.annotate(query_name, (tsne_query[0, 0], tsne_query[0, 1]), textcoords='offset points',
                xytext=(10, 8), fontsize=10, fontweight='bold')
    ax.set_xlabel('t-SNE 1'); ax.set_ylabel('t-SNE 2')
    ax.set_title(f'{query_name} in natural {locus}V CDR1/CDR2 feature space (t-SNE, query merged into fit)')
    ax.legend(bbox_to_anchor=(1.02, 1), loc='upper left', fontsize=6.5, title='Species')
    fig.tight_layout()
    fig.savefig(f'{out_dir}/{query_name}_tSNE_projection.png', dpi=150, bbox_inches='tight')
    plt.close('all')

    return species_summary, per_gene_out
