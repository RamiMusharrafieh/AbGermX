#!/usr/bin/env python3
"""
Compares a query antibody loop sequence (e.g. an RFdiffusion + ProteinMPNN
designed VHH, or any VH/VHH) against the natural cross-species IGHV CDR1/CDR2
repertoire, both biochemically (the same feature space as
run_dimred_cdr_only.py) and by direct sequence alignment.

Scope note: only CDR1 and CDR2 are compared here. Use
compare_query_cdr3_to_species.py for the CDR3 loop.
"""
import argparse
import re

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
from config import DEFAULT_AA_BULK, ensure_out_dir
from imgt_parser import parse_bulk_fasta
from build_annotation_table import build_table
from annotate_regions import annotate_sequence
from aa_properties import (mean_hydrophobicity, count_aromatic, count_polar,
                           count_acidic, count_basic, net_charge, aa_composition, ALL_20_AA)
from compare_species import nw_identity

FR4_MOTIF = re.compile(r'WG.G')  # conserved Trp-Gly-x-Gly at the CDR3->FR4 boundary
N_HIGHLIGHT = 3


def extract_query_cdrs(full_sequence):
    """Extract CDR1/CDR2 (and, if present, the full CDR3) from a VH/VHH sequence
    using the same anchor-residue logic as the germline annotator. Works on plain
    ungapped sequences: the anchor search tolerates the absence of IMGT gap
    characters in a real expressed sequence."""
    ann = annotate_sequence(full_sequence)
    cdr3_full = ''
    if ann.fr3 and ann.fr3 in full_sequence:
        fr3_end_idx = full_sequence.find(ann.fr3) + len(ann.fr3)
        m = FR4_MOTIF.search(full_sequence, fr3_end_idx)
        if m:
            cdr3_full = full_sequence[fr3_end_idx:m.start()]
    return dict(CDR1=ann.cdr1_ungapped, CDR2=ann.cdr2_ungapped, CDR3_full=cdr3_full,
                qc_pass=ann.qc_pass, cdr1_insertion=ann.cdr1_insertion,
                cdr2_insertion=ann.cdr2_insertion)


def _feature_vector(cdr1, cdr2):
    """The same CDR1/CDR2 feature schema as run_dimred_cdr_only.py, for one sequence."""
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


def _projection_figure(coords, query_coord, species, highlight, query_name,
                       title, subtitle, xlabel, ylabel, out_path, point_size=16):
    plot_style.use_style()
    fig, ax = plt.subplots(figsize=(9.6, 8.2))
    plot_style.query_projection(ax, coords, species, query_coord, highlight,
                                query_name, point_size=point_size)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    handles, labels = ax.get_legend_handles_labels()
    fig.subplots_adjust(left=0.085, right=0.975, top=0.855, bottom=0.20)
    plot_style.figure_title(fig, title, subtitle, x=0.085, y=0.945)
    plot_style.species_legend(fig, handles, labels, ncol=2, y=0.135,
                              title='Closest species by biochemical distance',
                              italic=set(highlight))
    fig.savefig(out_path)
    plt.close(fig)


def compare_query(query_name, query_cdr1, query_cdr2, locus='IGH', out_dir=None,
                  bulk_path=DEFAULT_AA_BULK):
    out_dir = ensure_out_dir(out_dir)
    records = parse_bulk_fasta(bulk_path)
    ann = build_table(records, locus=locus, species_list=None, functional_only=True,
                      allele_scope='representative')
    ann = ann[ann.qc_pass].reset_index(drop=True)

    # --- biochemical feature space (natural genes) ---
    feat_df = pd.DataFrame([_feature_vector(r.CDR1, r.CDR2) for r in ann.itertuples()])
    scaler = StandardScaler().fit(feat_df.values)
    X = scaler.transform(feat_df.values)

    query_feat = pd.DataFrame([_feature_vector(query_cdr1, query_cdr2)])[feat_df.columns]
    Xq = scaler.transform(query_feat.values)

    # per-gene Euclidean distance from the query, in standardized feature space
    ann = ann.copy()
    ann['feature_distance_to_query'] = np.linalg.norm(X - Xq, axis=1)

    # --- sequence identity (NW alignment) ---
    ann['CDR1_identity_to_query'] = [nw_identity(query_cdr1, c) if query_cdr1 and c else 0.0
                                     for c in ann.CDR1]
    ann['CDR2_identity_to_query'] = [nw_identity(query_cdr2, c) if query_cdr2 and c else 0.0
                                     for c in ann.CDR2]
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
                        'CDR1_identity_to_query', 'CDR2_identity_to_query',
                        'CDR1_CDR2_identity_to_query']].sort_values('feature_distance_to_query')

    species_summary.to_csv(f'{out_dir}/{query_name}_species_ranking.csv', index=False)
    per_gene_out.to_csv(f'{out_dir}/{query_name}_per_gene_comparison.csv', index=False)

    highlight = species_summary.species.head(N_HIGHLIGHT).tolist()
    subtitle = (f'{len(ann)} germline genes across {ann.species.nunique()} species; '
                f'the {N_HIGHLIGHT} species closest to the query are picked out')

    # --- PCA projection (true out-of-sample transform of the query) ---
    pca = PCA(n_components=2, random_state=42)
    pca_coords = pca.fit_transform(X)
    _projection_figure(
        pca_coords, pca.transform(Xq), ann.species, highlight, query_name,
        title=f'{query_name} against natural {locus}V CDR1/CDR2 space',
        subtitle=subtitle,
        xlabel=f'PC1 ({pca.explained_variance_ratio_[0]:.1%} of variance)',
        ylabel=f'PC2 ({pca.explained_variance_ratio_[1]:.1%} of variance)',
        out_path=f'{out_dir}/{query_name}_PCA_projection.png',
    )

    # --- t-SNE: no out-of-sample extension, so the query joins the fit ---
    perp = min(30, max(5, len(ann) // 4))
    tsne = TSNE(n_components=2, perplexity=perp, random_state=42, init='pca',
                learning_rate='auto')
    tsne_coords = tsne.fit_transform(np.vstack([X, Xq]))
    _projection_figure(
        tsne_coords[:-1], tsne_coords[-1:], ann.species, highlight, query_name,
        title=f'{query_name} against natural {locus}V CDR1/CDR2 space',
        subtitle=f'{subtitle}. t-SNE (perplexity {perp}) with the query included in the fit',
        xlabel='t-SNE 1', ylabel='t-SNE 2',
        out_path=f'{out_dir}/{query_name}_tSNE_projection.png',
    )

    return species_summary, per_gene_out


def main():
    from example_designed_loops import EXAMPLE_LOOPS
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--name', default='VHH_flu_01', help='Label for the query, used in filenames')
    p.add_argument('--sequence', default=None,
                   help='Full VH/VHH sequence (default: the bundled VHH_flu_01 example)')
    p.add_argument('--locus', default='IGH', help='Locus to compare against (default: IGH)')
    p.add_argument('--bulk_path', default=DEFAULT_AA_BULK, help='IMGT/GENE-DB AA bulk FASTA')
    p.add_argument('--out_dir', default=None, help='Output directory')
    args = p.parse_args()

    sequence = args.sequence or EXAMPLE_LOOPS[args.name]['sequence']
    cdrs = extract_query_cdrs(sequence)
    print(f"Query CDR1: {cdrs['CDR1']}\nQuery CDR2: {cdrs['CDR2']}")
    summary, _ = compare_query(args.name, cdrs['CDR1'], cdrs['CDR2'],
                               locus=args.locus, out_dir=args.out_dir, bulk_path=args.bulk_path)
    print(summary.head(10).to_string(index=False))


if __name__ == '__main__':
    main()
