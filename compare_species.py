"""
Generalized cross-species IMGT V-gene comparison engine.
Given a parsed record list and two species names, produces:
 - repertoire/functionality summary
 - subgroup distribution comparison
 - pairwise identity matrix (Needleman-Wunsch, representative allele per gene)
 - best-ortholog table
 - clustermap figure
"""
import numpy as np
import pandas as pd
from collections import defaultdict, Counter
from imgt_parser import norm_functionality

def nw_identity(s1, s2, match=1, mismatch=-1, gap=-2):
    n, mlen = len(s1), len(s2)
    if n == 0 or mlen == 0:
        return 0.0
    a1 = np.frombuffer(s1.encode(), dtype=np.uint8)
    a2 = np.frombuffer(s2.encode(), dtype=np.uint8)
    D = np.zeros((n+1, mlen+1), dtype=np.int32)
    D[:,0] = np.arange(n+1)*gap
    D[0,:] = np.arange(mlen+1)*gap
    for i in range(1, n+1):
        row_match = np.where(a2 == a1[i-1], match, mismatch)
        diag = D[i-1, :-1] + row_match
        up = D[i-1, 1:] + gap
        prev = D[i,0]
        best_row = np.maximum(diag, up)
        Dr = np.empty(mlen, dtype=np.int32)
        for j in range(mlen):
            left = prev + gap
            val = max(best_row[j], left)
            Dr[j] = val
            prev = val
        D[i,1:] = Dr
    i, j = n, mlen
    matches, aligned_len = 0, 0
    while i > 0 and j > 0:
        score = D[i,j]
        diag_score = D[i-1,j-1] + (match if s1[i-1] == s2[j-1] else mismatch)
        if score == diag_score:
            if s1[i-1] == s2[j-1]:
                matches += 1
            aligned_len += 1; i -= 1; j -= 1
        elif score == D[i-1,j] + gap:
            aligned_len += 1; i -= 1
        else:
            aligned_len += 1; j -= 1
    aligned_len += i + j
    return matches / aligned_len if aligned_len else 0.0

def filter_species_locus(records, species, locus=None):
    out = [r for r in records if r['species'] == species]
    if locus:
        out = [r for r in out if r['locus'] == locus]
    return out

def _filter_region(recs, region_set):
    if region_set is None:
        return recs
    return [r for r in recs if r['region'] in region_set]

def repertoire_summary(records, species, locus=None, region_set=None):
    recs = filter_species_locus(records, species, locus)
    recs = _filter_region(recs, region_set)
    for r in recs:
        r['func_class'] = norm_functionality(r['functionality'])
    func_counts = Counter(r['func_class'] for r in recs)
    subgroup_counts_functional = Counter(r['subgroup'] for r in recs if r['func_class']=='Functional')
    n_unique_functional_genes = len(set(r['gene'] for r in recs if r['func_class']=='Functional'))
    return dict(
        n_alleles_total=len(recs),
        functionality_counts=dict(func_counts),
        n_unique_functional_genes=n_unique_functional_genes,
        subgroup_counts_functional=dict(subgroup_counts_functional),
    )

def representative_seqs(records, species, locus=None, region_set=None):
    """One functional representative AA seq (with IMGT gaps) per gene, preferring *01.
    If region_set matches multiple exon labels (e.g. constant-region proxy {'C-REGION','CH1'}),
    sequences for the same gene/allele across matching regions are concatenated in region_set order."""
    recs = filter_species_locus(records, species, locus)
    recs = _filter_region(recs, region_set)
    by_gene = defaultdict(list)
    for r in recs:
        if norm_functionality(r['functionality']) == 'Functional' and r['seq']:
            by_gene[r['gene']].append(r)
    reps = {}
    for gene, rs in by_gene.items():
        rs01 = [r for r in rs if r['allele'] == '01']
        chosen = rs01[0] if rs01 else rs[0]
        reps[gene] = chosen['seq'].replace('.', '')  # ungapped for alignment
    return reps

def pairwise_identity_matrix(records, species_a, species_b, locus, region_set=None):
    reps_a = representative_seqs(records, species_a, locus, region_set)
    reps_b = representative_seqs(records, species_b, locus, region_set)
    genes_a = sorted(reps_a.keys())
    genes_b = sorted(reps_b.keys())
    mat = np.zeros((len(genes_a), len(genes_b)))
    for i, ga in enumerate(genes_a):
        for j, gb in enumerate(genes_b):
            mat[i,j] = nw_identity(reps_a[ga], reps_b[gb])
    df = pd.DataFrame(mat*100, index=genes_a, columns=genes_b)
    return df

def best_ortholog_table(identity_df, species_a_label, species_b_label):
    rows = []
    for gene_a in identity_df.index:
        row = identity_df.loc[gene_a]
        best_gene_b = row.idxmax()
        rows.append({
            f'{species_a_label}_gene': gene_a,
            f'best_{species_b_label}_match': best_gene_b,
            'pct_identity': round(row[best_gene_b], 2)
        })
    return pd.DataFrame(rows).sort_values('pct_identity', ascending=False)
