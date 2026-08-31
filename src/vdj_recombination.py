"""
Stochastic germline V(D)J recombination simulator for HCDR3 (IGH) generation.

Combines, per species: the germline-encoded V-gene CDR3 stub (from the 2nd-Cys
onward), a D gene (in whatever reading frame emerges from random 5'/3' trimming
-- not a forced "frame 1/2/3" choice, since real frame is determined by total
nucleotide count, not gene identity), and a J gene's CDR3-contributing prefix
(up to the conserved FR4 WG.G motif) -- with random exonuclease trimming at each
junction and random N-nucleotide insertion (TdT-like) between segments.

This models the COMBINATORIAL + rough JUNCTIONAL diversity ceiling of germline
HCDR3 generation. It is a simplified illustrative model, not a validated
biophysical simulation of real recombination -- see caveats in generate_hcdr3s().
"""
import re
import random
from collections import defaultdict
from imgt_parser import norm_functionality
from annotate_regions import annotate_sequence

FR4_MOTIF = re.compile(r'TGG(?:GG[ACGT]|[ACGT]G[ACGT])GG[ACGT]GG[ACGT]', re.IGNORECASE)  # not used; AA-level motif used instead
FR4_AA_MOTIF = re.compile(r'WG.G')

CODON_TABLE = {
    'TTT':'F','TTC':'F','TTA':'L','TTG':'L','CTT':'L','CTC':'L','CTA':'L','CTG':'L',
    'ATT':'I','ATC':'I','ATA':'I','ATG':'M','GTT':'V','GTC':'V','GTA':'V','GTG':'V',
    'TCT':'S','TCC':'S','TCA':'S','TCG':'S','CCT':'P','CCC':'P','CCA':'P','CCG':'P',
    'ACT':'T','ACC':'T','ACA':'T','ACG':'T','GCT':'A','GCC':'A','GCA':'A','GCG':'A',
    'TAT':'Y','TAC':'Y','TAA':'*','TAG':'*','CAT':'H','CAC':'H','CAA':'Q','CAG':'Q',
    'AAT':'N','AAC':'N','AAA':'K','AAG':'K','GAT':'D','GAC':'D','GAA':'E','GAG':'E',
    'TGT':'C','TGC':'C','TGA':'*','TGG':'W','CGT':'R','CGC':'R','CGA':'R','CGG':'R',
    'AGT':'S','AGC':'S','AGA':'R','AGG':'R','GGT':'G','GGC':'G','GGA':'G','GGG':'G',
}

def translate_dna(nt_seq):
    """Translate from position 0, in-frame, stopping at (and excluding) the first
    stop codon. Returns (aa_string, hit_stop_codon: bool)."""
    aa = []
    nt_seq = nt_seq.upper()
    for i in range(0, len(nt_seq) - 2, 3):
        codon = nt_seq[i:i+3]
        residue = CODON_TABLE.get(codon, 'X')
        if residue == '*':
            return ''.join(aa), True
        aa.append(residue)
    return ''.join(aa), False

def build_species_vdj_pools(aa_records, nt_records, species, locus='IGH'):
    """
    Returns dict(v_tails, d_segments, j_prefixes) of (gene_name, nt_sequence) lists
    -- the germline nucleotide material each gene type contributes to a CDR3 junction.
    """
    aa_by_key = {(r['gene_allele'], r['species']): r for r in aa_records if r['locus']==locus}
    nt_by_key = {(r['gene_allele'], r['species']): r for r in nt_records if r['locus']==locus}

    def functional_genes(region):
        seen_genes = {}
        for r in aa_records:
            if r['locus']==locus and r['species']==species and r['region']==region \
               and norm_functionality(r['functionality'])=='Functional':
                if r['allele']=='01' or r['gene'] not in seen_genes:
                    seen_genes[r['gene']] = r
        return list(seen_genes.values())

    v_tails = []
    for r in functional_genes('V-REGION'):
        nt_r = nt_by_key.get((r['gene_allele'], r['species']))
        if nt_r is None:
            continue
        ann = annotate_sequence(r['seq'])
        if not ann.fr3:
            continue
        fr3_end_aa_col = len(ann.fr1) + len(ann.cdr1) + len(ann.fr2) + len(ann.cdr2) + len(ann.fr3)
        # fr3_end_aa_col counts *gapped* columns consumed so far only if slices came from
        # the gapped string in order -- recompute directly from the gapped seq length instead:
        gapped = r['seq']
        # locate fr3 end by finding fr3 substring's end position in the gapped AA seq
        idx = gapped.find(ann.fr3)
        if idx == -1 or not ann.fr3:
            continue
        fr3_end_aa_col = idx + len(ann.fr3)
        nt_tail = nt_r['seq'][fr3_end_aa_col*3:].replace('.', '')
        if nt_tail:
            v_tails.append((r['gene'], nt_tail))

    d_segments = []
    for r in functional_genes('D-REGION'):
        nt_r = nt_by_key.get((r['gene_allele'], r['species']))
        if nt_r is None:
            continue
        d_nt = nt_r['seq'].replace('.', '')
        if d_nt:
            d_segments.append((r['gene'], d_nt))

    j_prefixes = []
    for r in functional_genes('J-REGION'):
        nt_r = nt_by_key.get((r['gene_allele'], r['species']))
        if nt_r is None:
            continue
        m = FR4_AA_MOTIF.search(r['seq'])
        if not m:
            continue
        j_prefix_nt = nt_r['seq'][:m.start()*3].replace('.', '')
        if j_prefix_nt:
            j_prefixes.append((r['gene'], j_prefix_nt))

    return dict(v_tails=v_tails, d_segments=d_segments, j_prefixes=j_prefixes)

def simulate_one_junction(pools, rng, max_v_trim=6, max_d_trim=8, max_j_trim=6, max_n_insert=12):
    """One stochastic VDJ recombination draw. Returns dict with the outcome, or
    None if no productive (stop-codon-free, non-empty) peptide resulted."""
    if not (pools['v_tails'] and pools['d_segments'] and pools['j_prefixes']):
        return None
    v_gene, v_nt = rng.choice(pools['v_tails'])
    d_gene, d_nt = rng.choice(pools['d_segments'])
    j_gene, j_nt = rng.choice(pools['j_prefixes'])

    v_trim = min(rng.randint(0, max_v_trim), max(0, len(v_nt)-1))
    v_kept = v_nt[:len(v_nt)-v_trim] if v_trim else v_nt

    d_trim5 = rng.randint(0, min(max_d_trim, len(d_nt)))
    d_trim3 = rng.randint(0, min(max_d_trim, max(0, len(d_nt)-d_trim5)))
    d_kept = d_nt[d_trim5: len(d_nt)-d_trim3] if len(d_nt) - d_trim5 - d_trim3 > 0 else ''

    j_trim = min(rng.randint(0, max_j_trim), max(0, len(j_nt)-1))
    j_kept = j_nt[j_trim:] if j_trim else j_nt

    n1 = ''.join(rng.choice('ACGT') for _ in range(rng.randint(0, max_n_insert)))
    n2 = ''.join(rng.choice('ACGT') for _ in range(rng.randint(0, max_n_insert)))

    junction_nt = v_kept + n1 + d_kept + n2 + j_kept
    aa, hit_stop = translate_dna(junction_nt)
    if hit_stop or len(aa) < 3:
        return None
    return dict(v_gene=v_gene, d_gene=d_gene, j_gene=j_gene, junction_nt=junction_nt,
                cdr3_aa=aa, v_trim=v_trim, d_trim5=d_trim5, d_trim3=d_trim3, j_trim=j_trim,
                n1_len=len(n1), n2_len=len(n2))

def generate_hcdr3s(pools, n_target, rng, max_attempts_factor=6, **kwargs):
    """
    Draws stochastic recombination events until n_target productive HCDR3s are
    collected (or max_attempts is exhausted). Returns list of result dicts.

    CAVEATS (illustrative model, not a calibrated repertoire simulator):
    - Trim-length and N-insertion-length distributions are uniform random over a
      fixed range, not fit to any real empirical distribution (real trimming/TdT
      insertion vary by species, locus, and even individual, and are known to be
      biased rather than uniform).
    - Gene usage is uniform random across functional V/D/J genes; real repertoires
      show strong, unequal gene usage biases (some V/D/J genes are used far more
      than others).
    - No modeling of P-nucleotides (palindromic additions at trimmed ends).
    - No selection beyond "no premature stop codon" -- real repertoires undergo
      further selection (structural viability, self-tolerance, antigen selection)
      not modeled here.
    """
    results = []
    attempts = 0
    max_attempts = n_target * max_attempts_factor
    while len(results) < n_target and attempts < max_attempts:
        r = simulate_one_junction(pools, rng, **kwargs)
        attempts += 1
        if r is not None:
            results.append(r)
    return results, attempts
