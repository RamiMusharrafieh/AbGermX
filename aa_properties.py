"""
Standard amino-acid property definitions used for CDR biophysical annotation.

Hydrophobicity: Kyte & Doolittle (1982) scale, "A simple method for displaying
the hydropathic character of a protein", J Mol Biol 157:105-132. Higher = more
hydrophobic. Reported as the mean over the residues in a region.

Aromatic residues: Phe (F), Trp (W), Tyr (Y)

Polar (uncharged) residues: Ser (S), Thr (T), Asn (N), Gln (Q), Cys (C), Tyr (Y)

KYTE_DOOLITTLE = {
    'A': 1.8, 'R': -4.5, 'N': -3.5, 'D': -3.5, 'C': 2.5,
    'Q': -3.5, 'E': -3.5, 'G': -0.4, 'H': -3.2, 'I': 4.5,
    'L': 3.8, 'K': -3.9, 'M': 1.9, 'F': 2.8, 'P': -1.6,
    'S': -0.8, 'T': -0.7, 'W': -0.9, 'Y': -1.3, 'V': 4.2,
}

AROMATIC_RESIDUES = {'F', 'W', 'Y'}
POLAR_RESIDUES = {'S', 'T', 'N', 'Q', 'C', 'Y'}

def mean_hydrophobicity(seq):
    """Mean Kyte-Doolittle hydrophobicity over a (ungapped) sequence. None if empty
    or if it contains no scorable residues (e.g. all gaps/unknowns)."""
    scores = [KYTE_DOOLITTLE[aa] for aa in seq if aa in KYTE_DOOLITTLE]
    return round(sum(scores) / len(scores), 3) if scores else None

def count_aromatic(seq):
    return sum(1 for aa in seq if aa in AROMATIC_RESIDUES)

def count_polar(seq):
    return sum(1 for aa in seq if aa in POLAR_RESIDUES)

# --- extended properties for dimensionality-reduction feature building ---

ACIDIC_RESIDUES = {'D', 'E'}
BASIC_RESIDUES = {'K', 'R', 'H'}
ALL_20_AA = list('ACDEFGHIKLMNPQRSTVWY')

def count_acidic(seq):
    return sum(1 for aa in seq if aa in ACIDIC_RESIDUES)

def count_basic(seq):
    return sum(1 for aa in seq if aa in BASIC_RESIDUES)

def net_charge(seq):
    return count_basic(seq) - count_acidic(seq)

def aa_composition(seq):
    """Returns a dict of the 20 standard amino acids -> frequency (0-1) in seq.
    Non-standard characters (gaps, X, etc.) are ignored in the denominator."""
    valid = [aa for aa in seq if aa in ALL_20_AA]
    n = len(valid)
    if n == 0:
        return {aa: 0.0 for aa in ALL_20_AA}
    from collections import Counter
    counts = Counter(valid)
    return {aa: counts.get(aa, 0) / n for aa in ALL_20_AA}
