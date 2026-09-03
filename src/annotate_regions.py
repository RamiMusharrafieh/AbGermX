"""
Annotates IMGT-gapped V-REGION amino acid sequences into IMGT framework (FR)
and complementarity-determining region (CDR) segments, using the standard
IMGT unique numbering boundaries (Lefranc et al.):

    FR1   : positions 1-26
    CDR1  : positions 27-38
    FR2   : positions 39-55
    CDR2  : positions 56-65
    FR3   : positions 66-104
    CDR3* : positions 105 onward (partial)

Because these are already IMGT-gapped ('.' = alignment gap at a numbered
position), boundaries can be read directly by column index for the large
majority of sequences. Three IMGT-invariant residues are used as a per-sequence
QC check: 1st-Cys (~pos 23), conserved Trp (~pos 41), 2nd-Cys (~pos 104).
"""
from collections import namedtuple

FR1_END = 26
CDR1_END_NOMINAL = 38
FR2_END_NOMINAL = 55
CDR2_END_NOMINAL = 65
FR3_END_NOMINAL = 104
FR2_WIDTH = FR2_END_NOMINAL - CDR1_END_NOMINAL      # 17, fixed width; never has insertions
FR3_WIDTH = FR3_END_NOMINAL - CDR2_END_NOMINAL      # 39, fixed width; never has insertions
TRP_NOMINAL = 41       # 3rd residue of FR2 in the conserved germline motif
CYS104_NOMINAL = 104   # last residue of FR3

Annotation = namedtuple('Annotation', [
    'fr1', 'cdr1', 'fr2', 'cdr2', 'fr3', 'cdr3_partial',
    'fr1_ungapped', 'cdr1_ungapped', 'fr2_ungapped', 'cdr2_ungapped',
    'fr3_ungapped', 'cdr3_partial_ungapped',
    'anchor_cys23_ok', 'anchor_trp41_ok', 'anchor_cys104_ok', 'qc_pass',
    'cdr1_insertion', 'cdr2_insertion',
])

def _slice(seq, start_1based, end_1based):
    """1-based inclusive slice; returns '' if the sequence doesn't reach that far."""
    if end_1based < start_1based or len(seq) < start_1based:
        return ''
    return seq[start_1based-1:min(end_1based, len(seq))]

def _find_near(seq, target_1based, letter, max_shift=6):
    """Search outward from target_1based (1-based) for the given residue letter,
    preferring the closest match. Returns the 1-based position found, or None."""
    idx0 = target_1based - 1
    if 0 <= idx0 < len(seq) and seq[idx0] == letter:
        return target_1based
    for d in range(1, max_shift+1):
        for cand in (idx0 + d, idx0 - d):
            if 0 <= cand < len(seq) and seq[cand] == letter:
                return cand + 1
    return None

def annotate_sequence(seq):
    """
    seq: IMGT-gapped AA string (with '.' gap characters).
    Boundaries are located per-sequence: FR1 is always the fixed first 26 columns
    (IMGT never places insertions there); the conserved Trp (~pos 41) and 2nd-Cys
    (~pos 104) are then located by local search to absorb any CDR1/CDR2 insertion
    shift before slicing FR2/CDR2/FR3, since insertions there shift everything downstream.
    """
    def ungap(s):
        return s.replace('.', '')

    cys23 = seq[22] if len(seq) > 22 else ''
    anchor_cys23_ok = (cys23 == 'C')

    # --- locate conserved Trp (3rd residue of FR2) to absorb any CDR1 insertion ---
    trp_pos = _find_near(seq, TRP_NOMINAL, 'W', max_shift=6)
    anchor_trp41_ok = (trp_pos == TRP_NOMINAL)
    if trp_pos is not None:
        fr2_start = trp_pos - 2          # Trp is the 3rd column of FR2
        cdr1_insertion = fr2_start - (CDR1_END_NOMINAL + 1)
    else:
        fr2_start = CDR1_END_NOMINAL + 1  # fall back to nominal boundary
        cdr1_insertion = 0
    cdr1_end = fr2_start - 1
    fr2_end = fr2_start + FR2_WIDTH - 1   # FR2 is always exactly 17 columns wide

    # --- locate conserved 2nd-Cys (last residue of FR3) to absorb any CDR2 insertion ---
    cys104_target = fr2_end + (CDR2_END_NOMINAL - FR2_END_NOMINAL) + FR3_WIDTH
    cys104_pos = _find_near(seq, cys104_target, 'C', max_shift=6)
    anchor_cys104_ok = (cys104_pos == CYS104_NOMINAL)
    if cys104_pos is not None:
        fr3_end = cys104_pos
        cdr2_insertion = fr3_end - FR3_WIDTH - (fr2_end + 1) - (CDR2_END_NOMINAL - FR2_END_NOMINAL - 1)
    else:
        fr3_end = cys104_target
        cdr2_insertion = 0
    cdr2_start = fr2_end + 1
    cdr2_end = fr3_end - FR3_WIDTH
    fr3_start = cdr2_end + 1

    fr1 = _slice(seq, 1, FR1_END)
    cdr1 = _slice(seq, FR1_END+1, cdr1_end)
    fr2 = _slice(seq, fr2_start, fr2_end)
    cdr2 = _slice(seq, cdr2_start, cdr2_end)
    fr3 = _slice(seq, fr3_start, fr3_end)
    cdr3_partial = seq[fr3_end:] if len(seq) > fr3_end else ''

    # QC passes when both downstream anchors were located at all, shifted or
    # not: a shifted anchor is an insertion the slicing above already absorbed,
    # whereas a missing one means the sequence is truncated and the FR/CDR
    # boundaries after that point are guesses.
    qc_pass = (trp_pos is not None) and (cys104_pos is not None)

    return Annotation(
        fr1=fr1, cdr1=cdr1, fr2=fr2, cdr2=cdr2, fr3=fr3, cdr3_partial=cdr3_partial,
        fr1_ungapped=ungap(fr1), cdr1_ungapped=ungap(cdr1), fr2_ungapped=ungap(fr2),
        cdr2_ungapped=ungap(cdr2), fr3_ungapped=ungap(fr3),
        cdr3_partial_ungapped=ungap(cdr3_partial),
        anchor_cys23_ok=anchor_cys23_ok, anchor_trp41_ok=anchor_trp41_ok,
        anchor_cys104_ok=anchor_cys104_ok, qc_pass=qc_pass,
        cdr1_insertion=cdr1_insertion, cdr2_insertion=cdr2_insertion,
    )
