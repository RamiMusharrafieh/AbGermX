"""
A small set of real, published de novo antibody designs (RFdiffusion + ProteinMPNN)
to use as test inputs for compare_query_to_species.py.

Source: Bennett, N.R., Watson, J.L., Ragotte, R.J. et al. "Atomically accurate de
novo design of antibodies with RFdiffusion." Nature 649, 183-193 (2026).
https://doi.org/10.1038/s41586-025-09721-5 (preprint: bioRxiv 2024.03.14.585103)

VHH_flu_01: RFdiffusion-generated backbone + ProteinMPNN-designed CDR loops,
targeting the influenza H1 hemagglutinin stem epitope. Cryo-EM confirmed
(PDB 9NH7, EMD-49405) atomic-level agreement between the computational design
and the experimentally determined structure (backbone RMSD 1.45 A, CDR3 RMSD
0.8 A). Binds with 78 nM affinity. Sequence pulled directly from the deposited
PDB entry (https://www.rcsb.org/structure/9NH7); includes the expression
construct's N-terminal cloning residues and C-terminal linker/His-tag, which
are not part of the natural VHH fold -- included as-deposited for traceability.
"""

EXAMPLE_LOOPS = {
    'VHH_flu_01': {
        'sequence': 'MSGQVQLVESGGGLVQPGGSLRLSCAASGKYVNLMSLGWFRQAPGQGLEAVAAISFDGKKTYYADSVKGRFTISRDNSKNTLYLQMNSLRAEDTAVYYCAASVVDSLGVGFYSYWGQGTLVTVSGSGSHHWGSTHHHHHH',
        'target': 'Influenza H1 hemagglutinin (stem epitope)',
        'pdb_id': '9NH7',
        'method': 'RFdiffusion (antibody fine-tune) + ProteinMPNN, cryo-EM validated',
        'affinity': '78 nM (SPR)',
        'source': 'Bennett et al. 2025, Nature 649:183-193, doi:10.1038/s41586-025-09721-5',
    },
}
