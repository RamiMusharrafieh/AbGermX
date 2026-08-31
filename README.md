# IMGT Antibody Germline Cross-Species Comparison Toolkit

Compares immunoglobulin (IG) germline genes between species in an
IMGT/GENE-DB bulk FASTA export. **TCR loci (TRA/TRB/TRG/TRD) are excluded
entirely** — only antibody genes (IGH, IGK, IGL) are loaded and analyzed.

This toolkit covers several related analyses:

1. **Pairwise species/locus/gene-type comparison** — identity clustermaps, best-ortholog tables (`workflow.py`)
2. **FR/CDR annotation** — splits V genes into framework/CDR regions with biophysical properties (`annotate_regions.py`, `build_annotation_table.py`)
3. **Cross-species clustering (PCA/t-SNE)** — visualizes how species group by CDR/framework biochemistry (`run_dimred.py`)
4. **Query sequence comparison** — projects an arbitrary antibody loop (e.g. a computationally designed binder) onto the natural cross-species feature space (`compare_query_to_species.py`)
5. **Stochastic HCDR3 simulation** — simulates V(D)J recombination per species from germline genes to explore CDR3 diversity (`vdj_recombination.py`, `run_vdj_simulation.py`)

See the [`examples/`](examples/) folder for sample output images.

## Quickstart

```bash
git clone <this-repo-url>
cd <repo-directory>
pip install -r requirements.txt

# Download IMGT bulk data yourself -- see data/README.md (not bundled here)

python3 workflow.py --query "macaque vs mouse, heavy chain, V genes" \
    --bulk_path /path/to/IMGTGENEDB-ReferenceSequences.fasta-AA-WithGaps-F+ORF+inframeP
```

```bash
# Freeform query
python3 workflow.py --query "macaque vs mouse, heavy chain, V genes"

# Structured
python3 workflow.py --species_a "macaque" --species_b "mouse" --chain heavy --gene V

# See what species/loci are available in your bulk file
python3 list_species.py
```

## Options

**Species** — common name (human, mouse, macaque, rabbit, cow, ...) or IMGT
scientific name, exact/partial matched against what's actually in your
file. Unresolvable names print the full species list so you can pick a valid
one. Some common names are ambiguous (e.g. "macaque" = several real species);
see `species_aliases.py` to see/adjust the defaults.

**Chain** — `heavy` (IGH) | `kappa` (IGK) | `lambda` (IGL) | `light` (runs
kappa AND lambda as two separate outputs)

**Gene type** — `V` | `D` | `J` | `C`. V/D/J use IMGT's single-exon region
labels directly. `C` (constant region) is straightforward for kappa/lambda
(single C-REGION exon) but heavy-chain constant genes span multiple exons
(CH1, hinge, CH2, CH3...). as a simplification, heavy-chain `C` comparisons
use the CH1 exon as a proxy rather than the full multi-domain constant region.

## What it produces (per locus, in the output directory)

- `<a>_vs_<b>_<locus>_<gene>_identity_matrix.csv` — full pairwise %identity
  matrix (Needleman-Wunsch global alignment, one representative functional
  allele per gene)
- `<a>_vs_<b>_<locus>_<gene>_best_orthologs.csv` — each species-A gene's best
  species-B match and %identity
- `<a>_vs_<b>_<locus>_<gene>_clustermap.png` — hierarchically clustered
  identity heatmap

Console output also prints repertoire size and functionality (Functional /
ORF / Pseudogene) counts for both species.

## Files

| File | Purpose |
|---|---|
| `imgt_parser.py` | Parses IMGT/GENE-DB bulk FASTA; restricts to IG loci |
| `species_aliases.py` | Common-name → IMGT species-name resolution (+ match) |
| `query_resolver.py` | Resolves chain/gene-type keywords and freeform queries |
| `compare_species.py` | Alignment engine: NW identity, repertoire summaries, ortholog tables |
| `run_comparison.py` | Orchestrates one species-pair/locus comparison + figure |
| `workflow.py` | Top-level CLI (freeform or structured input) |
| `list_species.py` | Prints species/locus availability in a bulk file |

## Known limitations

- Runtime scales with (genes in species A) × (genes in species B); a
  126×335-gene IGH V-gene comparison takes ~2-3 minutes. D/J genes are much
  faster (short sequences). Very large pairs (e.g. two 300+ gene species) may
  take several minutes.
- Species/gene coverage in IMGT/GENE-DB varies a lot and some species only have
  a handful of alleles for a given locus, or none at all for a given chain;
  the tool reports this rather than failing silently.
- Global (not local) alignment is used for identity scoring, appropriate for
  same-length germline domains; very short D genes can have noisy %identity
  simply because there's little sequence to align.

---

## FR/CDR annotation workflow

`annotate_regions.py` + `build_annotation_table.py` + `aa_properties.py` produce
a per-gene CSV that splits each IGHV sequence into its framework (FR) and
complementarity-determining (CDR) regions, plus biophysical summary columns for
CDR1/CDR2 (length, hydrophobicity, aromatic and polar residue counts).

### Region boundary strategy

Boundaries follow the **IMGT unique numbering scheme** for the V-DOMAIN
(Lefranc et al., "IMGT unique numbering for immunoglobulin and T cell receptor
variable domains and Ig superfamily V-like domains", *Dev Comp Immunol*
27(1):55-77, 2003. https://doi.org/10.1016/S0145-305X(02)00039-3):

| Region | IMGT positions |
|---|---|
| FR1 | 1–26 |
| CDR1 | 27–38 |
| FR2 | 39–55 |
| CDR2 | 56–65 |
| FR3 | 66–104 |
| CDR3 (partial)* | 105 onward |

\*Only the germline-V-encoded start of CDR3 is present — the full CDR3 loop
also depends on D and J gene segments and junctional (N/P nucleotide) addition,
none of which exist in a germline V-REGION reference sequence alone.

Because the input sequences are already **IMGT-gapped** (`.` = alignment gap
at a specific numbered position), boundary columns can mostly be read directly.
The complication: FR1/FR2/FR3 are fixed-width and never contain insertions,
but CDR1 and CDR2 can vary in length between genes via IMGT's insertion-position
codes — confirmed empirically in this dataset (several species, e.g. macaque,
orangutan, rat, trout, salmon, carry a consistent +1 to +4 residue CDR1
insertion relative to the human/mouse template). A single fixed column cutoff
handles ~49% of sequences correctly across all species; the rest are silently
misaligned by the insertion offset.

**Fix used here:** rather than trusting a fixed column index, three IMGT-invariant
anchor residues are located empirically per sequence (1st-Cys ~pos 23,
conserved Trp ~pos 41, 2nd-Cys ~pos 104), searched in a small local window
around their nominal position. The Trp anchors the CDR1→FR2 boundary and the
2nd-Cys anchors the CDR2→FR3 boundary, absorbing any insertion shift before
slicing. This raised the overall QC pass rate to 96.1% across all 24 species
in the reference file; the remaining failures are genuinely partial/truncated
sequences (confirmed by their much shorter length), not annotation errors.

Output columns `anchor_cys23_ok`, `anchor_trp41_ok`, `anchor_cys104_ok`, and
`qc_pass` flag exactly which entries to treat cautiously. `cdr1_insertion` and
`cdr2_insertion` report the detected insertion size (in residues) for that gene.

### Biophysical property definitions

**Hydrophobicity** — Kyte & Doolittle (1982) scale: "A simple method for
displaying the hydropathic character of a protein", *J Mol Biol* 157(1):105-132.
https://doi.org/10.1016/0022-2836(82)90515-0

The scale assigns each amino acid a value from **+4.5 (Ile, most hydrophobic)
to −4.5 (Arg, most hydrophilic)**, derived from water-vapor transfer free
energies and interior/exterior side-chain distributions in known protein
structures. `CDR1_hydrophobicity` / `CDR2_hydrophobicity` /
`CDR1_CDR2_hydrophobicity` report the **mean** score across the residues in
that region — positive values indicate a more hydrophobic loop on average,
negative values a more hydrophilic one.

**Aromatic residues** — Phe (F), Trp (W), Tyr (Y) Aromatic residues in CDR
loops are frequently implicated in antigen contact interfaces due to their
capacity for π-stacking and van der Waals packing.

**Polar (uncharged) residues** — Ser (S), Thr (T), Asn (N), Gln (Q), Cys (C),
Tyr (Y).


### Annotation output files

- `IGHV_FR_CDR_annotations_all_species_all_alleles.csv` — every functional
  IGHV allele, all species
- `IGHV_FR_CDR_annotations_all_species_representative.csv` — one representative
  allele per gene (preferring `*01`), more compact for cross-species work

**Columns:** `species, locus, gene, allele, gene_allele, accession,
functionality, gapped_len, FR1, CDR1, FR2, CDR2, FR3, CDR3_partial,
CDR1_length, CDR2_length, CDR1_hydrophobicity, CDR2_hydrophobicity,
CDR1_CDR2_hydrophobicity, CDR1_aromatic_count, CDR2_aromatic_count,
CDR1_CDR2_aromatic_count, CDR1_polar_count, CDR2_polar_count,
CDR1_CDR2_polar_count, anchor_cys23_ok, anchor_trp41_ok, anchor_cys104_ok,
qc_pass, cdr1_insertion, cdr2_insertion`

### Annotation files

| File | Purpose |
|---|---|
| `annotate_regions.py` | Per-sequence FR/CDR boundary annotator with anchor-residue QC |
| `aa_properties.py` | Hydrophobicity scale + aromatic/polar residue definitions |
| `build_annotation_table.py` | Builds the annotation CSV across species/loci |

---

## Cross-species clustering (PCA / t-SNE)

`feature_matrix.py` + `run_dimred.py` build a ~61-feature vector per V gene
(length, hydrophobicity, aromatic/polar/acidic/basic counts, net charge, and
full 20-amino-acid composition, computed separately for CDR1, CDR2, and
FR1+FR2+FR3 combined), standardize it, and run PCA + t-SNE, colored by species.

`run_dimred_cdr_only.py` is the same pipeline restricted to CDR1/CDR2-derived
features only (drops framework composition) — useful for isolating
functionally-relevant signal from phylogenetic/framework signal. In this
toolkit's own testing, the full feature set gave much tighter species
clustering than CDR-only, because framework composition is more strongly
phylogenetically conserved within a species than CDR composition is.

**A quantitative check worth doing before trusting any 2D projection**: train a
classifier (e.g. `RandomForestClassifier`) directly on the un-reduced,
standardized feature matrix with cross-validation, and compare its accuracy to
random-chance baseline (1/n_species). If accuracy is low, no projection method
— t-SNE, UMAP, or otherwise — will show meaningful clustering, because the
separable signal isn't there in the first place. This toolkit doesn't include
`umap-learn` as a hard dependency, but it's a one-line swap for the `TSNE` call
in `run_dimred.py` if you want to compare.

Example output: [`examples/IGH_Vgene_dimred_by_species.png`](examples/IGH_Vgene_dimred_by_species.png)

---

## Comparing a query sequence against the natural repertoire

`compare_query_to_species.py` takes an arbitrary antibody sequence — e.g. a
computationally designed binder (RFdiffusion + ProteinMPNN output), or any VH/VHH
sequence — and:

1. Extracts its CDR1/CDR2 (and, if a full VH/VHH domain, its full CDR3) using
   the same anchor-residue logic as `annotate_regions.py`. This works directly
   on plain (ungapped) expressed sequences, not just IMGT-gapped germline
   references.
2. Computes its biochemical feature vector in the same schema as the natural
   dataset, and projects it onto the existing PCA space (`sklearn.PCA.transform`
   supports true out-of-sample projection) and a merged t-SNE re-fit (t-SNE has
   no native out-of-sample extension, so the query is included in the fit).
3. Ranks every species by biochemical distance and sequence identity (NW
   alignment) to the query.

`example_designed_loops.py` bundles one real, published test case: `VHH_flu_01`,
an RFdiffusion+ProteinMPNN-designed nanobody targeting influenza hemagglutinin,
cryo-EM validated (PDB [9NH7](https://www.rcsb.org/structure/9NH7)), from
Bennett et al. 2025, *Nature* 649:183-193 (doi:10.1038/s41586-025-09721-5).

```bash
python3 -c "
from example_designed_loops import EXAMPLE_LOOPS
from compare_query_to_species import extract_query_cdrs, compare_query
vhh = EXAMPLE_LOOPS['VHH_flu_01']
cdrs = extract_query_cdrs(vhh['sequence'])
compare_query('VHH_flu_01', cdrs['CDR1'], cdrs['CDR2'])
"
```

**Important scope note:** only CDR1/CDR2 are compared this way. 

Example output: [`examples/VHH_flu_01_CDR3_PCA_projection.png`](examples/VHH_flu_01_CDR3_PCA_projection.png)

---

## Stochastic HCDR3 (V(D)J recombination) simulation

`vdj_recombination.py` simulates germline V(D)J recombination per species:
picks a random functional V/D/J gene, applies random exonuclease trimming at
each junction, inserts random N-nucleotides (TdT-like) between segments, then
translates the full concatenated nucleotide sequence in one continuous
reading frame (D-gene reading frame emerges from trim/insertion lengths rather
than being artificially fixed to "frame 1/2/3"). Sequences with a premature
stop codon are discarded as non-productive.

**Requires the nucleotide bulk file** (`fasta-nt-WithGaps`), not just the AA
one — needed to translate D genes in arbitrary reading frames.

`run_vdj_simulation.py` runs this across all species with usable V/D/J data,
computes the same biochemical feature schema as the CDR1/CDR2 analysis (applied
to the simulated CDR3s), and produces PCA/t-SNE plots + a per-species length
distribution violin plot.

**Explicit caveats (this is an illustrative model, not a calibrated repertoire
simulator):**
- Trim-length and N-insertion-length distributions are uniform-random over a
  fixed range, not fit to any empirical distribution.
- Gene usage is uniform-random; real repertoires show strong, unequal V/D/J
  usage biases.
- No P-nucleotides (palindromic additions at trimmed ends).
- No selection beyond "no premature stop codon" — real repertoires undergo
  further structural/self-tolerance/antigen selection not modeled here.

In this toolkit's own testing, species separated far less cleanly in simulated
CDR3 feature space than in germline CDR1/CDR2 space (~28% vs. ~79% cross-validated
classification accuracy) — consistent with CDR3 diversity being dominated by
stochastic combinatorial/junctional processes rather than phylogenetically
conserved sequence, unlike CDR1/CDR2. One clear biological signal did emerge
correctly despite the simulation's simplicity: *Bos taurus* (cattle) produced a
long right-tailed length distribution (up to 59 aa), consistent with cattle's
well-documented "ultralong" CDR-H3 biology (Wang et al. 2013, *Cell*), driven by
a genuinely long germline D gene already present in the IMGT data.

`compare_query_cdr3_to_species.py` extends the query-comparison tool above to
project a query CDR3 onto this simulated space.

Example output: [`examples/simulated_HCDR3_length_distribution_by_species.png`](examples/simulated_HCDR3_length_distribution_by_species.png)

---

## Citation

If this toolkit is useful in published work, please cite:
- IMGT/GENE-DB as the data source (see `data/README.md`)
- This repository
- Any literature example sequences used (e.g. the VHH_flu_01 citation above)

## Contributing

Issues and pull requests are welcome. 

## License

Code in this repository is MIT-licensed (see `LICENSE`). IMGT data (which you
download separately -- see `data/README.md`) is licensed by IMGT under CC BY
4.0 as of 1 July 2026, permitting commercial use and redistribution with
attribution.
