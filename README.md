# Data

This repository does **not bundle** IMGT/GENE-DB data files directly to keep the repo small and
to make sure you always start from the current IMGT release rather than a
potentially stale copy.

## License

As of **1 July 2026**, IMGT data (including GENE-DB) is provided under the
**Creative Commons Attribution 4.0 license (CC BY 4.0)** to both public and
private/commercial users on IMGT's own terms of use page:
https://www.imgt.org/about/termsofuse.php


**In practice, CC BY 4.0 means:**
- Commercial use is permitted.
- Derivative works (like the processed CSVs and comparisons this toolkit
  produces) are permitted.
- Redistribution is permitted.
- **Attribution is the only requirement.** Cite IMGT when you publish or
  redistribute results derived from their data (see https://www.imgt.org/about/CitingIMGT.php
  for their preferred citation format).

This toolkit's own code is separately MIT-licensed (see `../LICENSE`) and both now permit
commercial use, but attribution to IMGT specifically remains required by their
license regardless of what license covers the code that processes their data.

## Where to get it

Bulk downloads (all species combined) are available at:

https://www.imgt.org/download/GENE-DB/

Two files are used by this toolkit:

| File pattern | Used for |
|---|---|
| `IMGTGENEDB-ReferenceSequences.fasta-AA-WithGaps-F+ORF+inframeP` | Everything except the VDJ recombination simulator: V/D/J/C gene comparison, FR/CDR annotation, dimensionality reduction |
| `IMGTGENEDB-ReferenceSequences.fasta-nt-WithGaps-F+ORF+inframeP` | The VDJ recombination simulator (`vdj_recombination.py`, `run_vdj_simulation.py`); needs nucleotide sequence to translate D genes in multiple reading frames |

Download whichever you need, then either:
- pass the path explicitly: `python3 workflow.py --query "..." --bulk_path /path/to/file`, or
- edit `DEFAULT_BULK` in `workflow.py` (and the equivalent default arguments in
  `run_vdj_simulation.py`, `build_annotation_table.py`, etc.) to point at your
  local copy.

A useful third-party tool for splitting/filtering these bulk files by species
before downloading is [IMGTgeneDL](https://github.com/JamieHeather/IMGTgeneDL).

## Format expected

Standard IMGT/GENE-DB pipe-delimited FASTA headers, e.g.:

```
>M99641|IGHV1-18*01|Homo sapiens|F|V-REGION|188..483|296 nt|1| | | |98 AA|98+8=106| | |
QVQLVQSGA.EVKKPGASVKVSCKASGYTF....TSYGISWVRQAPGQGLEWMGWISAY.
```

`imgt_parser.py` parses this format directly; if IMGT changes their header
schema in the future, that's the one file to update.

