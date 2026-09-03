"""Default input/output locations, and the shared CLI flags for overriding them.

Nothing here is baked into the analysis code: every entry point takes explicit
``bulk_path`` / ``out_dir`` arguments and merely defaults to these values. Set
``ABGERMX_DATA_DIR`` and ``ABGERMX_OUT_DIR`` to point at your own copies without
editing any source.
"""
import os
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent

DATA_DIR = Path(os.environ.get('ABGERMX_DATA_DIR', REPO_ROOT / 'data'))
OUT_DIR = Path(os.environ.get('ABGERMX_OUT_DIR', REPO_ROOT / 'outputs'))

# IMGT/GENE-DB bulk export filenames, exactly as they download from
# https://www.imgt.org/download/GENE-DB/ (see data/README.md).
AA_BULK_NAME = 'IMGTGENEDB-ReferenceSequences.fasta-AA-WithGaps-F+ORF+inframeP'
NT_BULK_NAME = 'IMGTGENEDB-ReferenceSequences.fasta-nt-WithGaps-F+ORF+inframeP'

DEFAULT_AA_BULK = str(DATA_DIR / AA_BULK_NAME)
DEFAULT_NT_BULK = str(DATA_DIR / NT_BULK_NAME)


def ensure_out_dir(out_dir=None):
    """Resolve an output directory (defaulting to OUT_DIR) and create it."""
    path = Path(out_dir) if out_dir is not None else OUT_DIR
    path.mkdir(parents=True, exist_ok=True)
    return str(path)
