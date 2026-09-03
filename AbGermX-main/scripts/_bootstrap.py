"""Puts ``src/`` on sys.path so the scripts in this folder can import the
library modules directly (``from imgt_parser import ...``).

Every script here does ``import _bootstrap  # noqa: F401`` as its first import.
Python already puts a script's own directory on sys.path, so importing this
module works from any working directory.
"""
import sys
from pathlib import Path

SRC = Path(__file__).resolve().parent.parent / 'src'
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))
