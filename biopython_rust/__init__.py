"""
BioPython: High-performance bioinformatics sequence operations library.

This package provides Python bindings for the Rust-based bioseq library.
"""

__version__ = "0.1.0"

# Import submodules
from biopython_rust import seq
from biopython_rust import io

# Core functionality re-exported at the top level
from biopython_rust.seq import Sequence, DNASequence, RNASequence, ProteinSequence
from biopython_rust.io import read_fasta, write_fasta, FastaRecord

# Initialize the library
import biopython_rust._rust_bindings as _rust
_rust.init()

__all__ = [
    "seq",
    "io",
    "Sequence",
    "DNASequence",
    "RNASequence",
    "ProteinSequence",
    "read_fasta",
    "write_fasta",
    "FastaRecord",
]