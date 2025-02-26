"""
I/O module for reading and writing sequence files.

This module provides functions for reading and writing various sequence file formats.
"""

from typing import List, Optional, Dict, Union, Any, Tuple, Iterator, Iterable
from pathlib import Path
import biopython_rust._rust_bindings as _rust
from biopython_rust.seq import Sequence, DNASequence, RNASequence, ProteinSequence, SequenceError


class FastaRecord:
    """A record from a FASTA file."""
    
    def __init__(self, id: str, sequence: Sequence, description: Optional[str] = None) -> None:
        """
        Initialize a FASTA record.
        
        Args:
            id: The sequence identifier (without the '>')
            sequence: The sequence
            description: Optional description
        """
        self.id = id
        self.sequence = sequence
        self.description = description
    
    def __str__(self) -> str:
        """Get the string representation of the record in FASTA format."""
        header = f">{self.id}"
        if self.description:
            header += f" {self.description}"
        
        # Format sequence with 60 characters per line
        seq_str = str(self.sequence)
        formatted_seq = '\n'.join(seq_str[i:i+60] for i in range(0, len(seq_str), 60))
        
        return f"{header}\n{formatted_seq}"
    
    def __repr__(self) -> str:
        """Get a detailed representation of the record."""
        desc_str = f", description='{self.description}'" if self.description else ""
        return f"FastaRecord(id='{self.id}'{desc_str}, sequence={repr(self.sequence)})"


def read_fasta(path: Union[str, Path]) -> List[FastaRecord]:
    """
    Read sequences from a FASTA file.
    
    Args:
        path: Path to the FASTA file
        
    Returns:
        List of FastaRecord objects
        
    Raises:
        SequenceError: If the file cannot be read or parsed
    """
    if isinstance(path, str):
        path = Path(path)
    
    try:
        fasta_records = _rust.read_fasta(str(path))
    except RuntimeError as e:
        raise SequenceError(f"Failed to read FASTA file: {e}")
    
    records = []
    for rust_record in fasta_records:
        seq_id = rust_record["id"]
        description = rust_record.get("description")
        seq_data = rust_record["sequence"]
        
        # Try to detect the sequence type
        alphabet = _rust.get_sequence_alphabet(seq_data)
        
        if alphabet == "DNA":
            sequence = DNASequence(seq_data, id=seq_id, description=description)
        elif alphabet == "RNA":
            sequence = RNASequence(seq_data, id=seq_id, description=description)
        elif alphabet == "Protein":
            sequence = ProteinSequence(seq_data, id=seq_id, description=description)
        else:
            sequence = Sequence(seq_data, id=seq_id, description=description)
        
        records.append(FastaRecord(seq_id, sequence, description))
    
    return records


def write_fasta(records: List[FastaRecord], path: Union[str, Path]) -> None:
    """
    Write sequences to a FASTA file.
    
    Args:
        records: List of FastaRecord objects
        path: Path to write the FASTA file
        
    Raises:
        SequenceError: If the file cannot be written
    """
    if isinstance(path, str):
        path = Path(path)
    
    # Convert to Rust-compatible format
    rust_records = []
    for record in records:
        rust_records.append({
            "id": record.id,
            "description": record.description,
            "sequence": record.sequence._rust_seq
        })
    
    try:
        _rust.write_fasta(rust_records, str(path))
    except RuntimeError as e:
        raise SequenceError(f"Failed to write FASTA file: {e}")


def detect_format(path: Union[str, Path]) -> str:
    """
    Detect the format of a sequence file.
    
    Args:
        path: Path to the sequence file
        
    Returns:
        The format name ("FASTA", "FASTQ", etc.)
        
    Raises:
        SequenceError: If the format cannot be detected
    """
    if isinstance(path, str):
        path = Path(path)
    
    try:
        return _rust.detect_file_format(str(path))
    except RuntimeError as e:
        raise SequenceError(f"Failed to detect file format: {e}")


# Additional reader functions that could be implemented:
# - read_fastq
# - read_genbank
# - read_embl
# - read_clustal
# - etc.