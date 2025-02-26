"""
Sequence module for handling biological sequences.

This module provides classes and functions for working with DNA, RNA, and protein sequences.
"""

from typing import List, Optional, Dict, Union, Any, Tuple, Iterator, Iterable
import numpy as np
from pathlib import Path
import biopython_rust._rust_bindings as _rust

class SequenceError(Exception):
    """Exception raised for sequence-related errors."""
    pass


class Sequence:
    """Base class for biological sequences."""
    
    def __init__(self, sequence: Union[str, bytes, "Sequence"], id: Optional[str] = None, 
                 description: Optional[str] = None) -> None:
        """
        Initialize a sequence.
        
        Args:
            sequence: The sequence as a string, bytes, or another Sequence object
            id: Optional identifier for the sequence
            description: Optional description for the sequence
        
        Raises:
            SequenceError: If the sequence is invalid
        """
        # Convert input to appropriate format
        if isinstance(sequence, Sequence):
            self._rust_seq = sequence._rust_seq
        elif isinstance(sequence, str):
            try:
                self._rust_seq = _rust.create_sequence(sequence.encode('utf-8'))
            except RuntimeError as e:
                raise SequenceError(str(e))
        elif isinstance(sequence, bytes):
            try:
                self._rust_seq = _rust.create_sequence(sequence)
            except RuntimeError as e:
                raise SequenceError(str(e))
        else:
            raise TypeError(f"Sequence must be a string, bytes, or Sequence object, not {type(sequence)}")
        
        # Set ID and description if provided
        if id is not None:
            self._rust_seq = _rust.set_sequence_id(self._rust_seq, id)
        
        if description is not None:
            self._rust_seq = _rust.set_sequence_description(self._rust_seq, description)
    
    @property
    def id(self) -> Optional[str]:
        """Get the sequence identifier."""
        return _rust.get_sequence_id(self._rust_seq)
    
    @id.setter
    def id(self, value: str) -> None:
        """Set the sequence identifier."""
        self._rust_seq = _rust.set_sequence_id(self._rust_seq, value)
    
    @property
    def description(self) -> Optional[str]:
        """Get the sequence description."""
        return _rust.get_sequence_description(self._rust_seq)
    
    @description.setter
    def description(self, value: str) -> None:
        """Set the sequence description."""
        self._rust_seq = _rust.set_sequence_description(self._rust_seq, value)
    
    def __len__(self) -> int:
        """Get the length of the sequence."""
        return _rust.get_sequence_length(self._rust_seq)
    
    def __str__(self) -> str:
        """Get the string representation of the sequence."""
        return _rust.sequence_to_string(self._rust_seq)
    
    def __repr__(self) -> str:
        """Get a detailed representation of the sequence."""
        id_str = f" id='{self.id}'" if self.id else ""
        desc_str = f" description='{self.description}'" if self.description else ""
        alphabet = self.alphabet
        return f"{self.__class__.__name__}('{str(self)[:20]}...'{id_str}{desc_str} length={len(self)} alphabet={alphabet})"
    
    def __getitem__(self, key: Union[int, slice]) -> Union[str, "Sequence"]:
        """
        Get a subsequence using slice notation.
        
        Args:
            key: An integer index or slice
            
        Returns:
            A single character (for integer index) or a subsequence (for slice)
        """
        if isinstance(key, int):
            # Handle negative indices
            if key < 0:
                key = len(self) + key
            
            # Check bounds
            if key < 0 or key >= len(self):
                raise IndexError(f"Index {key} out of range for sequence of length {len(self)}")
            
            # Get a single character
            return _rust.get_sequence_base(self._rust_seq, key).decode('utf-8')
        
        elif isinstance(key, slice):
            # Convert slice to start/end indices
            start, stop, step = key.indices(len(self))
            
            if step != 1:
                # For non-contiguous slices, we need to handle it in Python
                result = ''.join(self[i] for i in range(start, stop, step))
                return self.__class__(result)
            
            # For contiguous slices, use the Rust implementation
            try:
                subsequence = _rust.get_subsequence(self._rust_seq, start, stop)
                return self.__class__(subsequence)
            except RuntimeError as e:
                raise SequenceError(str(e))
        
        else:
            raise TypeError(f"Invalid key type: {type(key)}")
    
    @property
    def alphabet(self) -> str:
        """Get the alphabet name for this sequence."""
        return _rust.get_sequence_alphabet(self._rust_seq)
    
    def count(self, pattern: Union[str, bytes]) -> int:
        """
        Count the occurrences of a pattern in the sequence.
        
        Args:
            pattern: The pattern to count
            
        Returns:
            The number of occurrences
        """
        if isinstance(pattern, str):
            pattern = pattern.encode('utf-8')
        
        return _rust.count_pattern(self._rust_seq, pattern)
    
    def find_all(self, pattern: Union[str, bytes]) -> List[int]:
        """
        Find all occurrences of a pattern in the sequence.
        
        Args:
            pattern: The pattern to find
            
        Returns:
            List of starting positions (0-based)
        """
        if isinstance(pattern, str):
            pattern = pattern.encode('utf-8')
        
        return _rust.find_all_patterns(self._rust_seq, pattern)
    
    def to_bytes(self) -> bytes:
        """Get the raw bytes of the sequence."""
        return _rust.sequence_to_bytes(self._rust_seq)
    
    def to_numpy(self) -> np.ndarray:
        """Get the sequence as a NumPy array of bytes."""
        return np.frombuffer(self.to_bytes(), dtype=np.uint8)


class DNASequence(Sequence):
    """DNA sequence class."""
    
    def __init__(self, sequence: Union[str, bytes, Sequence], id: Optional[str] = None, 
                 description: Optional[str] = None) -> None:
        """
        Initialize a DNA sequence.
        
        Args:
            sequence: The sequence as a string, bytes, or another Sequence object
            id: Optional identifier for the sequence
            description: Optional description for the sequence
        
        Raises:
            SequenceError: If the sequence is not a valid DNA sequence
        """
        # Convert to bytes if needed
        if isinstance(sequence, str):
            sequence = sequence.encode('utf-8')
        elif isinstance(sequence, Sequence):
            sequence = sequence.to_bytes()
        
        # Create a DNA sequence using the Rust implementation
        try:
            self._rust_seq = _rust.create_dna_sequence(sequence)
        except RuntimeError as e:
            raise SequenceError(str(e))
        
        # Set ID and description if provided
        if id is not None:
            self._rust_seq = _rust.set_sequence_id(self._rust_seq, id)
        
        if description is not None:
            self._rust_seq = _rust.set_sequence_description(self._rust_seq, description)
    
    def complement(self) -> "DNASequence":
        """
        Get the complement of the DNA sequence.
        
        Returns:
            A new DNASequence with the complementary bases
        """
        try:
            complemented = _rust.complement_sequence(self._rust_seq)
            return DNASequence(complemented)
        except RuntimeError as e:
            raise SequenceError(str(e))
    
    def reverse_complement(self) -> "DNASequence":
        """
        Get the reverse complement of the DNA sequence.
        
        Returns:
            A new DNASequence with the reverse complementary bases
        """
        try:
            rev_comp = _rust.reverse_complement_sequence(self._rust_seq)
            return DNASequence(rev_comp)
        except RuntimeError as e:
            raise SequenceError(str(e))
    
    def transcribe(self) -> "RNASequence":
        """
        Transcribe the DNA sequence to RNA.
        
        Returns:
            An RNASequence
        """
        try:
            rna = _rust.transcribe_sequence(self._rust_seq)
            return RNASequence(rna)
        except RuntimeError as e:
            raise SequenceError(str(e))
    
    def gc_content(self) -> float:
        """
        Calculate the GC content of the sequence.
        
        Returns:
            The percentage of G and C bases in the sequence
        """
        try:
            return _rust.gc_content(self._rust_seq)
        except RuntimeError as e:
            raise SequenceError(str(e))


class RNASequence(Sequence):
    """RNA sequence class."""
    
    def __init__(self, sequence: Union[str, bytes, Sequence], id: Optional[str] = None, 
                 description: Optional[str] = None) -> None:
        """
        Initialize an RNA sequence.
        
        Args:
            sequence: The sequence as a string, bytes, or another Sequence object
            id: Optional identifier for the sequence
            description: Optional description for the sequence
        
        Raises:
            SequenceError: If the sequence is not a valid RNA sequence
        """
        # Convert to bytes if needed
        if isinstance(sequence, str):
            sequence = sequence.encode('utf-8')
        elif isinstance(sequence, Sequence):
            sequence = sequence.to_bytes()
        
        # Create an RNA sequence using the Rust implementation
        try:
            self._rust_seq = _rust.create_rna_sequence(sequence)
        except RuntimeError as e:
            raise SequenceError(str(e))
        
        # Set ID and description if provided
        if id is not None:
            self._rust_seq = _rust.set_sequence_id(self._rust_seq, id)
        
        if description is not None:
            self._rust_seq = _rust.set_sequence_description(self._rust_seq, description)
    
    def complement(self) -> "RNASequence":
        """
        Get the complement of the RNA sequence.
        
        Returns:
            A new RNASequence with the complementary bases
        """
        try:
            complemented = _rust.complement_sequence(self._rust_seq)
            return RNASequence(complemented)
        except RuntimeError as e:
            raise SequenceError(str(e))
    
    def reverse_complement(self) -> "RNASequence":
        """
        Get the reverse complement of the RNA sequence.
        
        Returns:
            A new RNASequence with the reverse complementary bases
        """
        try:
            rev_comp = _rust.reverse_complement_sequence(self._rust_seq)
            return RNASequence(rev_comp)
        except RuntimeError as e:
            raise SequenceError(str(e))
    
    def reverse_transcribe(self) -> DNASequence:
        """
        Reverse transcribe the RNA sequence to DNA.
        
        Returns:
            A DNASequence
        """
        try:
            dna = _rust.reverse_transcribe_sequence(self._rust_seq)
            return DNASequence(dna)
        except RuntimeError as e:
            raise SequenceError(str(e))
    
    def gc_content(self) -> float:
        """
        Calculate the GC content of the sequence.
        
        Returns:
            The percentage of G and C bases in the sequence
        """
        try:
            return _rust.gc_content(self._rust_seq)
        except RuntimeError as e:
            raise SequenceError(str(e))


class ProteinSequence(Sequence):
    """Protein sequence class."""
    
    def __init__(self, sequence: Union[str, bytes, Sequence], id: Optional[str] = None, 
                 description: Optional[str] = None) -> None:
        """
        Initialize a protein sequence.
        
        Args:
            sequence: The sequence as a string, bytes, or another Sequence object
            id: Optional identifier for the sequence
            description: Optional description for the sequence
        
        Raises:
            SequenceError: If the sequence is not a valid protein sequence
        """
        # Convert to bytes if needed
        if isinstance(sequence, str):
            sequence = sequence.encode('utf-8')
        elif isinstance(sequence, Sequence):
            sequence = sequence.to_bytes()
        
        # Create a protein sequence using the Rust implementation
        try:
            self._rust_seq = _rust.create_protein_sequence(sequence)
        except RuntimeError as e:
            raise SequenceError(str(e))
        
        # Set ID and description if provided
        if id is not None:
            self._rust_seq = _rust.set_sequence_id(self._rust_seq, id)
        
        if description is not None:
            self._rust_seq = _rust.set_sequence_description(self._rust_seq, description)
    
    def molecular_weight(self) -> float:
        """
        Calculate the molecular weight of the protein.
        
        Returns:
            The molecular weight in Daltons
        """
        try:
            return _rust.protein_molecular_weight(self._rust_seq)
        except RuntimeError as e:
            raise SequenceError(str(e))


# Utility functions
def random_dna_sequence(length: int) -> DNASequence:
    """
    Generate a random DNA sequence of the specified length.
    
    Args:
        length: The length of the sequence to generate
    
    Returns:
        A random DNASequence
    """
    if length <= 0:
        raise ValueError("Length must be positive")
    
    try:
        random_seq = _rust.random_dna_sequence(length)
        return DNASequence(random_seq)
    except RuntimeError as e:
        raise SequenceError(str(e))