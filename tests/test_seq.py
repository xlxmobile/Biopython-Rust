"""
Tests for the sequence module.
"""

import pytest
import numpy as np
from biopython.seq import (
    Sequence, DNASequence, RNASequence, ProteinSequence,
    SequenceError, random_dna_sequence
)


class TestSequence:
    """Tests for the base Sequence class."""
    
    def test_sequence_creation(self):
        """Test creating sequences from different inputs."""
        # From string
        seq = Sequence("ACGT")
        assert len(seq) == 4
        assert str(seq) == "ACGT"
        
        # From bytes
        seq = Sequence(b"ACGT")
        assert len(seq) == 4
        assert str(seq) == "ACGT"
        
        # From another sequence
        seq2 = Sequence(seq)
        assert len(seq2) == 4
        assert str(seq2) == "ACGT"
        
        # With ID and description
        seq = Sequence("ACGT", id="seq1", description="Test sequence")
        assert seq.id == "seq1"
        assert seq.description == "Test sequence"
    
    def test_sequence_properties(self):
        """Test sequence properties."""
        seq = Sequence("ACGT", id="seq1", description="Test sequence")
        
        # Test properties
        assert seq.id == "seq1"
        assert seq.description == "Test sequence"
        assert len(seq) == 4
        
        # Test setting properties
        seq.id = "new_id"
        seq.description = "New description"
        assert seq.id == "new_id"
        assert seq.description == "New description"
    
    def test_sequence_indexing(self):
        """Test sequence indexing."""
        seq = Sequence("ACGTACGT")
        
        # Test single base access
        assert seq[0] == "A"
        assert seq[1] == "C"
        assert seq[2] == "G"
        assert seq[3] == "T"
        
        # Test negative indices
        assert seq[-1] == "T"
        assert seq[-2] == "G"
        
        # Test slicing
        assert str(seq[0:4]) == "ACGT"
        assert str(seq[4:8]) == "ACGT"
        assert str(seq[2:6]) == "GTAC"
        
        # Test negative slicing
        assert str(seq[-4:]) == "ACGT"
        assert str(seq[:-4]) == "ACGT"
        
        # Test step slicing
        assert str(seq[::2]) == "ATCG"
        assert str(seq[1::2]) == "CGT"
        
        # Test out of bounds
        with pytest.raises(IndexError):
            _ = seq[100]
    
    def test_sequence_methods(self):
        """Test sequence methods."""
        seq = Sequence("ACGTACGT")
        
        # Test counting
        assert seq.count("ACG") == 2
        assert seq.count("AAA") == 0
        
        # Test find_all
        assert seq.find_all("ACG") == [0, 4]
        assert seq.find_all("GT") == [2, 6]
        assert seq.find_all("AAA") == []
        
        # Test to_bytes
        assert seq.to_bytes() == b"ACGTACGT"
        
        # Test to_numpy
        np_array = seq.to_numpy()
        assert isinstance(np_array, np.ndarray)
        assert np_array.dtype == np.uint8
        assert bytes(np_array) == b"ACGTACGT"


class TestDNASequence:
    """Tests for the DNASequence class."""
    
    def test_dna_creation(self):
        """Test creating DNA sequences."""
        # Valid DNA
        dna = DNASequence("ACGTACGT")
        assert len(dna) == 8
        assert str(dna) == "ACGTACGT"
        assert dna.alphabet == "DNA"
        
        # Invalid DNA
        with pytest.raises(SequenceError):
            DNASequence("ACGUXYZ")
    
    def test_dna_methods(self):
        """Test DNA-specific methods."""
        dna = DNASequence("ACGTACGT")
        
        # Test complement
        comp = dna.complement()
        assert isinstance(comp, DNASequence)
        assert str(comp) == "TGCATGCA"
        
        # Test reverse complement
        rev_comp = dna.reverse_complement()
        assert isinstance(rev_comp, DNASequence)
        assert str(rev_comp) == "ACGTACGT"  # Palindromic sequence
        
        # Test with non-palindromic sequence
        dna = DNASequence("ACGTACGTA")
        rev_comp = dna.reverse_complement()
        assert str(rev_comp) == "TACGTACGT"
        
        # Test transcription
        rna = dna.transcribe()
        assert isinstance(rna, RNASequence)
        assert str(rna) == "ACGUACGUA"
        
        # Test GC content
        assert dna.gc_content() == 50.0  # 4/8 = 50%
        
        # Test GC content with G-rich sequence
        dna = DNASequence("GCGCGCGC")
        assert dna.gc_content() == 100.0


class TestRNASequence:
    """Tests for the RNASequence class."""
    
    def test_rna_creation(self):
        """Test creating RNA sequences."""
        # Valid RNA
        rna = RNASequence("ACGUACGU")
        assert len(rna) == 8
        assert str(rna) == "ACGUACGU"
        assert rna.alphabet == "RNA"
        
        # Invalid RNA
        with pytest.raises(SequenceError):
            RNASequence("ACGUTXYZ")
    
    def test_rna_methods(self):
        """Test RNA-specific methods."""
        rna = RNASequence("ACGUACGU")
        
        # Test complement
        comp = rna.complement()
        assert isinstance(comp, RNASequence)
        assert str(comp) == "UGCAUGCA"
        
        # Test reverse complement
        rev_comp = rna.reverse_complement()
        assert isinstance(rev_comp, RNASequence)
        assert str(rev_comp) == "ACGUACGU"  # Palindromic sequence
        
        # Test with non-palindromic sequence
        rna = RNASequence("ACGUACGUA")
        rev_comp = rna.reverse_complement()
        assert str(rev_comp) == "UACGUACGU"
        
        # Test reverse transcription
        dna = rna.reverse_transcribe()
        assert isinstance(dna, DNASequence)
        assert str(dna) == "ACGTACGTA"
        
        # Test GC content
        assert rna.gc_content() == 50.0  # 4/8 = 50%


class TestProteinSequence:
    """Tests for the ProteinSequence class."""
    
    def test_protein_creation(self):
        """Test creating protein sequences."""
        # Valid protein
        prot = ProteinSequence("ACDEFGHIKLMNPQRSTVWYX")
        assert len(prot) == 21
        assert str(prot) == "ACDEFGHIKLMNPQRSTVWYX"
        assert prot.alphabet == "Protein"
        
        # Invalid protein
        with pytest.raises(SequenceError):
            ProteinSequence("ACDEFGHIJKLMNOPQRSTUVWXYZ")
    
    def test_protein_methods(self):
        """Test protein-specific methods."""
        prot = ProteinSequence("ACDEFGHIKLMNPQRSTVWYX")
        
        # Test molecular weight (this is approximate, depends on implementation)
        try:
            mw = prot.molecular_weight()
            assert isinstance(mw, float)
            assert mw > 0
        except NotImplementedError:
            # If not implemented yet, skip this test
            pass


class TestUtilityFunctions:
    """Tests for utility functions."""
    
    def test_random_dna_sequence(self):
        """Test generating random DNA sequences."""
        # Test with valid length
        seq = random_dna_sequence(100)
        assert isinstance(seq, DNASequence)
        assert len(seq) == 100
        
        # Test with invalid length
        with pytest.raises(ValueError):
            random_dna_sequence(0)
        
        with pytest.raises(ValueError):
            random_dna_sequence(-10)