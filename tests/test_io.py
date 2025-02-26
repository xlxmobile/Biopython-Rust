"""
Tests for the I/O module.
"""

import pytest
import os
import tempfile
from pathlib import Path

from biopython.seq import DNASequence, RNASequence, ProteinSequence
from biopython.io import FastaRecord, read_fasta, write_fasta, detect_format, SequenceError


class TestFastaRecord:
    """Tests for the FastaRecord class."""
    
    def test_fasta_record_creation(self):
        """Test creating FASTA records."""
        # Create a DNA sequence
        dna = DNASequence("ACGTACGT")
        
        # Create a FASTA record
        record = FastaRecord("seq1", dna, "Test sequence")
        
        # Test properties
        assert record.id == "seq1"
        assert record.description == "Test sequence"
        assert str(record.sequence) == "ACGTACGT"
    
    def test_fasta_record_string_representation(self):
        """Test string representation of FASTA records."""
        # Create a DNA sequence
        dna = DNASequence("ACGTACGT")
        
        # Create a FASTA record
        record = FastaRecord("seq1", dna, "Test sequence")
        
        # Test string representation
        fasta_str = str(record)
        assert fasta_str.startswith(">seq1 Test sequence\n")
        assert "ACGTACGT" in fasta_str
        
        # Create a record without description
        record = FastaRecord("seq1", dna)
        
        # Test string representation
        fasta_str = str(record)
        assert fasta_str.startswith(">seq1\n")
        assert "ACGTACGT" in fasta_str


class TestFastaIO:
    """Tests for FASTA I/O functions."""
    
    def test_read_write_fasta(self):
        """Test reading and writing FASTA files."""
        # Create temporary file
        with tempfile.NamedTemporaryFile(suffix=".fasta", delete=False) as temp:
            temp_path = temp.name
            
            # Write a simple FASTA file
            temp.write(b">seq1 First sequence\nACGTACGT\n")
            temp.write(b">seq2 Second sequence\nGTACGTAC\n")
            temp.close()
            
            try:
                # Read the FASTA file
                records = read_fasta(temp_path)
                
                # Test the records
                assert len(records) == 2
                assert records[0].id == "seq1"
                assert records[0].description == "First sequence"
                assert str(records[0].sequence) == "ACGTACGT"
                
                assert records[1].id == "seq2"
                assert records[1].description == "Second sequence"
                assert str(records[1].sequence) == "GTACGTAC"
                
                # Create a new temp file for writing
                write_path = temp_path + ".out"
                
                # Write the records to a new file
                write_fasta(records, write_path)
                
                # Read the written file
                new_records = read_fasta(write_path)
                
                # Test the new records
                assert len(new_records) == 2
                assert new_records[0].id == "seq1"
                assert new_records[0].description == "First sequence"
                assert str(new_records[0].sequence) == "ACGTACGT"
                
                assert new_records[1].id == "seq2"
                assert new_records[1].description == "Second sequence"
                assert str(new_records[1].sequence) == "GTACGTAC"
                
            finally:
                # Clean up
                os.unlink(temp_path)
                if os.path.exists(write_path):
                    os.unlink(write_path)
    
    def test_read_fasta_errors(self):
        """Test error handling when reading FASTA files."""
        # Test with non-existent file
        with pytest.raises(SequenceError):
            read_fasta("non_existent_file.fasta")
        
        # Test with invalid FASTA
        with tempfile.NamedTemporaryFile(suffix=".fasta", delete=False) as temp:
            temp_path = temp.name
            
            # Write an invalid FASTA file (missing header)
            temp.write(b"ACGTACGT\n")
            temp.close()
            
            try:
                with pytest.raises(SequenceError):
                    read_fasta(temp_path)
            finally:
                os.unlink(temp_path)
    
    def test_detect_format(self):
        """Test detecting file formats."""
        # Create temporary FASTA file
        with tempfile.NamedTemporaryFile(suffix=".fasta", delete=False) as temp:
            temp_path = temp.name
            
            # Write a simple FASTA file
            temp.write(b">seq1\nACGTACGT\n")
            temp.close()
            
            try:
                # Detect format
                format_name = detect_format(temp_path)
                assert format_name == "FASTA"
            finally:
                os.unlink(temp_path)
        
        # Test with non-existent file
        with pytest.raises(SequenceError):
            detect_format("non_existent_file.txt")