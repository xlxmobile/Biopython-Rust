//! FASTA format I/O
//!
//! This module provides functions for reading and writing FASTA files.

use std::path::Path;
use std::io::{self, BufRead, Write};
use std::fs::File;

use crate::engines::storage::formats::{FastaParser, FastaWriter, SequenceRecord};
use crate::modules::seq::{Sequence, SequenceError};

/// A FASTA record
#[derive(Debug, Clone)]
pub struct FastaRecord {
    /// Sequence identifier
    pub id: String,
    /// Optional sequence description
    pub description: Option<String>,
    /// The sequence
    pub sequence: Sequence,
}

impl FastaRecord {
    /// Create a new FASTA record
    pub fn new(id: &str, description: Option<&str>, sequence: Sequence) -> Self {
        Self {
            id: id.to_string(),
            description: description.map(|s| s.to_string()),
            sequence,
        }
    }
    
    /// Convert to a string in FASTA format
    pub fn to_string(&self) -> String {
        let mut output = String::new();
        
        // Write header
        match &self.description {
            Some(desc) => output.push_str(&format!(">{} {}\n", self.id, desc)),
            None => output.push_str(&format!(">{}\n", self.id)),
        };
        
        // Write sequence with line wrapping (60 characters per line)
        let seq_str = self.sequence.as_string();
        for i in (0..seq_str.len()).step_by(60) {
            let end = (i + 60).min(seq_str.len());
            output.push_str(&seq_str[i..end]);
            output.push('\n');
        }
        
        output
    }
}

/// Read sequences from a FASTA file
pub fn read_fasta<P: AsRef<Path>>(path: P) -> Result<Vec<FastaRecord>, SequenceError> {
    let parser = FastaParser::new();
    let engine_records = parser.parse_file(path)
        .map_err(|e| SequenceError::EngineError(e))?;
    
    let mut records = Vec::with_capacity(engine_records.len());
    
    for record in engine_records {
        // Convert engine record to Sequence
        let sequence = Sequence::new(&record.sequence_as_vec())?
            .with_id(&record.id);
        
        let sequence = if let Some(desc) = &record.description {
            sequence.with_description(desc)
        } else {
            sequence
        };
        
        records.push(FastaRecord {
            id: record.id.clone(),
            description: record.description.clone(),
            sequence,
        });
    }
    
    Ok(records)
}

/// Write sequences to a FASTA file
pub fn write_fasta<P: AsRef<Path>>(records: &[FastaRecord], path: P) -> Result<(), SequenceError> {
    let writer = FastaWriter::new();
    
    // Convert FastaRecord to SequenceRecord
    let engine_records: Vec<SequenceRecord> = records.iter().map(|record| {
        SequenceRecord::new(
            record.id.clone(),
            record.description.clone(),
            record.sequence.as_bytes().to_vec(),
        )
    }).collect();
    
    // Write using the engine writer
    writer.write_file(&engine_records, path)
        .map_err(|e| SequenceError::EngineError(e))?;
    
    Ok(())
}

/// Read sequences from a FASTA string
pub fn read_fasta_string(content: &str) -> Result<Vec<FastaRecord>, SequenceError> {
    let parser = FastaParser::new();
    let engine_records = parser.parse_string(content)
        .map_err(|e| SequenceError::EngineError(e))?;
    
    let mut records = Vec::with_capacity(engine_records.len());
    
    for record in engine_records {
        // Convert engine record to Sequence
        let sequence = Sequence::new(&record.sequence_as_vec())?
            .with_id(&record.id);
        
        let sequence = if let Some(desc) = &record.description {
            sequence.with_description(desc)
        } else {
            sequence
        };
        
        records.push(FastaRecord {
            id: record.id.clone(),
            description: record.description.clone(),
            sequence,
        });
    }
    
    Ok(records)
}

/// Write sequences to a FASTA string
pub fn write_fasta_string(records: &[FastaRecord]) -> Result<String, SequenceError> {
    let writer = FastaWriter::new();
    
    // Convert FastaRecord to SequenceRecord
    let engine_records: Vec<SequenceRecord> = records.iter().map(|record| {
        SequenceRecord::new(
            record.id.clone(),
            record.description.clone(),
            record.sequence.as_bytes().to_vec(),
        )
    }).collect();
    
    // Write using the engine writer
    let content = writer.write_string(&engine_records)
        .map_err(|e| SequenceError::EngineError(e))?;
    
    Ok(content)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;
    use tempfile::tempdir;
    
    #[test]
    fn test_fasta_record() {
        // Create a DNA sequence
        let seq = Sequence::new_dna(b"ACGTACGT").unwrap();
        
        // Create a FASTA record
        let record = FastaRecord::new("seq1", Some("Test sequence"), seq);
        
        // Test properties
        assert_eq!(record.id, "seq1");
        assert_eq!(record.description, Some("Test sequence".to_string()));
        assert_eq!(record.sequence.as_bytes().as_ref(), b"ACGTACGT");
        
        // Test to_string
        let fasta_str = record.to_string();
        assert!(fasta_str.starts_with(">seq1 Test sequence\n"));
        assert!(fasta_str.contains("ACGTACGT"));
    }
    
    #[test]
    fn test_read_write_fasta() -> std::io::Result<()> {
        // Create a temporary directory
        let dir = tempdir()?;
        let file_path = dir.path().join("test.fasta");
        
        // Create a FASTA file
        let fasta_content = ">seq1 First sequence\nACGTACGT\n>seq2 Second sequence\nGTACGTAC\n";
        {
            let mut file = File::create(&file_path)?;
            file.write_all(fasta_content.as_bytes())?;
        }
        
        // Read the FASTA file
        let records = read_fasta(&file_path).unwrap();
        
        // Check the records
        assert_eq!(records.len(), 2);
        assert_eq!(records[0].id, "seq1");
        assert_eq!(records[0].description, Some("First sequence".to_string()));
        assert_eq!(records[0].sequence.as_bytes().as_ref(), b"ACGTACGT");
        
        assert_eq!(records[1].id, "seq2");
        assert_eq!(records[1].description, Some("Second sequence".to_string()));
        assert_eq!(records[1].sequence.as_bytes().as_ref(), b"GTACGTAC");
        
        // Write the records to a new file
        let output_path = dir.path().join("output.fasta");
        write_fasta(&records, &output_path).unwrap();
        
        // Read the new file and check if the records are the same
        let records2 = read_fasta(&output_path).unwrap();
        
        assert_eq!(records2.len(), 2);
        assert_eq!(records2[0].id, "seq1");
        assert_eq!(records2[0].description, Some("First sequence".to_string()));
        assert_eq!(records2[0].sequence.as_bytes().as_ref(), b"ACGTACGT");
        
        assert_eq!(records2[1].id, "seq2");
        assert_eq!(records2[1].description, Some("Second sequence".to_string()));
        assert_eq!(records2[1].sequence.as_bytes().as_ref(), b"GTACGTAC");
        
        Ok(())
    }
    
    #[test]
    fn test_read_write_fasta_string() {
        // Create a FASTA string
        let fasta_content = ">seq1 First sequence\nACGTACGT\n>seq2 Second sequence\nGTACGTAC\n";
        
        // Read the FASTA string
        let records = read_fasta_string(fasta_content).unwrap();
        
        // Check the records
        assert_eq!(records.len(), 2);
        assert_eq!(records[0].id, "seq1");
        assert_eq!(records[0].description, Some("First sequence".to_string()));
        assert_eq!(records[0].sequence.as_bytes().as_ref(), b"ACGTACGT");
        
        // Write the records to a string
        let output = write_fasta_string(&records).unwrap();
        
        // Check the output string
        assert!(output.contains(">seq1 First sequence"));
        assert!(output.contains("ACGTACGT"));
        assert!(output.contains(">seq2 Second sequence"));
        assert!(output.contains("GTACGTAC"));
    }
}