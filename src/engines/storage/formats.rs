//! Parsers and writers for common bioinformatics file formats
//!
//! This module provides high-performance parsers and writers for various
//! bioinformatics file formats, including FASTA, FASTQ, etc.

use std::path::Path;
use std::io::{self, BufRead, Write};
use std::collections::HashMap;
use crate::engines::EngineResult;
use crate::engines::EngineError;
use crate::engines::core::io::{FastReader, FastWriter};
use crate::engines::storage::{StorableSequence, InMemoryStorage, StorageFactory, StorageMode};

/// Trait for sequence record parsers
pub trait SequenceParser: Send + Sync {
    /// Parse a file and create sequence records
    fn parse_file<P: AsRef<Path>>(&self, path: P) -> EngineResult<Vec<SequenceRecord>>;
    
    /// Parse a string and create sequence records
    fn parse_string(&self, content: &str) -> EngineResult<Vec<SequenceRecord>>;
    
    /// Get the format name
    fn format_name(&self) -> &str;
}

/// Trait for sequence record writers
pub trait SequenceWriter: Send + Sync {
    /// Write sequence records to a file
    fn write_file<P: AsRef<Path>>(&self, records: &[SequenceRecord], path: P) -> EngineResult<()>;
    
    /// Write sequence records to a string
    fn write_string(&self, records: &[SequenceRecord]) -> EngineResult<String>;
    
    /// Get the format name
    fn format_name(&self) -> &str;
}

/// A sequence record with ID, description, and sequence data
#[derive(Debug, Clone)]
pub struct SequenceRecord {
    /// Sequence identifier
    pub id: String,
    /// Optional sequence description
    pub description: Option<String>,
    /// The sequence data storage
    pub sequence: Box<dyn StorableSequence>,
    /// Optional quality scores (for formats like FASTQ)
    pub quality: Option<Box<dyn StorableSequence>>,
    /// Optional metadata as key-value pairs
    pub metadata: HashMap<String, String>,
}

impl SequenceRecord {
    /// Create a new in-memory sequence record
    pub fn new(
        id: String,
        description: Option<String>,
        sequence: Vec<u8>,
    ) -> Self {
        Self {
            id,
            description,
            sequence: Box::new(InMemoryStorage::new(sequence)),
            quality: None,
            metadata: HashMap::new(),
        }
    }
    
    /// Create a new in-memory sequence record with quality scores
    pub fn with_quality(
        id: String,
        description: Option<String>,
        sequence: Vec<u8>,
        quality: Vec<u8>,
    ) -> Self {
        Self {
            id,
            description,
            sequence: Box::new(InMemoryStorage::new(sequence)),
            quality: Some(Box::new(InMemoryStorage::new(quality))),
            metadata: HashMap::new(),
        }
    }
    
    /// Get the length of the sequence
    pub fn len(&self) -> usize {
        self.sequence.len()
    }
    
    /// Check if the sequence is empty
    pub fn is_empty(&self) -> bool {
        self.sequence.is_empty()
    }
    
    /// Get a subsequence
    pub fn subsequence(&self, start: usize, end: usize) -> Vec<u8> {
        self.sequence.subsequence(start, end)
    }
    
    /// Get the full sequence as a vector
    pub fn sequence_as_vec(&self) -> Vec<u8> {
        match self.sequence.as_slice() {
            Some(slice) => slice.to_vec(),
            None => self.sequence.subsequence(0, self.sequence.len()),
        }
    }
    
    /// Get the full quality scores as a vector (if available)
    pub fn quality_as_vec(&self) -> Option<Vec<u8>> {
        self.quality.as_ref().map(|q| {
            match q.as_slice() {
                Some(slice) => slice.to_vec(),
                None => q.subsequence(0, q.len()),
            }
        })
    }
    
    /// Add metadata to the record
    pub fn add_metadata(&mut self, key: &str, value: &str) {
        self.metadata.insert(key.to_string(), value.to_string());
    }
    
    /// Get metadata value
    pub fn get_metadata(&self, key: &str) -> Option<&String> {
        self.metadata.get(key)
    }
}

/// FASTA format parser
#[derive(Debug, Clone)]
pub struct FastaParser {
    /// Storage mode to use for sequences
    storage_mode: StorageMode,
    /// Buffer size for reading
    buffer_size: usize,
}

impl FastaParser {
    /// Create a new FASTA parser with the default storage mode
    pub fn new() -> Self {
        Self {
            storage_mode: StorageMode::default(),
            buffer_size: 1024 * 1024, // 1MB
        }
    }
    
    /// Create a new FASTA parser with the specified storage mode
    pub fn with_storage_mode(storage_mode: StorageMode) -> Self {
        Self {
            storage_mode,
            buffer_size: 1024 * 1024, // 1MB
        }
    }
    
    /// Set the buffer size
    pub fn with_buffer_size(mut self, buffer_size: usize) -> Self {
        self.buffer_size = buffer_size;
        self
    }
}

impl Default for FastaParser {
    fn default() -> Self {
        Self::new()
    }
}

impl SequenceParser for FastaParser {
    fn parse_file<P: AsRef<Path>>(&self, path: P) -> EngineResult<Vec<SequenceRecord>> {
        let mut reader = FastReader::new(path.as_ref(), Some(self.buffer_size))?;
        
        let mut records = Vec::new();
        let mut current_id = String::new();
        let mut current_desc = None;
        let mut current_seq = Vec::new();
        
        for line_result in reader.read_lines() {
            let line = line_result?;
            
            // Skip empty lines
            if line.is_empty() {
                continue;
            }
            
            // Header line
            if line.starts_with('>') {
                // Save the previous record if any
                if !current_id.is_empty() && !current_seq.is_empty() {
                    // Create storage according to the chosen mode
                    let sequence = StorageFactory::create_storage(
                        Some(current_seq.clone()),
                        Some(path.as_ref()),
                        Some(current_seq.len()),
                        Some(self.storage_mode),
                    )?;
                    
                    records.push(SequenceRecord {
                        id: current_id.clone(),
                        description: current_desc.clone(),
                        sequence,
                        quality: None,
                        metadata: HashMap::new(),
                    });
                    
                    current_seq.clear();
                }
                
                // Parse header
                let header = &line[1..];
                let parts: Vec<&str> = header.splitn(2, ' ').collect();
                
                current_id = parts[0].to_string();
                current_desc = parts.get(1).map(|s| s.to_string());
            } else {
                // Sequence line (add to current sequence)
                current_seq.extend(line.trim().as_bytes());
            }
        }
        
        // Add the last record if any
        if !current_id.is_empty() && !current_seq.is_empty() {
            // Create storage according to the chosen mode
            let sequence = StorageFactory::create_storage(
                Some(current_seq.clone()),
                Some(path.as_ref()),
                Some(current_seq.len()),
                Some(self.storage_mode),
            )?;
            
            records.push(SequenceRecord {
                id: current_id,
                description: current_desc,
                sequence,
                quality: None,
                metadata: HashMap::new(),
            });
        }
        
        Ok(records)
    }
    
    fn parse_string(&self, content: &str) -> EngineResult<Vec<SequenceRecord>> {
        let mut records = Vec::new();
        let mut current_id = String::new();
        let mut current_desc = None;
        let mut current_seq = Vec::new();
        
        for line in content.lines() {
            let line = line.trim();
            
            // Skip empty lines
            if line.is_empty() {
                continue;
            }
            
            // Header line
            if line.starts_with('>') {
                // Save the previous record if any
                if !current_id.is_empty() && !current_seq.is_empty() {
                    records.push(SequenceRecord::new(
                        current_id.clone(),
                        current_desc.clone(),
                        current_seq.clone(),
                    ));
                    
                    current_seq.clear();
                }
                
                // Parse header
                let header = &line[1..];
                let parts: Vec<&str> = header.splitn(2, ' ').collect();
                
                current_id = parts[0].to_string();
                current_desc = parts.get(1).map(|s| s.to_string());
            } else {
                // Sequence line (add to current sequence)
                current_seq.extend(line.as_bytes());
            }
        }
        
        // Add the last record if any
        if !current_id.is_empty() && !current_seq.is_empty() {
            records.push(SequenceRecord::new(
                current_id,
                current_desc,
                current_seq,
            ));
        }
        
        Ok(records)
    }
    
    fn format_name(&self) -> &str {
        "FASTA"
    }
}

/// FASTA format writer
#[derive(Debug, Clone)]
pub struct FastaWriter {
    /// Line width for sequence output
    line_width: usize,
    /// Buffer size for writing
    buffer_size: usize,
}

impl FastaWriter {
    /// Create a new FASTA writer with the default line width
    pub fn new() -> Self {
        Self {
            line_width: 60,
            buffer_size: 1024 * 1024, // 1MB
        }
    }
    
    /// Create a new FASTA writer with the specified line width
    pub fn with_line_width(line_width: usize) -> Self {
        Self {
            line_width,
            buffer_size: 1024 * 1024, // 1MB
        }
    }
    
    /// Set the buffer size
    pub fn with_buffer_size(mut self, buffer_size: usize) -> Self {
        self.buffer_size = buffer_size;
        self
    }
}

impl Default for FastaWriter {
    fn default() -> Self {
        Self::new()
    }
}

impl SequenceWriter for FastaWriter {
    fn write_file<P: AsRef<Path>>(&self, records: &[SequenceRecord], path: P) -> EngineResult<()> {
        let mut writer = FastWriter::new(path, Some(self.buffer_size))?;
        
        for record in records {
            // Write header
            let header = match &record.description {
                Some(desc) => format!(">{} {}\n", record.id, desc),
                None => format!(">{}\n", record.id),
            };
            writer.write(header.as_bytes())?;
            
            // Write sequence with line wrapping
            for chunk in record.sequence_as_vec().chunks(self.line_width) {
                writer.write(chunk)?;
                writer.write(b"\n")?;
            }
        }
        
        writer.flush()?;
        Ok(())
    }
    
    fn write_string(&self, records: &[SequenceRecord]) -> EngineResult<String> {
        let mut output = String::new();
        
        for record in records {
            // Write header
            match &record.description {
                Some(desc) => output.push_str(&format!(">{} {}\n", record.id, desc)),
                None => output.push_str(&format!(">{}\n", record.id)),
            };
            
            // Write sequence with line wrapping
            let sequence = record.sequence_as_vec();
            for i in (0..sequence.len()).step_by(self.line_width) {
                let end = (i + self.line_width).min(sequence.len());
                output.push_str(&String::from_utf8_lossy(&sequence[i..end]));
                output.push('\n');
            }
        }
        
        Ok(output)
    }
    
    fn format_name(&self) -> &str {
        "FASTA"
    }
}

/// FASTQ format parser
#[derive(Debug, Clone)]
pub struct FastqParser {
    /// Storage mode to use for sequences
    storage_mode: StorageMode,
    /// Buffer size for reading
    buffer_size: usize,
}

impl FastqParser {
    /// Create a new FASTQ parser with the default storage mode
    pub fn new() -> Self {
        Self {
            storage_mode: StorageMode::default(),
            buffer_size: 1024 * 1024, // 1MB
        }
    }
    
    /// Create a new FASTQ parser with the specified storage mode
    pub fn with_storage_mode(storage_mode: StorageMode) -> Self {
        Self {
            storage_mode,
            buffer_size: 1024 * 1024, // 1MB
        }
    }
    
    /// Set the buffer size
    pub fn with_buffer_size(mut self, buffer_size: usize) -> Self {
        self.buffer_size = buffer_size;
        self
    }
}

impl Default for FastqParser {
    fn default() -> Self {
        Self::new()
    }
}

impl SequenceParser for FastqParser {
    fn parse_file<P: AsRef<Path>>(&self, path: P) -> EngineResult<Vec<SequenceRecord>> {
        let mut reader = FastReader::new(path.as_ref(), Some(self.buffer_size))?;
        
        let mut records = Vec::new();
        let mut line_counter = 0;
        
        let mut current_id = String::new();
        let mut current_desc = None;
        let mut current_seq = Vec::new();
        let mut current_qual = Vec::new();
        
        for line_result in reader.read_lines() {
            let line = line_result?;
            let phase = line_counter % 4;
            
            match phase {
                0 => {
                    // Header line
                    if !line.starts_with('@') {
                        return Err(EngineError::InvalidSequenceData(
                            format!("Invalid FASTQ header: {}", line)
                        ));
                    }
                    
                    // Parse header
                    let header = &line[1..];
                    let parts: Vec<&str> = header.splitn(2, ' ').collect();
                    
                    current_id = parts[0].to_string();
                    current_desc = parts.get(1).map(|s| s.to_string());
                },
                1 => {
                    // Sequence line
                    current_seq = line.as_bytes().to_vec();
                },
                2 => {
                    // Separator line (should start with '+')
                    if !line.starts_with('+') {
                        return Err(EngineError::InvalidSequenceData(
                            format!("Invalid FASTQ separator: {}", line)
                        ));
                    }
                },
                3 => {
                    // Quality line
                    current_qual = line.as_bytes().to_vec();
                    
                    // Validate quality length
                    if current_qual.len() != current_seq.len() {
                        return Err(EngineError::InvalidSequenceData(
                            format!(
                                "Quality length ({}) does not match sequence length ({}) for record {}",
                                current_qual.len(), current_seq.len(), current_id
                            )
                        ));
                    }
                    
                    // Create sequence storages
                    let sequence = StorageFactory::create_storage(
                        Some(current_seq.clone()),
                        Some(path.as_ref()),
                        Some(current_seq.len()),
                        Some(self.storage_mode),
                    )?;
                    
                    let quality = StorageFactory::create_storage(
                        Some(current_qual.clone()),
                        Some(path.as_ref()),
                        Some(current_qual.len()),
                        Some(self.storage_mode),
                    )?;
                    
                    // Add the record
                    records.push(SequenceRecord {
                        id: current_id.clone(),
                        description: current_desc.clone(),
                        sequence,
                        quality: Some(quality),
                        metadata: HashMap::new(),
                    });
                },
                _ => unreachable!(),
            }
            
            line_counter += 1;
        }
        
        // Validate that we have complete records
        if line_counter % 4 != 0 {
            return Err(EngineError::InvalidSequenceData(
                "Incomplete FASTQ record at end of file".to_string()
            ));
        }
        
        Ok(records)
    }
    
    fn parse_string(&self, content: &str) -> EngineResult<Vec<SequenceRecord>> {
        let mut records = Vec::new();
        let mut lines = content.lines();
        
        loop {
            // Header line
            let header = match lines.next() {
                Some(line) => line,
                None => break,
            };
            
            if !header.starts_with('@') {
                return Err(EngineError::InvalidSequenceData(
                    format!("Invalid FASTQ header: {}", header)
                ));
            }
            
            // Parse header
            let header = &header[1..];
            let parts: Vec<&str> = header.splitn(2, ' ').collect();
            
            let id = parts[0].to_string();
            let desc = parts.get(1).map(|s| s.to_string());
            
            // Sequence line
            let seq = match lines.next() {
                Some(line) => line.as_bytes().to_vec(),
                None => return Err(EngineError::InvalidSequenceData(
                    "Incomplete FASTQ record (missing sequence)".to_string()
                )),
            };
            
            // Separator line
            let separator = match lines.next() {
                Some(line) => line,
                None => return Err(EngineError::InvalidSequenceData(
                    "Incomplete FASTQ record (missing separator)".to_string()
                )),
            };
            
            if !separator.starts_with('+') {
                return Err(EngineError::InvalidSequenceData(
                    format!("Invalid FASTQ separator: {}", separator)
                ));
            }
            
            // Quality line
            let qual = match lines.next() {
                Some(line) => line.as_bytes().to_vec(),
                None => return Err(EngineError::InvalidSequenceData(
                    "Incomplete FASTQ record (missing quality)".to_string()
                )),
            };
            
            // Validate quality length
            if qual.len() != seq.len() {
                return Err(EngineError::InvalidSequenceData(
                    format!(
                        "Quality length ({}) does not match sequence length ({}) for record {}",
                        qual.len(), seq.len(), id
                    )
                ));
            }
            
            // Add the record
            records.push(SequenceRecord::with_quality(
                id,
                desc,
                seq,
                qual,
            ));
        }
        
        Ok(records)
    }
    
    fn format_name(&self) -> &str {
        "FASTQ"
    }
}

/// FASTQ format writer
#[derive(Debug, Clone)]
pub struct FastqWriter {
    /// Buffer size for writing
    buffer_size: usize,
}

impl FastqWriter {
    /// Create a new FASTQ writer
    pub fn new() -> Self {
        Self {
            buffer_size: 1024 * 1024, // 1MB
        }
    }
    
    /// Set the buffer size
    pub fn with_buffer_size(mut self, buffer_size: usize) -> Self {
        self.buffer_size = buffer_size;
        self
    }
}

impl Default for FastqWriter {
    fn default() -> Self {
        Self::new()
    }
}

impl SequenceWriter for FastqWriter {
    fn write_file<P: AsRef<Path>>(&self, records: &[SequenceRecord], path: P) -> EngineResult<()> {
        let mut writer = FastWriter::new(path, Some(self.buffer_size))?;
        
        for record in records {
            // Check if record has quality scores
            let quality = match record.quality_as_vec() {
                Some(q) => q,
                None => return Err(EngineError::InvalidSequenceData(
                    format!("Record {} does not have quality scores (required for FASTQ)", record.id)
                )),
            };
            
            // Write header
            let header = match &record.description {
                Some(desc) => format!("@{} {}\n", record.id, desc),
                None => format!("@{}\n", record.id),
            };
            writer.write(header.as_bytes())?;
            
            // Write sequence
            writer.write(&record.sequence_as_vec())?;
            writer.write(b"\n")?;
            
            // Write separator
            writer.write(b"+\n")?;
            
            // Write quality
            writer.write(&quality)?;
            writer.write(b"\n")?;
        }
        
        writer.flush()?;
        Ok(())
    }
    
    fn write_string(&self, records: &[SequenceRecord]) -> EngineResult<String> {
        let mut output = String::new();
        
        for record in records {
            // Check if record has quality scores
            let quality = match record.quality_as_vec() {
                Some(q) => q,
                None => return Err(EngineError::InvalidSequenceData(
                    format!("Record {} does not have quality scores (required for FASTQ)", record.id)
                )),
            };
            
            // Write header
            match &record.description {
                Some(desc) => output.push_str(&format!("@{} {}\n", record.id, desc)),
                None => output.push_str(&format!("@{}\n", record.id)),
            };
            
            // Write sequence
            output.push_str(&String::from_utf8_lossy(&record.sequence_as_vec()));
            output.push('\n');
            
            // Write separator
            output.push_str("+\n");
            
            // Write quality
            output.push_str(&String::from_utf8_lossy(&quality));
            output.push('\n');
        }
        
        Ok(output)
    }
    
    fn format_name(&self) -> &str {
        "FASTQ"
    }
}

/// Detect the format of a sequence file based on its content
pub fn detect_format<P: AsRef<Path>>(path: P) -> EngineResult<&'static str> {
    let mut reader = FastReader::new(path.as_ref(), None)?;
    
    // Read the first line to determine the format
    if let Some(first_line) = reader.read_lines().next() {
        let line = first_line?;
        
        if line.starts_with('>') {
            return Ok("FASTA");
        } else if line.starts_with('@') {
            // Check the third line to confirm it's FASTQ
            for (i, line_result) in reader.read_lines().enumerate() {
                if i == 1 { // Third line (after we've already read the first)
                    let line = line_result?;
                    if line.starts_with('+') {
                        return Ok("FASTQ");
                    }
                    break;
                }
            }
        }
    }
    
    Err(EngineError::InvalidSequenceData(
        "Could not determine file format".to_string()
    ))
}

/// Create a parser for the specified format
pub fn create_parser(format: &str) -> EngineResult<Box<dyn SequenceParser>> {
    match format.to_uppercase().as_str() {
        "FASTA" => Ok(Box::new(FastaParser::new())),
        "FASTQ" => Ok(Box::new(FastqParser::new())),
        _ => Err(EngineError::UnsupportedOperation(
            format!("Unsupported format: {}", format)
        )),
    }
}

/// Create a writer for the specified format
pub fn create_writer(format: &str) -> EngineResult<Box<dyn SequenceWriter>> {
    match format.to_uppercase().as_str() {
        "FASTA" => Ok(Box::new(FastaWriter::new())),
        "FASTQ" => Ok(Box::new(FastqWriter::new())),
        _ => Err(EngineError::UnsupportedOperation(
            format!("Unsupported format: {}", format)
        )),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;
    use tempfile::tempdir;
    
    #[test]
    fn test_fasta_parsing() -> std::io::Result<()> {
        // Create a temporary FASTA file
        let dir = tempdir()?;
        let file_path = dir.path().join("test.fasta");
        
        let fasta_content = ">seq1 First sequence\nACGTACGT\n>seq2 Second sequence\nGTACGTAC\n";
        {
            let mut file = std::fs::File::create(&file_path)?;
            file.write_all(fasta_content.as_bytes())?;
        }
        
        // Parse the file
        let parser = FastaParser::new();
        let records = parser.parse_file(&file_path).unwrap();
        
        // Check the parsed records
        assert_eq!(records.len(), 2);
        assert_eq!(records[0].id, "seq1");
        assert_eq!(records[0].description, Some("First sequence".to_string()));
        assert_eq!(records[0].sequence_as_vec(), b"ACGTACGT");
        
        assert_eq!(records[1].id, "seq2");
        assert_eq!(records[1].description, Some("Second sequence".to_string()));
        assert_eq!(records[1].sequence_as_vec(), b"GTACGTAC");
        
        Ok(())
    }
    
    #[test]
    fn test_fasta_writing() -> std::io::Result<()> {
        // Create records
        let records = vec![
            SequenceRecord::new(
                "seq1".to_string(),
                Some("First sequence".to_string()),
                b"ACGTACGT".to_vec(),
            ),
            SequenceRecord::new(
                "seq2".to_string(),
                Some("Second sequence".to_string()),
                b"GTACGTAC".to_vec(),
            ),
        ];
        
        // Create a temporary file
        let dir = tempdir()?;
        let file_path = dir.path().join("output.fasta");
        
        // Write the records
        let writer = FastaWriter::new();
        writer.write_file(&records, &file_path).unwrap();
        
        // Read the file and check the content
        let content = std::fs::read_to_string(&file_path)?;
        
        assert!(content.contains(">seq1 First sequence"));
        assert!(content.contains("ACGTACGT"));
        assert!(content.contains(">seq2 Second sequence"));
        assert!(content.contains("GTACGTAC"));
        
        Ok(())
    }
    
    #[test]
    fn test_fastq_parsing() -> std::io::Result<()> {
        // Create a temporary FASTQ file
        let dir = tempdir()?;
        let file_path = dir.path().join("test.fastq");
        
        let fastq_content = "@seq1 First sequence\nACGT\n+\nHHHH\n@seq2 Second sequence\nGTAC\n+\nIIII\n";
        {
            let mut file = std::fs::File::create(&file_path)?;
            file.write_all(fastq_content.as_bytes())?;
        }
        
        // Parse the file
        let parser = FastqParser::new();
        let records = parser.parse_file(&file_path).unwrap();
        
        // Check the parsed records
        assert_eq!(records.len(), 2);
        assert_eq!(records[0].id, "seq1");
        assert_eq!(records[0].description, Some("First sequence".to_string()));
        assert_eq!(records[0].sequence_as_vec(), b"ACGT");
        assert_eq!(records[0].quality_as_vec().unwrap(), b"HHHH");
        
        assert_eq!(records[1].id, "seq2");
        assert_eq!(records[1].description, Some("Second sequence".to_string()));
        assert_eq!(records[1].sequence_as_vec(), b"GTAC");
        assert_eq!(records[1].quality_as_vec().unwrap(), b"IIII");
        
        Ok(())
    }
    
    #[test]
    fn test_format_detection() -> std::io::Result<()> {
        // Create temporary files
        let dir = tempdir()?;
        let fasta_path = dir.path().join("test.fasta");
        let fastq_path = dir.path().join("test.fastq");
        
        // Create a FASTA file
        {
            let mut file = std::fs::File::create(&fasta_path)?;
            file.write_all(b">seq1\nACGT\n")?;
        }
        
        // Create a FASTQ file
        {
            let mut file = std::fs::File::create(&fastq_path)?;
            file.write_all(b"@seq1\nACGT\n+\nHHHH\n")?;
        }
        
        // Detect formats
        assert_eq!(detect_format(&fasta_path).unwrap(), "FASTA");
        assert_eq!(detect_format(&fastq_path).unwrap(), "FASTQ");
        
        Ok(())
    }
    
    #[test]
    fn test_format_factory() {
        // Create parsers
        let fasta_parser = create_parser("FASTA").unwrap();
        let fastq_parser = create_parser("FASTQ").unwrap();
        
        assert_eq!(fasta_parser.format_name(), "FASTA");
        assert_eq!(fastq_parser.format_name(), "FASTQ");
        
        // Create writers
        let fasta_writer = create_writer("FASTA").unwrap();
        let fastq_writer = create_writer("FASTQ").unwrap();
        
        assert_eq!(fasta_writer.format_name(), "FASTA");
        assert_eq!(fastq_writer.format_name(), "FASTQ");
        
        // Test error case with unsupported format
        assert!(create_parser("UNKNOWN").is_err());
        assert!(create_writer("UNKNOWN").is_err());
    }
    
    #[test]
    fn test_fastq_writing() -> std::io::Result<()> {
        // Create records with quality scores
        let records = vec![
            SequenceRecord::with_quality(
                "seq1".to_string(),
                Some("First sequence".to_string()),
                b"ACGT".to_vec(),
                b"HHHH".to_vec(),
            ),
            SequenceRecord::with_quality(
                "seq2".to_string(),
                Some("Second sequence".to_string()),
                b"GTAC".to_vec(),
                b"IIII".to_vec(),
            ),
        ];
        
        // Create a temporary file
        let dir = tempdir()?;
        let file_path = dir.path().join("output.fastq");
        
        // Write the records
        let writer = FastqWriter::new();
        writer.write_file(&records, &file_path).unwrap();
        
        // Read the file and check the content
        let content = std::fs::read_to_string(&file_path)?;
        
        assert!(content.contains("@seq1 First sequence"));
        assert!(content.contains("ACGT"));
        assert!(content.contains("+"));
        assert!(content.contains("HHHH"));
        assert!(content.contains("@seq2 Second sequence"));
        assert!(content.contains("GTAC"));
        assert!(content.contains("IIII"));
        
        Ok(())
    }
    
    #[test]
    fn test_sequence_record_methods() {
        // Create a record
        let mut record = SequenceRecord::new(
            "test".to_string(),
            Some("Test sequence".to_string()),
            b"ACGTACGT".to_vec(),
        );
        
        // Test basic properties
        assert_eq!(record.id, "test");
        assert_eq!(record.description, Some("Test sequence".to_string()));
        assert_eq!(record.len(), 8);
        assert!(!record.is_empty());
        
        // Test subsequence
        let sub = record.subsequence(2, 6);
        assert_eq!(sub, b"GTAC");
        
        // Test metadata
        record.add_metadata("source", "test data");
        record.add_metadata("date", "2023-01-01");
        
        assert_eq!(record.get_metadata("source"), Some(&"test data".to_string()));
        assert_eq!(record.get_metadata("date"), Some(&"2023-01-01".to_string()));
        assert_eq!(record.get_metadata("missing"), None);
    }