//! Core sequence type
//!
//! This module provides the core sequence types and operations for bioinformatics.

use std::fmt;
use std::ops::{Index, Range, RangeBounds}

/// A view into a sequence
pub struct SequenceView<'a> {
    sequence: &'a Sequence,
    start: usize,
    end: usize,
}

impl<'a> SequenceView<'a> {
    /// Create a new sequence view
    pub fn new(sequence: &'a Sequence, start: usize, end: usize) -> SequenceResult<Self> {
        if start > end || end > sequence.len() {
            return Err(SequenceError::IndexOutOfBounds(
                format!("Invalid range {}..{} for sequence of length {}", start, end, sequence.len())
            ));
        }
        
        Ok(Self {
            sequence,
            start,
            end,
        })
    }
    
    /// Get the length of the view
    pub fn len(&self) -> usize {
        self.end - self.start
    }
    
    /// Check if the view is empty
    pub fn is_empty(&self) -> bool {
        self.start == self.end
    }
    
    /// Get the sequence as a string
    pub fn as_string(&self) -> String {
        String::from_utf8_lossy(&self.as_bytes()).to_string()
    }
    
    /// Get the sequence as bytes
    pub fn as_bytes(&self) -> Vec<u8> {
        self.sequence.data.subsequence(self.start, self.end)
    }
    
    /// Slide the view to a new position
    pub fn slide(&self, offset: isize) -> SequenceResult<Self> {
        let new_start = if offset >= 0 {
            self.start.saturating_add(offset as usize)
        } else {
            self.start.saturating_sub((-offset) as usize)
        };
        
        let new_end = new_start + self.len();
        
        if new_end > self.sequence.len() {
            return Err(SequenceError::IndexOutOfBounds(
                format!("Sliding by {} would exceed sequence bounds", offset)
            ));
        }
        
        Ok(Self {
            sequence: self.sequence,
            start: new_start,
            end: new_end,
        })
    }
    
    /// Resize the view
    pub fn resize(&self, new_length: usize) -> SequenceResult<Self> {
        let new_end = self.start + new_length;
        
        if new_end > self.sequence.len() {
            return Err(SequenceError::IndexOutOfBounds(
                format!("Resizing to length {} would exceed sequence bounds", new_length)
            ));
        }
        
        Ok(Self {
            sequence: self.sequence,
            start: self.start,
            end: new_end,
        })
    }
    
    /// Convert the view to a full sequence
    pub fn to_sequence(&self) -> Sequence {
        Sequence {
            data: Box::new(InMemoryStorage::new(self.as_bytes())),
            alphabet: self.sequence.alphabet.clone(),
            id: self.sequence.id.clone(),
            description: self.sequence.description.clone().map(|desc| 
                format!("{} (view {}..{})", desc, self.start, self.end)
            ),
        }
    }
    
    /// Find all occurrences of a pattern in the view
    pub fn find_all(&self, pattern: &[u8]) -> Vec<usize> {
        if pattern.is_empty() || pattern.len() > self.len() {
            return Vec::new();
        }
        
        // Use the KMP algorithm for searching
        match string_ops::kmp_search(&self.as_bytes(), pattern) {
            Ok(matches) => matches,
            Err(_) => Vec::new(),
        }
    }
}

impl fmt::Display for Sequence {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        // Format ID and description if available
        if let Some(id) = &self.id {
            write!(f, ">{}", id)?;
            
            if let Some(desc) = &self.description {
                write!(f, " {}", desc)?;
            }
            
            writeln!(f)?;
        }
        
        // Format the sequence in lines of 60 characters
        let seq_str = self.as_string();
        for chunk in seq_str.as_bytes().chunks(60) {
            writeln!(f, "{}", std::str::from_utf8(chunk).unwrap_or("invalid UTF-8"))?;
        }
        
        Ok(())
    }
}

impl Index<usize> for Sequence {
    type Output = u8;
    
    fn index(&self, index: usize) -> &Self::Output {
        if let Some(slice) = self.data.as_slice() {
            &slice[index]
        } else {
            panic!("Index out of bounds: {}", index)
        }
    }
}

impl Index<Range<usize>> for Sequence {
    type Output = [u8];
    
    fn index(&self, range: Range<usize>) -> &Self::Output {
        if let Some(slice) = self.data.as_slice() {
            &slice[range]
        } else {
            panic!("Range out of bounds: {:?}", range)
        }
    }
}

impl<'a> Index<usize> for SequenceView<'a> {
    type Output = u8;
    
    fn index(&self, index: usize) -> &Self::Output {
        if index >= self.len() {
            panic!("Index out of bounds: {}", index);
        }
        
        if let Some(slice) = self.sequence.data.as_slice() {
            &slice[self.start + index]
        } else {
            panic!("Cannot access slice");
        }
    }
}

/// Get a sequence from string
impl From<&str> for Sequence {
    fn from(s: &str) -> Self {
        Self::new(s.as_bytes()).unwrap_or_else(|_| {
            // Default to DNA if we can't detect the alphabet
            Self {
                data: Box::new(InMemoryStorage::new(s.as_bytes().to_vec())),
                alphabet: Box::new(DNAAlphabet::default()),
                id: None,
                description: None,
            }
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    
    #[test]
    fn test_sequence_creation() {
        // Create a DNA sequence
        let dna = Sequence::new_dna(b"ACGTACGT").unwrap();
        assert_eq!(dna.len(), 8);
        assert_eq!(dna.alphabet_name(), "DNA");
        
        // Create an RNA sequence
        let rna = Sequence::new_rna(b"ACGUACGU").unwrap();
        assert_eq!(rna.len(), 8);
        assert_eq!(rna.alphabet_name(), "RNA");
        
        // Create a protein sequence
        let protein = Sequence::new_protein(b"ACDEFGHIKLMNPQRSTVWYX").unwrap();
        assert_eq!(protein.len(), 21);
        assert_eq!(protein.alphabet_name(), "Protein");
        
        // Auto-detect DNA
        let auto_dna = Sequence::new(b"ACGTACGT").unwrap();
        assert_eq!(auto_dna.alphabet_name(), "DNA");
        
        // Test with ID and description
        let seq_with_meta = Sequence::new_dna(b"ACGTACGT").unwrap()
            .with_id("seq1")
            .with_description("Test sequence");
        
        assert_eq!(seq_with_meta.id(), Some("seq1"));
        assert_eq!(seq_with_meta.description(), Some("Test sequence"));
    }
    
    #[test]
    fn test_sequence_operations() {
        // Create a DNA sequence
        let dna = Sequence::new_dna(b"ACGTACGT").unwrap();
        
        // Test reverse
        let reversed = dna.reverse();
        assert_eq!(reversed.as_bytes().as_ref(), b"TGCATGCA");
        
        // Test complement
        let complemented = dna.complement().unwrap();
        assert_eq!(complemented.as_bytes().as_ref(), b"TGCATGCA");
        
        // Test reverse complement
        let rev_comp = dna.reverse_complement().unwrap();
        assert_eq!(rev_comp.as_bytes().as_ref(), b"ACGTACGT");
        
        // Test transcription
        let rna = dna.transcribe().unwrap();
        assert_eq!(rna.as_bytes().as_ref(), b"ACGUACGU");
        assert_eq!(rna.alphabet_name(), "RNA");
        
        // Test concatenation
        let dna2 = Sequence::new_dna(b"TGCATGCA").unwrap();
        let concatenated = dna.concatenate(&dna2).unwrap();
        assert_eq!(concatenated.as_bytes().as_ref(), b"ACGTACGTTGCATGCA");
        
        // Test subsequence
        let subseq = dna.subsequence(2, 6).unwrap();
        assert_eq!(subseq.as_bytes().as_ref(), b"GTAC");
        
        // Test GC content
        let gc = dna.gc_content().unwrap();
        assert_eq!(gc, 50.0); // 4/8 = 50%
        
        // Test find all
        let positions = dna.find_all(b"AC");
        assert_eq!(positions, vec![0, 4]);
        
        // Test count
        let count = dna.count(b"AC");
        assert_eq!(count, 2);
    }
    
    #[test]
    fn test_sequence_view() {
        // Create a DNA sequence
        let dna = Sequence::new_dna(b"ACGTACGTACGT").unwrap();
        
        // Create a view
        let view = dna.view().subsequence(2, 10).unwrap();
        assert_eq!(view.as_bytes().as_ref(), b"GTACGTAC");
        assert_eq!(view.len(), 8);
        
        // Test slide
        let slid = view.slide(2).unwrap();
        assert_eq!(slid.as_bytes().as_ref(), b"TACGTACG");
        
        // Test resize
        let resized = view.resize(4).unwrap();
        assert_eq!(resized.as_bytes().as_ref(), b"GTAC");
        
        // Test to_sequence
        let new_seq = view.to_sequence();
        assert_eq!(new_seq.as_bytes().as_ref(), b"GTACGTAC");
    }
};
use std::borrow::Cow;
use thiserror::Error;

use crate::engines::core::memory::PackedDnaStorage;
use crate::engines::storage::{StorableSequence, InMemoryStorage};
use crate::engines::compute::string_ops;
use super::alphabet::{Alphabet, DNAAlphabet, RNAAlphabet, ProteinAlphabet};

/// Error type for sequence operations
#[derive(Error, Debug)]
pub enum SequenceError {
    #[error("Invalid sequence: {0}")]
    InvalidSequence(String),
    
    #[error("Invalid alphabet: {0}")]
    InvalidAlphabet(String),
    
    #[error("Index out of bounds: {0}")]
    IndexOutOfBounds(String),
    
    #[error("Operation not supported: {0}")]
    UnsupportedOperation(String),
    
    #[error("Engine error: {0}")]
    EngineError(#[from] crate::engines::EngineError),
}

/// Result type for sequence operations
pub type SequenceResult<T> = Result<T, SequenceError>;

/// Common sequence type for all biological sequences
#[derive(Clone)]
pub struct Sequence {
    /// The sequence data
    data: Box<dyn StorableSequence>,
    /// The alphabet used for this sequence
    alphabet: Box<dyn Alphabet>,
    /// Identifier for the sequence (optional)
    id: Option<String>,
    /// Description of the sequence (optional)
    description: Option<String>,
}

impl Sequence {
    /// Create a new sequence from raw bytes
    pub fn new(data: &[u8]) -> SequenceResult<Self> {
        // Detect alphabet
        let alphabet = super::alphabet::detect_alphabet(data)
            .ok_or_else(|| SequenceError::InvalidSequence(
                "Could not detect alphabet for sequence".to_string()
            ))?;
        
        Ok(Self {
            data: Box::new(InMemoryStorage::new(data.to_vec())),
            alphabet,
            id: None,
            description: None,
        })
    }
    
    /// Create a new sequence with a specific alphabet
    pub fn with_alphabet<A: Alphabet + 'static>(data: &[u8], alphabet: A) -> SequenceResult<Self> {
        // Validate sequence against alphabet
        if !alphabet.is_valid_sequence(data) {
            return Err(SequenceError::InvalidSequence(
                format!("Sequence contains invalid characters for {} alphabet", alphabet.name())
            ));
        }
        
        Ok(Self {
            data: Box::new(InMemoryStorage::new(data.to_vec())),
            alphabet: Box::new(alphabet),
            id: None,
            description: None,
        })
    }
    
    /// Create a new DNA sequence
    pub fn new_dna(data: &[u8]) -> SequenceResult<Self> {
        Self::with_alphabet(data, DNAAlphabet::default())
    }
    
    /// Create a new RNA sequence
    pub fn new_rna(data: &[u8]) -> SequenceResult<Self> {
        Self::with_alphabet(data, RNAAlphabet::default())
    }
    
    /// Create a new protein sequence
    pub fn new_protein(data: &[u8]) -> SequenceResult<Self> {
        Self::with_alphabet(data, ProteinAlphabet::default())
    }
    
    /// Set the sequence identifier
    pub fn with_id(mut self, id: &str) -> Self {
        self.id = Some(id.to_string());
        self
    }
    
    /// Set the sequence description
    pub fn with_description(mut self, description: &str) -> Self {
        self.description = Some(description.to_string());
        self
    }
    
    /// Get the sequence length
    pub fn len(&self) -> usize {
        self.data.len()
    }
    
    /// Check if the sequence is empty
    pub fn is_empty(&self) -> bool {
        self.data.len() == 0
    }
    
    /// Get the sequence as a string
    pub fn as_string(&self) -> String {
        String::from_utf8_lossy(self.as_bytes()).to_string()
    }
    
    /// Get the sequence as bytes
    pub fn as_bytes(&self) -> Cow<[u8]> {
        if let Some(slice) = self.data.as_slice() {
            Cow::Borrowed(slice)
        } else {
            Cow::Owned(self.data.subsequence(0, self.data.len()))
        }
    }
    
    /// Get a subsequence
    pub fn subsequence(&self, start: usize, end: usize) -> SequenceResult<Self> {
        if start > end || end > self.len() {
            return Err(SequenceError::IndexOutOfBounds(
                format!("Invalid range {}..{} for sequence of length {}", start, end, self.len())
            ));
        }
        
        let subseq = self.data.subsequence(start, end);
        
        Ok(Self {
            data: Box::new(InMemoryStorage::new(subseq)),
            alphabet: self.alphabet.clone(),
            id: self.id.clone(),
            description: self.description.clone().map(|desc| format!("{} (subsequence {}..{})", desc, start, end)),
        })
    }
    
    /// Get a view of the sequence
    pub fn view(&self) -> SequenceView {
        SequenceView {
            sequence: self,
            start: 0,
            end: self.len(),
        }
    }
    
    /// Get the alphabet name
    pub fn alphabet_name(&self) -> &str {
        self.alphabet.name()
    }
    
    /// Get the identifier (if any)
    pub fn id(&self) -> Option<&str> {
        self.id.as_deref()
    }
    
    /// Get the description (if any)
    pub fn description(&self) -> Option<&str> {
        self.description.as_deref()
    }
    
    /// Get the base composition
    pub fn base_composition(&self) -> std::collections::HashMap<u8, usize> {
        let mut counts = std::collections::HashMap::new();
        
        // Count occurrences of each base
        for &base in self.as_bytes().iter() {
            *counts.entry(base).or_insert(0) += 1;
        }
        
        counts
    }
    
    /// Get the GC content (for DNA/RNA sequences)
    pub fn gc_content(&self) -> SequenceResult<f64> {
        if self.alphabet_name() != "DNA" && self.alphabet_name() != "RNA" {
            return Err(SequenceError::UnsupportedOperation(
                format!("GC content calculation not supported for {} alphabet", self.alphabet_name())
            ));
        }
        
        let composition = self.base_composition();
        let total = self.len() as f64;
        
        if total == 0.0 {
            return Ok(0.0);
        }
        
        // Count G and C bases (both upper and lowercase)
        let gc_count = 
            composition.get(&b'G').unwrap_or(&0) +
            composition.get(&b'g').unwrap_or(&0) +
            composition.get(&b'C').unwrap_or(&0) +
            composition.get(&b'c').unwrap_or(&0);
        
        Ok((gc_count as f64) / total * 100.0)
    }
    
    /// Get the reverse of the sequence
    pub fn reverse(&self) -> Self {
        let mut reversed = self.as_bytes().to_vec();
        string_ops::reverse_in_place(&mut reversed);
        
        Self {
            data: Box::new(InMemoryStorage::new(reversed)),
            alphabet: self.alphabet.clone(),
            id: self.id.clone(),
            description: self.description.clone().map(|desc| format!("{} (reversed)", desc)),
        }
    }
    
    /// Get the complement of the sequence (for DNA/RNA)
    pub fn complement(&self) -> SequenceResult<Self> {
        if self.alphabet_name() != "DNA" && self.alphabet_name() != "RNA" {
            return Err(SequenceError::UnsupportedOperation(
                format!("Complement operation not supported for {} alphabet", self.alphabet_name())
            ));
        }
        
        let seq_bytes = self.as_bytes();
        let complemented = self.alphabet.complement_sequence(&seq_bytes)
            .ok_or_else(|| SequenceError::UnsupportedOperation(
                "Failed to compute complement".to_string()
            ))?;
        
        Ok(Self {
            data: Box::new(InMemoryStorage::new(complemented)),
            alphabet: self.alphabet.clone(),
            id: self.id.clone(),
            description: self.description.clone().map(|desc| format!("{} (complement)", desc)),
        })
    }
    
    /// Get the reverse complement of the sequence (for DNA/RNA)
    pub fn reverse_complement(&self) -> SequenceResult<Self> {
        let mut rev_comp = self.as_bytes().to_vec();
        
        // First complement
        let complemented = self.alphabet.complement_sequence(&rev_comp)
            .ok_or_else(|| SequenceError::UnsupportedOperation(
                "Failed to compute complement".to_string()
            ))?;
        
        // Then reverse
        let mut reversed = complemented;
        string_ops::reverse_in_place(&mut reversed);
        
        Ok(Self {
            data: Box::new(InMemoryStorage::new(reversed)),
            alphabet: self.alphabet.clone(),
            id: self.id.clone(),
            description: self.description.clone().map(|desc| format!("{} (reverse complement)", desc)),
        })
    }
    
    /// Transcribe a DNA sequence to RNA
    pub fn transcribe(&self) -> SequenceResult<Self> {
        if self.alphabet_name() != "DNA" {
            return Err(SequenceError::UnsupportedOperation(
                "Transcription operation only supported for DNA alphabet".to_string()
            ));
        }
        
        let dna = self.as_bytes();
        let rna = string_ops::transcribe(&dna);
        
        Ok(Self {
            data: Box::new(InMemoryStorage::new(rna)),
            alphabet: Box::new(RNAAlphabet::default()),
            id: self.id.clone(),
            description: self.description.clone().map(|desc| format!("{} (transcribed)", desc)),
        })
    }
    
    /// Find all occurrences of a subsequence
    pub fn find_all(&self, pattern: &[u8]) -> Vec<usize> {
        if pattern.is_empty() || pattern.len() > self.len() {
            return Vec::new();
        }
        
        // Use the KMP algorithm for searching
        match string_ops::kmp_search(self.as_bytes().as_ref(), pattern) {
            Ok(matches) => matches,
            Err(_) => Vec::new(),
        }
    }
    
    /// Count the occurrences of a subsequence
    pub fn count(&self, pattern: &[u8]) -> usize {
        self.find_all(pattern).len()
    }
    
    /// Convert to a specific storage format
    pub fn to_packed_storage(&self) -> SequenceResult<Self> {
        if self.alphabet_name() != "DNA" {
            return Err(SequenceError::UnsupportedOperation(
                "Packed storage is only supported for DNA sequences".to_string()
            ));
        }
        
        let seq_data = self.as_bytes();
        let mut packed = PackedDnaStorage::with_capacity(seq_data.len());
        packed.pack(&seq_data);
        
        // Create a sequence with the packed storage
        Ok(Self {
            data: Box::new(InMemoryStorage::new(seq_data.to_vec())), // We'd use packed storage here in a real implementation
            alphabet: self.alphabet.clone(),
            id: self.id.clone(),
            description: self.description.clone(),
        })
    }
    
    /// Create a mask from the given positions
    pub fn mask(&self, positions: &[usize], mask_char: u8) -> SequenceResult<Self> {
        let mut masked = self.as_bytes().to_vec();
        
        for &pos in positions {
            if pos >= self.len() {
                return Err(SequenceError::IndexOutOfBounds(
                    format!("Position {} is out of bounds for sequence of length {}", pos, self.len())
                ));
            }
            
            masked[pos] = mask_char;
        }
        
        Ok(Self {
            data: Box::new(InMemoryStorage::new(masked)),
            alphabet: self.alphabet.clone(),
            id: self.id.clone(),
            description: self.description.clone().map(|desc| format!("{} (masked)", desc)),
        })
    }
    
    /// Concatenate with another sequence
    pub fn concatenate(&self, other: &Self) -> SequenceResult<Self> {
        if self.alphabet_name() != other.alphabet_name() {
            return Err(SequenceError::InvalidAlphabet(
                format!("Cannot concatenate sequences with different alphabets: {} and {}", 
                        self.alphabet_name(), other.alphabet_name())
            ));
        }
        
        // Combine the sequences
        let mut combined = self.as_bytes().to_vec();
        combined.extend_from_slice(&other.as_bytes());
        
        Ok(Self {
            data: Box::new(InMemoryStorage::new(combined)),
            alphabet: self.alphabet.clone(),
            id: self.id.clone().or_else(|| other.id.clone()),
            description: match (self.description.clone(), other.description.clone()) {
                (Some(desc1), Some(desc2)) => Some(format!("{} + {}", desc1, desc2)),
                (Some(desc1), None) => Some(format!("{} + [unnamed]", desc1)),
                (None, Some(desc2)) => Some(format!("[unnamed] + {}", desc2)),
                (None, None) => None,
            },
        })
    }
}