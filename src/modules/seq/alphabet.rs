//! Sequence alphabets
//!
//! This module provides definitions and validation for various
//! biological sequence alphabets.

use std::collections::HashSet;

/// Alphabet trait for sequence validation
pub trait Alphabet: Send + Sync {
    /// Get the name of the alphabet
    fn name(&self) -> &str;
    
    /// Check if a character is valid in this alphabet
    fn is_valid_char(&self, c: u8) -> bool;
    
    /// Check if a sequence is valid in this alphabet
    fn is_valid_sequence(&self, seq: &[u8]) -> bool {
        seq.iter().all(|&c| self.is_valid_char(c))
    }
    
    /// Get all valid characters in this alphabet
    fn valid_chars(&self) -> &[u8];
    
    /// Get the number of valid characters in this alphabet
    fn size(&self) -> usize {
        self.valid_chars().len()
    }
    
    /// Sanitize a sequence by replacing invalid characters
    fn sanitize(&self, seq: &[u8], replacement: u8) -> Vec<u8> {
        seq.iter()
            .map(|&c| if self.is_valid_char(c) { c } else { replacement })
            .collect()
    }
    
    /// Get the complement of a character (if applicable)
    fn complement(&self, c: u8) -> Option<u8>;
    
    /// Get the complement of a sequence (if applicable)
    fn complement_sequence(&self, seq: &[u8]) -> Option<Vec<u8>> {
        let mut result = Vec::with_capacity(seq.len());
        
        for &c in seq {
            if let Some(comp) = self.complement(c) {
                result.push(comp);
            } else {
                return None;
            }
        }
        
        Some(result)
    }
}

/// DNA alphabet (A, C, G, T, N and lowercase)
#[derive(Debug, Clone)]
pub struct DNAAlphabet {
    valid_chars: Vec<u8>,
    valid_set: HashSet<u8>,
    complement_map: [u8; 256],
}

impl Default for DNAAlphabet {
    fn default() -> Self {
        let mut obj = Self {
            valid_chars: b"ACGTNacgtn".to_vec(),
            valid_set: HashSet::from([b'A', b'C', b'G', b'T', b'N', b'a', b'c', b'g', b't', b'n']),
            complement_map: [0; 256],
        };
        
        // Initialize complement map
        for i in 0..256 {
            obj.complement_map[i] = i as u8;
        }
        
        // Set up complements
        obj.complement_map[b'A' as usize] = b'T';
        obj.complement_map[b'C' as usize] = b'G';
        obj.complement_map[b'G' as usize] = b'C';
        obj.complement_map[b'T' as usize] = b'A';
        obj.complement_map[b'N' as usize] = b'N';
        
        obj.complement_map[b'a' as usize] = b't';
        obj.complement_map[b'c' as usize] = b'g';
        obj.complement_map[b'g' as usize] = b'c';
        obj.complement_map[b't' as usize] = b'a';
        obj.complement_map[b'n' as usize] = b'n';
        
        obj
    }
}

impl Alphabet for DNAAlphabet {
    fn name(&self) -> &str {
        "DNA"
    }
    
    fn is_valid_char(&self, c: u8) -> bool {
        self.valid_set.contains(&c)
    }
    
    fn valid_chars(&self) -> &[u8] {
        &self.valid_chars
    }
    
    fn complement(&self, c: u8) -> Option<u8> {
        if self.is_valid_char(c) {
            Some(self.complement_map[c as usize])
        } else {
            None
        }
    }
}

/// RNA alphabet (A, C, G, U, N and lowercase)
#[derive(Debug, Clone)]
pub struct RNAAlphabet {
    valid_chars: Vec<u8>,
    valid_set: HashSet<u8>,
    complement_map: [u8; 256],
}

impl Default for RNAAlphabet {
    fn default() -> Self {
        let mut obj = Self {
            valid_chars: b"ACGUNacgun".to_vec(),
            valid_set: HashSet::from([b'A', b'C', b'G', b'U', b'N', b'a', b'c', b'g', b'u', b'n']),
            complement_map: [0; 256],
        };
        
        // Initialize complement map
        for i in 0..256 {
            obj.complement_map[i] = i as u8;
        }
        
        // Set up complements
        obj.complement_map[b'A' as usize] = b'U';
        obj.complement_map[b'C' as usize] = b'G';
        obj.complement_map[b'G' as usize] = b'C';
        obj.complement_map[b'U' as usize] = b'A';
        obj.complement_map[b'N' as usize] = b'N';
        
        obj.complement_map[b'a' as usize] = b'u';
        obj.complement_map[b'c' as usize] = b'g';
        obj.complement_map[b'g' as usize] = b'c';
        obj.complement_map[b'u' as usize] = b'a';
        obj.complement_map[b'n' as usize] = b'n';
        
        obj
    }
}

impl Alphabet for RNAAlphabet {
    fn name(&self) -> &str {
        "RNA"
    }
    
    fn is_valid_char(&self, c: u8) -> bool {
        self.valid_set.contains(&c)
    }
    
    fn valid_chars(&self) -> &[u8] {
        &self.valid_chars
    }
    
    fn complement(&self, c: u8) -> Option<u8> {
        if self.is_valid_char(c) {
            Some(self.complement_map[c as usize])
        } else {
            None
        }
    }
}

/// Protein alphabet (standard amino acids and X for unknown)
#[derive(Debug, Clone)]
pub struct ProteinAlphabet {
    valid_chars: Vec<u8>,
    valid_set: HashSet<u8>,
}

impl Default for ProteinAlphabet {
    fn default() -> Self {
        Self {
            valid_chars: b"ACDEFGHIKLMNPQRSTVWYXacdefghiklmnpqrstvwyx*".to_vec(),
            valid_set: HashSet::from([
                b'A', b'C', b'D', b'E', b'F', b'G', b'H', b'I', b'K', b'L',
                b'M', b'N', b'P', b'Q', b'R', b'S', b'T', b'V', b'W', b'Y', b'X',
                b'a', b'c', b'd', b'e', b'f', b'g', b'h', b'i', b'k', b'l',
                b'm', b'n', b'p', b'q', b'r', b's', b't', b'v', b'w', b'y', b'x',
                b'*',
            ]),
        }
    }
}

impl Alphabet for ProteinAlphabet {
    fn name(&self) -> &str {
        "Protein"
    }
    
    fn is_valid_char(&self, c: u8) -> bool {
        self.valid_set.contains(&c)
    }
    
    fn valid_chars(&self) -> &[u8] {
        &self.valid_chars
    }
    
    fn complement(&self, _c: u8) -> Option<u8> {
        None // Proteins don't have complements
    }
}

/// Detect the alphabet of a sequence
pub fn detect_alphabet(seq: &[u8]) -> Option<Box<dyn Alphabet>> {
    // Check for DNA
    let dna_alphabet = DNAAlphabet::default();
    if dna_alphabet.is_valid_sequence(seq) {
        return Some(Box::new(dna_alphabet));
    }
    
    // Check for RNA
    let rna_alphabet = RNAAlphabet::default();
    if rna_alphabet.is_valid_sequence(seq) {
        return Some(Box::new(rna_alphabet));
    }
    
    // Check for protein
    let protein_alphabet = ProteinAlphabet::default();
    if protein_alphabet.is_valid_sequence(seq) {
        return Some(Box::new(protein_alphabet));
    }
    
    None
}

/// Convert a DNA sequence to RNA
pub fn dna_to_rna(seq: &[u8]) -> Vec<u8> {
    seq.iter()
        .map(|&c| match c {
            b'T' => b'U',
            b't' => b'u',
            _ => c,
        })
        .collect()
}

/// Convert an RNA sequence to DNA
pub fn rna_to_dna(seq: &[u8]) -> Vec<u8> {
    seq.iter()
        .map(|&c| match c {
            b'U' => b'T',
            b'u' => b't',
            _ => c,
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    
    #[test]
    fn test_dna_alphabet() {
        let alphabet = DNAAlphabet::default();
        
        // Test valid characters
        assert!(alphabet.is_valid_char(b'A'));
        assert!(alphabet.is_valid_char(b'C'));
        assert!(alphabet.is_valid_char(b'G'));
        assert!(alphabet.is_valid_char(b'T'));
        assert!(alphabet.is_valid_char(b'N'));
        assert!(alphabet.is_valid_char(b'a'));
        assert!(alphabet.is_valid_char(b'c'));
        assert!(alphabet.is_valid_char(b'g'));
        assert!(alphabet.is_valid_char(b't'));
        assert!(alphabet.is_valid_char(b'n'));
        
        // Test invalid characters
        assert!(!alphabet.is_valid_char(b'U'));
        assert!(!alphabet.is_valid_char(b'X'));
        assert!(!alphabet.is_valid_char(b'Z'));
        
        // Test valid sequence
        assert!(alphabet.is_valid_sequence(b"ACGTACGT"));
        assert!(alphabet.is_valid_sequence(b"acgtACGT"));
        
        // Test invalid sequence
        assert!(!alphabet.is_valid_sequence(b"ACGUXCGT"));
        
        // Test complement
        assert_eq!(alphabet.complement(b'A'), Some(b'T'));
        assert_eq!(alphabet.complement(b'C'), Some(b'G'));
        assert_eq!(alphabet.complement(b'G'), Some(b'C'));
        assert_eq!(alphabet.complement(b'T'), Some(b'A'));
        assert_eq!(alphabet.complement(b'N'), Some(b'N'));
        
        // Test complement sequence
        assert_eq!(
            alphabet.complement_sequence(b"ACGT"),
            Some(b"TGCA".to_vec())
        );
        
        // Test sanitize
        assert_eq!(
            alphabet.sanitize(b"ACGUX", b'N'),
            b"ACGNN".to_vec()
        );
    }
    
    #[test]
    fn test_rna_alphabet() {
        let alphabet = RNAAlphabet::default();
        
        // Test valid characters
        assert!(alphabet.is_valid_char(b'A'));
        assert!(alphabet.is_valid_char(b'C'));
        assert!(alphabet.is_valid_char(b'G'));
        assert!(alphabet.is_valid_char(b'U'));
        assert!(alphabet.is_valid_char(b'N'));
        
        // Test invalid characters
        assert!(!alphabet.is_valid_char(b'T'));
        assert!(!alphabet.is_valid_char(b'X'));
        
        // Test valid sequence
        assert!(alphabet.is_valid_sequence(b"ACGUACGU"));
        
        // Test invalid sequence
        assert!(!alphabet.is_valid_sequence(b"ACGUTCGU"));
        
        // Test complement
        assert_eq!(alphabet.complement(b'A'), Some(b'U'));
        assert_eq!(alphabet.complement(b'C'), Some(b'G'));
        assert_eq!(alphabet.complement(b'G'), Some(b'C'));
        assert_eq!(alphabet.complement(b'U'), Some(b'A'));
    }
    
    #[test]
    fn test_protein_alphabet() {
        let alphabet = ProteinAlphabet::default();
        
        // Test valid characters
        assert!(alphabet.is_valid_char(b'A'));
        assert!(alphabet.is_valid_char(b'C'));
        assert!(alphabet.is_valid_char(b'D'));
        assert!(alphabet.is_valid_char(b'X'));
        assert!(alphabet.is_valid_char(b'*'));
        
        // Test invalid characters
        assert!(!alphabet.is_valid_char(b'J'));
        assert!(!alphabet.is_valid_char(b'O'));
        assert!(!alphabet.is_valid_char(b'U'));
        assert!(!alphabet.is_valid_char(b'Z'));
        
        // Test valid sequence
        assert!(alphabet.is_valid_sequence(b"ACDEFGHIKLMNPQRSTVWYX"));
        
        // Test invalid sequence
        assert!(!alphabet.is_valid_sequence(b"ACDJFGHIKLMNPQRSTVWYX"));
        
        // Test complement (should be None for proteins)
        assert_eq!(alphabet.complement(b'A'), None);
    }
    
    #[test]
    fn test_detect_alphabet() {
        // Test DNA detection
        let dna_seq = b"ACGTACGT";
        let alphabet = detect_alphabet(dna_seq).unwrap();
        assert_eq!(alphabet.name(), "DNA");
        
        // Test RNA detection
        let rna_seq = b"ACGUACGU";
        let alphabet = detect_alphabet(rna_seq).unwrap();
        assert_eq!(alphabet.name(), "RNA");
        
        // Test protein detection
        let protein_seq = b"ACDEFGHIKLMNPQRSTVWYX";
        let alphabet = detect_alphabet(protein_seq).unwrap();
        assert_eq!(alphabet.name(), "Protein");
        
        // Test unknown sequence
        let unknown_seq = b"ACGTJ123";
        assert!(detect_alphabet(unknown_seq).is_none());
    }
    
    #[test]
    fn test_dna_rna_conversion() {
        // Test DNA to RNA
        assert_eq!(dna_to_rna(b"ACGT"), b"ACGU");
        assert_eq!(dna_to_rna(b"acgt"), b"acgu");
        
        // Test RNA to DNA
        assert_eq!(rna_to_dna(b"ACGU"), b"ACGT");
        assert_eq!(rna_to_dna(b"acgu"), b"acgt");
    }
}