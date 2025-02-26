//! High-performance string operations for sequence data
//!
//! This module provides optimized implementations of common string operations
//! used in biological sequence processing.

use super::{ComputeResult, ComputeError};
use crate::engines::core::simd;
use std::cmp;

/// Perform a substring search with the Knuth-Morris-Pratt algorithm
///
/// KMP is efficient for searching for occurrences of a pattern within a text,
/// with linear time complexity O(n + m) where n is the text length and m is the pattern length.
pub fn kmp_search(text: &[u8], pattern: &[u8]) -> ComputeResult<Vec<usize>> {
    if pattern.is_empty() {
        return Err(ComputeError::InvalidInput("Pattern cannot be empty".to_string()));
    }
    
    if text.is_empty() {
        return Ok(Vec::new());
    }
    
    // Compute the failure function (partial match table)
    let failure_table = compute_kmp_failure_table(pattern);
    
    // Perform the search
    let mut matches = Vec::new();
    let mut j = 0; // position in pattern
    
    for (i, &c) in text.iter().enumerate() {
        // Try to match current character
        while j > 0 && pattern[j] != c {
            j = failure_table[j - 1];
        }
        
        // If characters match, advance pattern position
        if pattern[j] == c {
            j += 1;
        }
        
        // If we reached the end of the pattern, we found a match
        if j == pattern.len() {
            matches.push(i + 1 - j);
            j = failure_table[j - 1];
        }
    }
    
    Ok(matches)
}

/// Compute the failure function table for KMP algorithm
fn compute_kmp_failure_table(pattern: &[u8]) -> Vec<usize> {
    let m = pattern.len();
    let mut failure = vec![0; m];
    let mut j = 0;
    
    for i in 1..m {
        // If there's a mismatch, backtrack using previous entries
        while j > 0 && pattern[j] != pattern[i] {
            j = failure[j - 1];
        }
        
        // If there's a match, advance to next position in pattern
        if pattern[j] == pattern[i] {
            j += 1;
        }
        
        // Set the failure function value for current position
        failure[i] = j;
    }
    
    failure
}

/// Perform a substring search with the Boyer-Moore algorithm
///
/// Boyer-Moore is often faster than KMP for larger alphabets like DNA/protein sequences,
/// particularly when the pattern is relatively long.
pub fn boyer_moore_search(text: &[u8], pattern: &[u8]) -> ComputeResult<Vec<usize>> {
    if pattern.is_empty() {
        return Err(ComputeError::InvalidInput("Pattern cannot be empty".to_string()));
    }
    
    if text.is_empty() {
        return Ok(Vec::new());
    }
    
    // Compute the bad character heuristic
    let bad_char = compute_bad_char_table(pattern);
    
    // Compute the good suffix heuristic
    let good_suffix = compute_good_suffix_table(pattern);
    
    // Perform the search
    let mut matches = Vec::new();
    let m = pattern.len();
    let n = text.len();
    
    let mut i = 0;
    while i <= n - m {
        let mut j = m - 1;
        
        // Compare pattern and text from right to left
        while j < m && pattern[j] == text[i + j] {
            if j == 0 {
                break;
            }
            j -= 1;
        }
        
        // If we matched the entire pattern
        if j == 0 && pattern[0] == text[i] {
            matches.push(i);
            i += good_suffix[0];
        } else {
            // Otherwise, shift based on the maximum of bad character and good suffix heuristics
            let bc_shift = bad_char[text[i + j] as usize].max(1);
            let gs_shift = good_suffix[j];
            i += cmp::max(bc_shift, gs_shift);
        }
    }
    
    Ok(matches)
}

/// Compute the bad character table for Boyer-Moore algorithm
fn compute_bad_char_table(pattern: &[u8]) -> Vec<usize> {
    let m = pattern.len();
    let mut bad_char = vec![m; 256]; // default shift is pattern length
    
    // Populate bad character table
    for (i, &c) in pattern.iter().enumerate().take(m - 1) {
        bad_char[c as usize] = m - 1 - i;
    }
    
    bad_char
}

/// Compute the good suffix table for Boyer-Moore algorithm
fn compute_good_suffix_table(pattern: &[u8]) -> Vec<usize> {
    let m = pattern.len();
    let mut good_suffix = vec![0; m];
    let mut suffix_table = compute_suffix_table(pattern);
    
    // Initialize with default value
    for i in 0..m {
        good_suffix[i] = m;
    }
    
    // Case 1: pattern substring matches a suffix of pattern
    let mut j = 0;
    for i in (0..m-1).rev() {
        if suffix_table[i] == i + 1 {
            while j < m - 1 - i {
                if good_suffix[j] == m {
                    good_suffix[j] = m - 1 - i;
                }
                j += 1;
            }
        }
    }
    
    // Case 2: suffix of pattern occurs as prefix of pattern
    for i in 0..m - 1 {
        good_suffix[m - 1 - suffix_table[i]] = m - 1 - i;
    }
    
    good_suffix
}

/// Compute the suffix table for Boyer-Moore algorithm
fn compute_suffix_table(pattern: &[u8]) -> Vec<usize> {
    let m = pattern.len();
    let mut suffix = vec![0; m];
    
    suffix[m - 1] = m;
    let mut g = m - 1;
    
    for i in (0..m - 1).rev() {
        if i > g && suffix[m - 1 - (m - 1 - i)] < i - g {
            suffix[i] = suffix[m - 1 - (m - 1 - i)];
        } else {
            if i < g {
                g = i;
            }
            let mut j = 0;
            while g >= 0 && pattern[g] == pattern[m - 1 - j] {
                g -= 1;
                j += 1;
            }
            suffix[i] = j;
        }
    }
    
    suffix
}

/// Reverse a sequence in-place
pub fn reverse_in_place(sequence: &mut [u8]) {
    let len = sequence.len();
    if len <= 1 {
        return;
    }
    
    for i in 0..len / 2 {
        sequence.swap(i, len - 1 - i);
    }
}

/// Reverse a sequence, returning a new vector
pub fn reverse(sequence: &[u8]) -> Vec<u8> {
    let mut result = sequence.to_vec();
    reverse_in_place(&mut result);
    result
}

/// Complement a DNA sequence in-place
pub fn complement_dna_in_place(sequence: &mut [u8]) {
    for base in sequence.iter_mut() {
        *base = match *base {
            b'A' | b'a' => b'T',
            b'C' | b'c' => b'G',
            b'G' | b'g' => b'C',
            b'T' | b't' => b'A',
            b'U' | b'u' => b'A', // Handle RNA as well
            b'N' | b'n' => b'N',
            _ => *base, // Keep other characters unchanged
        };
    }
}

/// Complement a DNA sequence, returning a new vector
pub fn complement_dna(sequence: &[u8]) -> Vec<u8> {
    let mut result = sequence.to_vec();
    complement_dna_in_place(&mut result);
    result
}

/// Reverse-complement a DNA sequence in-place
pub fn reverse_complement_dna_in_place(sequence: &mut [u8]) {
    complement_dna_in_place(sequence);
    reverse_in_place(sequence);
}

/// Reverse-complement a DNA sequence, returning a new vector
pub fn reverse_complement_dna(sequence: &[u8]) -> Vec<u8> {
    let mut result = sequence.to_vec();
    reverse_complement_dna_in_place(&mut result);
    result
}

/// Count occurrences of each base in a DNA sequence
pub fn count_bases(sequence: &[u8]) -> [usize; 5] {
    let mut counts = [0, 0, 0, 0, 0]; // A, C, G, T, N/Other
    
    // Use SIMD-accelerated counting if available
    if simd::has_avx2() || simd::has_sse41() {
        counts[0] = simd::count_byte(sequence, b'A') + simd::count_byte(sequence, b'a');
        counts[1] = simd::count_byte(sequence, b'C') + simd::count_byte(sequence, b'c');
        counts[2] = simd::count_byte(sequence, b'G') + simd::count_byte(sequence, b'g');
        counts[3] = simd::count_byte(sequence, b'T') + simd::count_byte(sequence, b't') +
                    simd::count_byte(sequence, b'U') + simd::count_byte(sequence, b'u');
        
        // Count Ns and others
        let total_bases = counts[0] + counts[1] + counts[2] + counts[3];
        counts[4] = sequence.len() - total_bases;
    } else {
        // Fallback to scalar counting
        for &base in sequence {
            match base {
                b'A' | b'a' => counts[0] += 1,
                b'C' | b'c' => counts[1] += 1,
                b'G' | b'g' => counts[2] += 1,
                b'T' | b't' | b'U' | b'u' => counts[3] += 1,
                _ => counts[4] += 1,
            }
        }
    }
    
    counts
}

/// Calculate GC content of a DNA/RNA sequence
pub fn gc_content(sequence: &[u8]) -> f64 {
    if sequence.is_empty() {
        return 0.0;
    }
    
    let counts = count_bases(sequence);
    let gc_count = counts[1] + counts[2]; // C + G
    let total_bases = sequence.len() - counts[4]; // Exclude N/other
    
    if total_bases == 0 {
        return 0.0;
    }
    
    (gc_count as f64) / (total_bases as f64) * 100.0
}

/// Transcribe DNA to RNA (T -> U)
pub fn transcribe(dna: &[u8]) -> Vec<u8> {
    let mut rna = Vec::with_capacity(dna.len());
    
    for &base in dna {
        match base {
            b'T' => rna.push(b'U'),
            b't' => rna.push(b'u'),
            _ => rna.push(base),
        }
    }
    
    rna
}

/// Reverse-transcribe RNA to DNA (U -> T)
pub fn reverse_transcribe(rna: &[u8]) -> Vec<u8> {
    let mut dna = Vec::with_capacity(rna.len());
    
    for &base in rna {
        match base {
            b'U' => dna.push(b'T'),
            b'u' => dna.push(b't'),
            _ => dna.push(base),
        }
    }
    
    dna
}

/// Generate random DNA sequence of given length
pub fn random_dna(length: usize) -> Vec<u8> {
    use rand::prelude::*;
    let bases = [b'A', b'C', b'G', b'T'];
    let mut rng = rand::thread_rng();
    
    (0..length)
        .map(|_| *bases.choose(&mut rng).unwrap())
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    
    #[test]
    fn test_kmp_search() {
        // Test simple patterns
        let text = b"ACGTACGTACGTACGT";
        let pattern = b"ACGT";
        let matches = kmp_search(text, pattern).unwrap();
        assert_eq!(matches, vec![0, 4, 8, 12]);
        
        // Test pattern not found
        let pattern = b"AAAA";
        let matches = kmp_search(text, pattern).unwrap();
        assert_eq!(matches, Vec::<usize>::new());
        
        // Test empty text
        let text = b"";
        let pattern = b"A";
        let matches = kmp_search(text, pattern).unwrap();
        assert_eq!(matches, Vec::<usize>::new());
        
        // Test empty pattern
        let text = b"ACGT";
        let pattern = b"";
        let result = kmp_search(text, pattern);
        assert!(result.is_err());
    }
    
    #[test]
    fn test_boyer_moore_search() {
        // Test simple patterns
        let text = b"ACGTACGTACGTACGT";
        let pattern = b"ACGT";
        let matches = boyer_moore_search(text, pattern).unwrap();
        assert_eq!(matches, vec![0, 4, 8, 12]);
        
        // Test pattern not found
        let pattern = b"AAAA";
        let matches = boyer_moore_search(text, pattern).unwrap();
        assert_eq!(matches, Vec::<usize>::new());
        
        // Test empty text
        let text = b"";
        let pattern = b"A";
        let matches = boyer_moore_search(text, pattern).unwrap();
        assert_eq!(matches, Vec::<usize>::new());
        
        // Test empty pattern
        let text = b"ACGT";
        let pattern = b"";
        let result = boyer_moore_search(text, pattern);
        assert!(result.is_err());
    }
    
    #[test]
    fn test_reverse() {
        let seq = b"ACGT";
        let reversed = reverse(seq);
        assert_eq!(reversed, b"TGCA");
        
        // Test in-place reversal
        let mut seq_mut = b"ACGT".to_vec();
        reverse_in_place(&mut seq_mut);
        assert_eq!(seq_mut, b"TGCA");
    }
    
    #[test]
    fn test_complement_dna() {
        let seq = b"ACGT";
        let complemented = complement_dna(seq);
        assert_eq!(complemented, b"TGCA");
        
        // Test in-place complement
        let mut seq_mut = b"ACGT".to_vec();
        complement_dna_in_place(&mut seq_mut);
        assert_eq!(seq_mut, b"TGCA");
    }
    
    #[test]
    fn test_reverse_complement_dna() {
        let seq = b"ACGT";
        let rev_comp = reverse_complement_dna(seq);
        assert_eq!(rev_comp, b"ACGT"); // ACGT is its own reverse complement
        
        let seq = b"AACGTT";
        let rev_comp = reverse_complement_dna(seq);
        assert_eq!(rev_comp, b"AACGTT"); // AACGTT is its own reverse complement
        
        // Test in-place reverse complement
        let mut seq_mut = b"ACGT".to_vec();
        reverse_complement_dna_in_place(&mut seq_mut);
        assert_eq!(seq_mut, b"ACGT");
    }
    
    #[test]
    fn test_count_bases() {
        let seq = b"ACGTACGTNNACGT";
        let counts = count_bases(seq);
        assert_eq!(counts, [4, 3, 3, 4, 2]); // A, C, G, T, N/Other
        
        // Test with lower case
        let seq = b"acgtACGTnnACGT";
        let counts = count_bases(seq);
        assert_eq!(counts, [4, 3, 3, 4, 2]); // A, C, G, T, N/Other
    }
    
    #[test]
    fn test_gc_content() {
        // Test basic GC content
        let seq = b"ACGT";
        let gc = gc_content(seq);
        assert_eq!(gc, 50.0); // 2/4 = 50%
        
        // Test with Ns
        let seq = b"ACGTNNNNNNNN";
        let gc = gc_content(seq);
        assert_eq!(gc, 50.0); // Ns are excluded, so still 2/4 = 50%
        
        // Test empty sequence
        let seq = b"";
        let gc = gc_content(seq);
        assert_eq!(gc, 0.0);
        
        // Test 100% GC
        let seq = b"GCGC";
        let gc = gc_content(seq);
        assert_eq!(gc, 100.0);
    }
    
    #[test]
    fn test_transcription() {
        // Test DNA to RNA
        let dna = b"ACGT";
        let rna = transcribe(dna);
        assert_eq!(rna, b"ACGU");
        
        // Test with lower case
        let dna = b"ACGt";
        let rna = transcribe(dna);
        assert_eq!(rna, b"ACGu");
        
        // Test RNA to DNA
        let rna = b"ACGU";
        let dna = reverse_transcribe(rna);
        assert_eq!(dna, b"ACGT");
    }
    
    #[test]
    fn test_random_dna() {
        // Test random sequence generation
        let dna = random_dna(100);
        assert_eq!(dna.len(), 100);
        
        // All bases should be valid
        for &base in &dna {
            assert!(base == b'A' || base == b'C' || base == b'G' || base == b'T');
        }
    }
}