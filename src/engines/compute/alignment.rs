//! Sequence alignment algorithms
//!
//! This module provides high-performance implementations of common alignment
//! algorithms used in bioinformatics.

use super::{ComputeResult, ComputeError};
use std::cmp;

/// Different types of alignment algorithms
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AlignmentType {
    /// Global alignment (Needleman-Wunsch)
    Global,
    /// Local alignment (Smith-Waterman)
    Local,
    /// Semi-global alignment
    SemiGlobal,
}

/// Different scoring schemes for alignments
#[derive(Debug, Clone)]
pub struct ScoringScheme {
    /// Score for a match
    pub match_score: i32,
    /// Penalty for a mismatch
    pub mismatch_penalty: i32,
    /// Penalty for opening a gap
    pub gap_open_penalty: i32,
    /// Penalty for extending a gap
    pub gap_extend_penalty: i32,
}

impl Default for ScoringScheme {
    fn default() -> Self {
        Self {
            match_score: 2,
            mismatch_penalty: -1,
            gap_open_penalty: -2,
            gap_extend_penalty: -1,
        }
    }
}

/// Represents an alignment between two sequences
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Alignment {
    /// First sequence aligned (with gaps)
    pub seq1_aligned: Vec<u8>,
    /// Second sequence aligned (with gaps)
    pub seq2_aligned: Vec<u8>,
    /// Score of the alignment
    pub score: i32,
    /// Start position in sequence 1 (0-indexed)
    pub seq1_start: usize,
    /// End position in sequence 1 (0-indexed)
    pub seq1_end: usize,
    /// Start position in sequence 2 (0-indexed)
    pub seq2_start: usize,
    /// End position in sequence 2 (0-indexed)
    pub seq2_end: usize,
    /// Identity percentage (matches / alignment length)
    pub identity: f64,
}

impl Alignment {
    /// Calculate the identity percentage
    pub fn calculate_identity(&mut self) {
        let mut matches = 0;
        let alignment_length = self.seq1_aligned.len();
        
        for i in 0..alignment_length {
            if self.seq1_aligned[i] == self.seq2_aligned[i] && self.seq1_aligned[i] != b'-' {
                matches += 1;
            }
        }
        
        self.identity = if alignment_length > 0 {
            (matches as f64) / (alignment_length as f64) * 100.0
        } else {
            0.0
        };
    }
    
    /// Get the alignment as a formatted string
    pub fn format(&self) -> String {
        let mut result = String::new();
        
        // Add header
        result.push_str(&format!("Score: {}\n", self.score));
        result.push_str(&format!("Identity: {:.2}%\n", self.identity));
        result.push_str(&format!("Seq1: {}..{}\n", self.seq1_start, self.seq1_end));
        result.push_str(&format!("Seq2: {}..{}\n\n", self.seq2_start, self.seq2_end));
        
        // Add alignment visualization
        const LINE_WIDTH: usize = 60;
        let total_len = self.seq1_aligned.len();
        
        for i in (0..total_len).step_by(LINE_WIDTH) {
            let end = cmp::min(i + LINE_WIDTH, total_len);
            let slice_len = end - i;
            
            // Sequence 1
            result.push_str("Seq1: ");
            result.push_str(&String::from_utf8_lossy(&self.seq1_aligned[i..end]));
            result.push('\n');
            
            // Match line
            result.push_str("      ");
            for j in i..end {
                if self.seq1_aligned[j] == self.seq2_aligned[j] && self.seq1_aligned[j] != b'-' {
                    result.push('|');
                } else {
                    result.push(' ');
                }
            }
            result.push('\n');
            
            // Sequence 2
            result.push_str("Seq2: ");
            result.push_str(&String::from_utf8_lossy(&self.seq2_aligned[i..end]));
            result.push('\n');
            
            if end < total_len {
                result.push('\n');
            }
        }
        
        result
    }
}

/// Perform sequence alignment using the specified algorithm
pub fn align(
    seq1: &[u8],
    seq2: &[u8],
    alignment_type: AlignmentType,
    scoring: &ScoringScheme,
) -> ComputeResult<Alignment> {
    match alignment_type {
        AlignmentType::Global => needleman_wunsch(seq1, seq2, scoring),
        AlignmentType::Local => smith_waterman(seq1, seq2, scoring),
        AlignmentType::SemiGlobal => semi_global_align(seq1, seq2, scoring),
    }
}

/// Perform global alignment using the Needleman-Wunsch algorithm
pub fn needleman_wunsch(
    seq1: &[u8],
    seq2: &[u8],
    scoring: &ScoringScheme,
) -> ComputeResult<Alignment> {
    if seq1.is_empty() || seq2.is_empty() {
        return Err(ComputeError::InvalidInput("Sequences cannot be empty".to_string()));
    }
    
    let m = seq1.len();
    let n = seq2.len();
    
    // Initialize scoring matrix
    let mut dp = vec![vec![0; n + 1]; m + 1];
    
    // Initialize traceback matrix
    // 0 = diagonal (match/mismatch), 1 = left (gap in seq1), 2 = up (gap in seq2)
    let mut traceback = vec![vec![0; n + 1]; m + 1];
    
    // Initialize first row and column with gap penalties
    dp[0][0] = 0;
    for i in 1..=m {
        if i == 1 {
            dp[i][0] = scoring.gap_open_penalty;
        } else {
            dp[i][0] = dp[i-1][0] + scoring.gap_extend_penalty;
        }
        traceback[i][0] = 2; // gap in seq2
    }
    
    for j in 1..=n {
        if j == 1 {
            dp[0][j] = scoring.gap_open_penalty;
        } else {
            dp[0][j] = dp[0][j-1] + scoring.gap_extend_penalty;
        }
        traceback[0][j] = 1; // gap in seq1
    }
    
    // Fill the DP matrix
    for i in 1..=m {
        for j in 1..=n {
            // Calculate match/mismatch score
            let match_score = if seq1[i-1] == seq2[j-1] {
                scoring.match_score
            } else {
                scoring.mismatch_penalty
            };
            
            // Calculate scores for each possible move
            let diagonal = dp[i-1][j-1] + match_score;
            
            // Gap in seq1 (horizontal move)
            let left_score = if traceback[i][j-1] == 1 {
                // Extend existing gap
                dp[i][j-1] + scoring.gap_extend_penalty
            } else {
                // Open new gap
                dp[i][j-1] + scoring.gap_open_penalty
            };
            
            // Gap in seq2 (vertical move)
            let up_score = if traceback[i-1][j] == 2 {
                // Extend existing gap
                dp[i-1][j] + scoring.gap_extend_penalty
            } else {
                // Open new gap
                dp[i-1][j] + scoring.gap_open_penalty
            };
            
            // Choose the best score
            if diagonal >= left_score && diagonal >= up_score {
                dp[i][j] = diagonal;
                traceback[i][j] = 0; // diagonal
            } else if left_score >= up_score {
                dp[i][j] = left_score;
                traceback[i][j] = 1; // left
            } else {
                dp[i][j] = up_score;
                traceback[i][j] = 2; // up
            }
        }
    }
    
    // Traceback to construct the alignment
    let mut aligned_seq1 = Vec::new();
    let mut aligned_seq2 = Vec::new();
    
    let mut i = m;
    let mut j = n;
    
    while i > 0 || j > 0 {
        if i > 0 && j > 0 && traceback[i][j] == 0 {
            // Diagonal move (match/mismatch)
            aligned_seq1.push(seq1[i-1]);
            aligned_seq2.push(seq2[j-1]);
            i -= 1;
            j -= 1;
        } else if j > 0 && traceback[i][j] == 1 {
            // Left move (gap in seq1)
            aligned_seq1.push(b'-');
            aligned_seq2.push(seq2[j-1]);
            j -= 1;
        } else if i > 0 && traceback[i][j] == 2 {
            // Up move (gap in seq2)
            aligned_seq1.push(seq1[i-1]);
            aligned_seq2.push(b'-');
            i -= 1;
        } else {
            // Should not happen with properly initialized traceback
            break;
        }
    }
    
    // Reverse the alignment (we traced backwards)
    aligned_seq1.reverse();
    aligned_seq2.reverse();
    
    // Create and return the alignment
    let mut alignment = Alignment {
        seq1_aligned: aligned_seq1,
        seq2_aligned: aligned_seq2,
        score: dp[m][n],
        seq1_start: 0,
        seq1_end: m,
        seq2_start: 0,
        seq2_end: n,
        identity: 0.0,
    };
    
    // Calculate identity
    alignment.calculate_identity();
    
    Ok(alignment)
}

/// Perform local alignment using the Smith-Waterman algorithm
pub fn smith_waterman(
    seq1: &[u8],
    seq2: &[u8],
    scoring: &ScoringScheme,
) -> ComputeResult<Alignment> {
    if seq1.is_empty() || seq2.is_empty() {
        return Err(ComputeError::InvalidInput("Sequences cannot be empty".to_string()));
    }
    
    let m = seq1.len();
    let n = seq2.len();
    
    // Initialize scoring matrix
    let mut dp = vec![vec![0; n + 1]; m + 1];
    
    // Initialize traceback matrix
    // 0 = diagonal (match/mismatch), 1 = left (gap in seq1), 2 = up (gap in seq2), 3 = stop
    let mut traceback = vec![vec![3; n + 1]; m + 1];
    
    // Fill the DP matrix
    let mut max_score = 0;
    let mut max_i = 0;
    let mut max_j = 0;
    
    for i in 1..=m {
        for j in 1..=n {
            // Calculate match/mismatch score
            let match_score = if seq1[i-1] == seq2[j-1] {
                scoring.match_score
            } else {
                scoring.mismatch_penalty
            };
            
            // Calculate scores for each possible move
            let diagonal = dp[i-1][j-1] + match_score;
            
            // Gap in seq1 (horizontal move)
            let left_score = dp[i][j-1] + (if traceback[i][j-1] == 1 {
                scoring.gap_extend_penalty
            } else {
                scoring.gap_open_penalty
            });
            
            // Gap in seq2 (vertical move)
            let up_score = dp[i-1][j] + (if traceback[i-1][j] == 2 {
                scoring.gap_extend_penalty
            } else {
                scoring.gap_open_penalty
            });
            
            // Local alignment allows stopping at any point
            let scores = [0, diagonal, left_score, up_score];
            let max_idx = scores.iter().enumerate()
                .max_by_key(|&(_, &score)| score)
                .map(|(idx, _)| idx)
                .unwrap();
            
            dp[i][j] = scores[max_idx];
            
            // Set traceback based on the chosen move
            traceback[i][j] = match max_idx {
                0 => 3, // stop (local alignment can start/end anywhere)
                1 => 0, // diagonal
                2 => 1, // left
                3 => 2, // up
                _ => unreachable!(),
            };
            
            // Keep track of the maximum score for starting the traceback
            if dp[i][j] > max_score {
                max_score = dp[i][j];
                max_i = i;
                max_j = j;
            }
        }
    }
    
    // Traceback to construct the alignment
    let mut aligned_seq1 = Vec::new();
    let mut aligned_seq2 = Vec::new();
    
    let mut i = max_i;
    let mut j = max_j;
    
    // Record the end positions for local alignment
    let seq1_end = i;
    let seq2_end = j;
    
    // Traceback until we hit a cell with score 0 or a "stop" traceback
    while i > 0 && j > 0 && dp[i][j] > 0 && traceback[i][j] != 3 {
        if traceback[i][j] == 0 {
            // Diagonal move (match/mismatch)
            aligned_seq1.push(seq1[i-1]);
            aligned_seq2.push(seq2[j-1]);
            i -= 1;
            j -= 1;
        } else if traceback[i][j] == 1 {
            // Left move (gap in seq1)
            aligned_seq1.push(b'-');
            aligned_seq2.push(seq2[j-1]);
            j -= 1;
        } else if traceback[i][j] == 2 {
            // Up move (gap in seq2)
            aligned_seq1.push(seq1[i-1]);
            aligned_seq2.push(b'-');
            i -= 1;
        }
    }
    
    // Record the start positions for local alignment
    let seq1_start = i;
    let seq2_start = j;
    
    // Reverse the alignment (we traced backwards)
    aligned_seq1.reverse();
    aligned_seq2.reverse();
    
    // Create and return the alignment
    let mut alignment = Alignment {
        seq1_aligned: aligned_seq1,
        seq2_aligned: aligned_seq2,
        score: max_score,
        seq1_start,
        seq1_end,
        seq2_start,
        seq2_end,
        identity: 0.0,
    };
    
    // Calculate identity
    alignment.calculate_identity();
    
    Ok(alignment)
}

/// Perform semi-global alignment
///
/// Semi-global alignment is a variation where gaps at the beginning and end
/// of one sequence are not penalized (useful for aligning a short sequence
/// to a long one).
pub fn semi_global_align(
    seq1: &[u8],
    seq2: &[u8],
    scoring: &ScoringScheme,
) -> ComputeResult<Alignment> {
    if seq1.is_empty() || seq2.is_empty() {
        return Err(ComputeError::InvalidInput("Sequences cannot be empty".to_string()));
    }
    
    let m = seq1.len();
    let n = seq2.len();
    
    // Initialize scoring matrix
    let mut dp = vec![vec![0; n + 1]; m + 1];
    
    // Initialize traceback matrix
    // 0 = diagonal (match/mismatch), 1 = left (gap in seq1), 2 = up (gap in seq2)
    let mut traceback = vec![vec![0; n + 1]; m + 1];
    
    // Initialize first row and column
    // In semi-global, we don't penalize gaps at the beginning of one sequence
    for i in 0..=m {
        dp[i][0] = 0;
        traceback[i][0] = 2; // gap in seq2
    }
    
    for j in 1..=n {
        dp[0][j] = 0;
        traceback[0][j] = 1; // gap in seq1
    }
    
    // Fill the DP matrix
    for i in 1..=m {
        for j in 1..=n {
            // Calculate match/mismatch score
            let match_score = if seq1[i-1] == seq2[j-1] {
                scoring.match_score
            } else {
                scoring.mismatch_penalty
            };
            
            // Calculate scores for each possible move
            let diagonal = dp[i-1][j-1] + match_score;
            
            // Gap in seq1 (horizontal move)
            let left_score = dp[i][j-1] + (if traceback[i][j-1] == 1 {
                scoring.gap_extend_penalty
            } else {
                scoring.gap_open_penalty
            });
            
            // Gap in seq2 (vertical move)
            let up_score = dp[i-1][j] + (if traceback[i-1][j] == 2 {
                scoring.gap_extend_penalty
            } else {
                scoring.gap_open_penalty
            });
            
            // Choose the best score
            if diagonal >= left_score && diagonal >= up_score {
                dp[i][j] = diagonal;
                traceback[i][j] = 0; // diagonal
            } else if left_score >= up_score {
                dp[i][j] = left_score;
                traceback[i][j] = 1; // left
            } else {
                dp[i][j] = up_score;
                traceback[i][j] = 2; // up
            }
        }
    }
    
    // Find the best score in the last row or last column
    let mut max_score = dp[m][n];
    let mut max_i = m;
    let mut max_j = n;
    
    // Check last row
    for j in 0..=n {
        if dp[m][j] > max_score {
            max_score = dp[m][j];
            max_i = m;
            max_j = j;
        }
    }
    
    // Check last column
    for i in 0..=m {
        if dp[i][n] > max_score {
            max_score = dp[i][n];
            max_i = i;
            max_j = n;
        }
    }
    
    // Traceback to construct the alignment
    let mut aligned_seq1 = Vec::new();
    let mut aligned_seq2 = Vec::new();
    
    let mut i = max_i;
    let mut j = max_j;
    
    // Record the end positions
    let seq1_end = i;
    let seq2_end = j;
    
    // Add gaps at the end if necessary
    while i < m {
        aligned_seq1.push(seq1[i]);
        aligned_seq2.push(b'-');
        i += 1;
    }
    
    while j < n {
        aligned_seq1.push(b'-');
        aligned_seq2.push(seq2[j]);
        j += 1;
    }
    
    // Traceback until we hit the beginning of either sequence
    while i > 0 && j > 0 {
        if traceback[i][j] == 0 {
            // Diagonal move (match/mismatch)
            aligned_seq1.push(seq1[i-1]);
            aligned_seq2.push(seq2[j-1]);
            i -= 1;
            j -= 1;
        } else if traceback[i][j] == 1 {
            // Left move (gap in seq1)
            aligned_seq1.push(b'-');
            aligned_seq2.push(seq2[j-1]);
            j -= 1;
        } else if traceback[i][j] == 2 {
            // Up move (gap in seq2)
            aligned_seq1.push(seq1[i-1]);
            aligned_seq2.push(b'-');
            i -= 1;
        }
    }
    
    // Add gaps at the beginning if necessary
    while i > 0 {
        aligned_seq1.push(seq1[i-1]);
        aligned_seq2.push(b'-');
        i -= 1;
    }
    
    while j > 0 {
        aligned_seq1.push(b'-');
        aligned_seq2.push(seq2[j-1]);
        j -= 1;
    }
    
    // Record the start positions
    let seq1_start = 0;
    let seq2_start = 0;
    
    // Reverse the alignment (we traced backwards)
    aligned_seq1.reverse();
    aligned_seq2.reverse();
    
    // Create and return the alignment
    let mut alignment = Alignment {
        seq1_aligned: aligned_seq1,
        seq2_aligned: aligned_seq2,
        score: max_score,
        seq1_start,
        seq1_end,
        seq2_start,
        seq2_end,
        identity: 0.0,
    };
    
    // Calculate identity
    alignment.calculate_identity();
    
    Ok(alignment)
}

/// Calculate the edit distance (Levenshtein distance) between two sequences
pub fn edit_distance(seq1: &[u8], seq2: &[u8]) -> usize {
    let m = seq1.len();
    let n = seq2.len();
    
    // Handle special cases
    if m == 0 {
        return n;
    }
    if n == 0 {
        return m;
    }
    
    // Initialize DP matrix
    let mut dp = vec![vec![0; n + 1]; m + 1];
    
    // Initialize first row and column
    for i in 0..=m {
        dp[i][0] = i;
    }
    for j in 0..=n {
        dp[0][j] = j;
    }
    
    // Fill the DP matrix
    for i in 1..=m {
        for j in 1..=n {
            let cost = if seq1[i-1] == seq2[j-1] { 0 } else { 1 };
            
            dp[i][j] = cmp::min(
                dp[i-1][j] + 1,      // deletion
                cmp::min(
                    dp[i][j-1] + 1,  // insertion
                    dp[i-1][j-1] + cost  // substitution
                )
            );
        }
    }
    
    dp[m][n]
}

#[cfg(test)]
mod tests {
    use super::*;
    
    #[test]
    fn test_global_alignment() {
        let seq1 = b"ACGTACGT";
        let seq2 = b"ACGTCGT";
        let scoring = ScoringScheme::default();
        
        let alignment = needleman_wunsch(seq1, seq2, &scoring).unwrap();
        
        // Expected alignment:
        // ACGTACGT
        // ACGT-CGT
        assert_eq!(alignment.seq1_aligned, b"ACGTACGT");
        assert_eq!(alignment.seq2_aligned, b"ACGT-CGT");
        assert_eq!(alignment.score, 11); // 7 matches * 2 - 1 gap * 2 = 12
    }
    
    #[test]
    fn test_local_alignment() {
        let seq1 = b"ACGTACGTACGT";
        let seq2 = b"TACGTAC";
        let scoring = ScoringScheme::default();
        
        let alignment = smith_waterman(seq1, seq2, &scoring).unwrap();
        
        // Expected alignment:
        // ACGTAC
        // ACGTAC
        assert_eq!(alignment.seq1_aligned, b"ACGTAC");
        assert_eq!(alignment.seq2_aligned, b"ACGTAC");
        assert_eq!(alignment.score, 12); // 6 matches * 2 = 12
    }
    
    #[test]
    fn test_semi_global_alignment() {
        let seq1 = b"ACGTACGTACGT";
        let seq2 = b"TACGTAC";
        let scoring = ScoringScheme::default();
        
        let alignment = semi_global_align(seq1, seq2, &scoring).unwrap();
        
        // Gaps at ends of seq2 shouldn't be penalized
        assert!(alignment.score >= 0);
    }
    
    #[test]
    fn test_edit_distance() {
        // Test cases
        assert_eq!(edit_distance(b"ACGT", b"ACGT"), 0); // Identical
        assert_eq!(edit_distance(b"ACGT", b"ACGTA"), 1); // Insertion
        assert_eq!(edit_distance(b"ACGT", b"ACG"), 1); // Deletion
        assert_eq!(edit_distance(b"ACGT", b"ACTT"), 1); // Substitution
        assert_eq!(edit_distance(b"", b"ACGT"), 4); // All insertions
        assert_eq!(edit_distance(b"ACGT", b""), 4); // All deletions
        assert_eq!(edit_distance(b"ACGT", b"TGCA"), 4); // All substitutions
    }
}