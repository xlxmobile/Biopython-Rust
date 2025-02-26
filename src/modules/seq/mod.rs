//! Sequence module
//!
//! This module provides sequence types and operations for bioinformatics.

pub mod sequence;
pub mod alphabet;

use crate::engines;

/// Initialize the sequence module
pub fn initialize() {
    // Initialize any required state for the sequence module
}

/// Convenience re-exports
pub use sequence::{Sequence, SequenceView, SequenceError};
pub use alphabet::{Alphabet, DNAAlphabet, RNAAlphabet, ProteinAlphabet};

#[cfg(test)]
mod tests {
    use super::*;
    
    #[test]
    fn test_module_initialization() {
        initialize();
    }
}