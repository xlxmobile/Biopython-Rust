//! I/O module
//!
//! This module provides I/O operations for reading and writing
//! biological sequence files.

pub mod fasta;

use crate::engines;

/// Initialize the I/O module
pub fn initialize() {
    // Initialize any required state for the I/O module
}

/// Convenience re-exports
pub use fasta::{read_fasta, write_fasta, FastaRecord};

#[cfg(test)]
mod tests {
    use super::*;
    
    #[test]
    fn test_module_initialization() {
        initialize();
    }
}