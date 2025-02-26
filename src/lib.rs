//! High-performance bioinformatics sequence operations library
//!
//! This library provides optimized implementations for common bioinformatics
//! sequence operations, with a focus on performance and memory efficiency.
//!
//! # Features
//! - Memory-efficient sequence storage with 2-bit and 4-bit encoding
//! - Memory-mapped files for processing large sequences
//! - Parallel processing framework for sequence operations
//! - SIMD-accelerated implementations for common operations

#![cfg_attr(not(feature = "std"), no_std)]

pub mod engines;
pub mod modules;

// Re-export commonly used items
pub use modules::seq::sequence::{Sequence, SequenceView};
pub use modules::seq::alphabet::{Alphabet, DNAAlphabet, RNAAlphabet, ProteinAlphabet};

/// Version information
pub const VERSION: &str = env!("CARGO_PKG_VERSION");

/// Library name
pub const PKG_NAME: &str = env!("CARGO_PKG_NAME");

/// Initialize the library with default settings
#[inline]
pub fn init() {
    // Initialize global state, loggers, etc. if needed
    engines::core::parallel::initialize_thread_pool();
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_library_initialization() {
        init();
        // Simple sanity check
        assert_eq!(PKG_NAME, "biopython-rust");
    }
}