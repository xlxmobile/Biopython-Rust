//! Bioinformatics modules
//!
//! This module provides high-level bioinformatics functionality built on top of
//! the optimized engines.

pub mod seq;
pub mod io;

use crate::engines;

/// Initialize the modules
pub fn initialize() {
    // Initialize the underlying engines
    engines::core::initialize();
    
    // Initialize sequence-specific modules
    seq::initialize();
    
    // Initialize I/O-specific modules
    io::initialize();
}

#[cfg(test)]
mod tests {
    use super::*;
    
    #[test]
    fn test_initialization() {
        // Simple sanity check
        initialize();
    }
}