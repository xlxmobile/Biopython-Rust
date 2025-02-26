//! High-performance computation engines for bioinformatics operations
//!
//! This module contains low-level, optimized implementations of common
//! operations used in bioinformatics.

pub mod core;
pub mod compute;
pub mod storage;

/// Execution configuration for parallel operations
#[derive(Debug, Clone)]
pub struct ExecutionConfig {
    /// Number of threads to use for parallel execution
    pub num_threads: usize,
    /// Size of chunks for parallel processing
    pub chunk_size: usize,
    /// Whether to use SIMD operations when available
    pub use_simd: bool,
}

impl Default for ExecutionConfig {
    fn default() -> Self {
        Self {
            num_threads: core::parallel::default_num_threads(),
            chunk_size: 1024 * 1024, // 1MB default chunk size
            use_simd: true,
        }
    }
}

/// Error types for engine operations
#[derive(Debug, thiserror::Error)]
pub enum EngineError {
    #[error("I/O error: {0}")]
    IoError(#[from] std::io::Error),
    
    #[error("Memory allocation error: {0}")]
    MemoryError(String),
    
    #[error("Invalid sequence data: {0}")]
    InvalidSequenceData(String),
    
    #[error("Parallel execution error: {0}")]
    ParallelExecutionError(String),
    
    #[error("Unsupported operation: {0}")]
    UnsupportedOperation(String),
}

/// Result type for engine operations
pub type EngineResult<T> = Result<T, EngineError>;

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_default_execution_config() {
        let config = ExecutionConfig::default();
        assert!(config.num_threads > 0);
        assert_eq!(config.chunk_size, 1024 * 1024);
        assert!(config.use_simd);
    }
}