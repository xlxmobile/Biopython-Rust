//! Compute primitives for bioinformatics operations
//!
//! This module provides high-performance implementations of common
//! bioinformatics computation primitives.

pub mod string_ops;
pub mod alignment;

use crate::engines::core::parallel::ParallelChunkProcessor;

/// Compute operation result type
pub type ComputeResult<T> = Result<T, ComputeError>;

/// Error types for compute operations
#[derive(Debug, thiserror::Error)]
pub enum ComputeError {
    #[error("Invalid input data: {0}")]
    InvalidInput(String),
    
    #[error("Computation error: {0}")]
    ComputationError(String),
    
    #[error("Operation not supported: {0}")]
    UnsupportedOperation(String),
    
    #[error("Resource limit exceeded: {0}")]
    ResourceLimitExceeded(String),
}

/// Trait for parallelizable compute operations
pub trait ParallelCompute<T, R> {
    /// Execute the operation in parallel
    fn execute_parallel(&self, data: &[T], chunk_size: Option<usize>) -> ComputeResult<Vec<R>>;
    
    /// Check if the operation can be executed in parallel
    fn supports_parallel(&self) -> bool;
}

/// Base implementation for parallel compute operations
impl<T, R, F> ParallelCompute<T, R> for F
where
    T: Clone + Send + Sync + 'static,
    R: Send + 'static,
    F: Fn(&T) -> ComputeResult<R> + Sync + Send + Clone + 'static,
{
    fn execute_parallel(&self, data: &[T], chunk_size: Option<usize>) -> ComputeResult<Vec<R>> {
        if !self.supports_parallel() {
            return Err(ComputeError::UnsupportedOperation(
                "This operation does not support parallel execution".to_string(),
            ));
        }
        
        // Create chunks for parallel processing
        let chunk_size = chunk_size.unwrap_or(1024);
        let chunks: Vec<Vec<T>> = data
            .chunks(chunk_size)
            .map(|chunk| chunk.to_vec())
            .collect();
        
        // Create parallel processor
        let processor = ParallelChunkProcessor::new(chunks);
        
        // Process chunks in parallel
        let f = self.clone();
        let results = processor.process(move |chunk| {
            chunk
                .iter()
                .map(|item| f(item))
                .collect::<Vec<ComputeResult<R>>>()
        });
        
        // Flatten results
        let mut output = Vec::with_capacity(data.len());
        for chunk_results in results {
            for result in chunk_results {
                match result {
                    Ok(value) => output.push(value),
                    Err(e) => return Err(e),
                }
            }
        }
        
        Ok(output)
    }
    
    fn supports_parallel(&self) -> bool {
        true
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    
    // Test function to multiply a number by 2
    fn double(x: &i32) -> ComputeResult<i32> {
        Ok(x * 2)
    }
    
    #[test]
    fn test_parallel_compute() {
        // Initialize the parallel execution framework
        crate::engines::core::parallel::initialize_thread_pool();
        
        // Create test data
        let data: Vec<i32> = (1..1001).collect();
        
        // Execute parallel computation
        let results = double.execute_parallel(&data, Some(100)).unwrap();
        
        // Check results
        assert_eq!(results.len(), data.len());
        for (i, &result) in results.iter().enumerate() {
            assert_eq!(result, data[i] * 2);
        }
    }
    
    #[test]
    fn test_error_propagation() {
        // Function that fails for certain inputs
        let fail_on_negative = |x: &i32| -> ComputeResult<i32> {
            if *x < 0 {
                Err(ComputeError::InvalidInput(format!("Negative input: {}", x)))
            } else {
                Ok(x * 2)
            }
        };
        
        // Create test data with a negative number
        let data: Vec<i32> = vec![1, 2, 3, -4, 5];
        
        // Execute parallel computation
        let result = fail_on_negative.execute_parallel(&data, None);
        
        // Should propagate the error
        assert!(result.is_err());
        match result {
            Err(ComputeError::InvalidInput(msg)) => {
                assert!(msg.contains("Negative input"));
            }
            _ => panic!("Expected InvalidInput error"),
        }
    }
}