//! Core system optimization primitives
//!
//! This module provides low-level optimizations for memory, 
//! parallelism, IO, and SIMD operations.

pub mod parallel;
pub mod memory;
pub mod io;
pub mod simd;

/// Version feature detection for runtime optimization
pub fn detect_cpu_features() -> CpuFeatures {
    CpuFeatures {
        has_avx2: is_x86_feature_detected!("avx2"),
        has_avx512: is_x86_feature_detected!("avx512f"),
        has_sse41: is_x86_feature_detected!("sse4.1"),
        has_sse42: is_x86_feature_detected!("sse4.2"),
    }
}

/// CPU feature detection results
#[derive(Debug, Clone, Copy)]
pub struct CpuFeatures {
    /// Whether AVX2 instructions are available
    pub has_avx2: bool,
    /// Whether AVX-512 instructions are available
    pub has_avx512: bool,
    /// Whether SSE4.1 instructions are available
    pub has_sse41: bool,
    /// Whether SSE4.2 instructions are available
    pub has_sse42: bool,
}

/// Initialize the core engine with optimal settings for the current system
pub fn initialize() {
    // Detect CPU features
    let features = detect_cpu_features();
    
    // Initialize memory subsystem
    memory::initialize();
    
    // Initialize parallel execution
    parallel::initialize_thread_pool();
    
    // Set up optimal I/O configuration
    io::initialize();
    
    // Initialize SIMD settings based on detected features
    simd::initialize(features);
    
    log::info!("Core engine initialized with features: {:?}", features);
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_cpu_feature_detection() {
        let features = detect_cpu_features();
        // At least one of these should be true on modern CPUs
        assert!(features.has_sse41 || features.has_sse42 || features.has_avx2 || features.has_avx512);
    }
}