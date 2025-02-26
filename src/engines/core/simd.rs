//! SIMD operations for sequence processing
//!
//! This module provides SIMD-accelerated implementations of common
//! sequence operations, with runtime feature detection and fallbacks.

use std::arch::x86_64::*;
use std::sync::atomic::{AtomicBool, Ordering};
use crate::engines::core::CpuFeatures;

// Track whether SIMD is available
static AVX2_AVAILABLE: AtomicBool = AtomicBool::new(false);
static SSE41_AVAILABLE: AtomicBool = AtomicBool::new(false);

/// Initialize SIMD settings based on detected CPU features
pub fn initialize(features: CpuFeatures) {
    AVX2_AVAILABLE.store(features.has_avx2, Ordering::SeqCst);
    SSE41_AVAILABLE.store(features.has_sse41, Ordering::SeqCst);
    
    log::info!("SIMD initialized - AVX2: {}, SSE4.1: {}", 
               features.has_avx2, features.has_sse41);
}

/// Check if AVX2 instructions are available
#[inline]
pub fn has_avx2() -> bool {
    AVX2_AVAILABLE.load(Ordering::Relaxed)
}

/// Check if SSE4.1 instructions are available
#[inline]
pub fn has_sse41() -> bool {
    SSE41_AVAILABLE.load(Ordering::Relaxed)
}

/// Count occurrences of a byte in a slice using the most efficient
/// available SIMD instruction set
pub fn count_byte(slice: &[u8], byte: u8) -> usize {
    if has_avx2() {
        unsafe { count_byte_avx2(slice, byte) }
    } else if has_sse41() {
        unsafe { count_byte_sse41(slice, byte) }
    } else {
        count_byte_scalar(slice, byte)
    }
}

/// Find the first occurrence of a byte in a slice using the most efficient
/// available SIMD instruction set
pub fn find_byte(slice: &[u8], byte: u8) -> Option<usize> {
    if has_avx2() {
        unsafe { find_byte_avx2(slice, byte) }
    } else if has_sse41() {
        unsafe { find_byte_sse41(slice, byte) }
    } else {
        find_byte_scalar(slice, byte)
    }
}

/// Compare two slices for equality using the most efficient
/// available SIMD instruction set
pub fn compare_slices(a: &[u8], b: &[u8]) -> bool {
    if a.len() != b.len() {
        return false;
    }
    
    if has_avx2() {
        unsafe { compare_slices_avx2(a, b) }
    } else if has_sse41() {
        unsafe { compare_slices_sse41(a, b) }
    } else {
        compare_slices_scalar(a, b)
    }
}

/// Convert a DNA sequence to a 2-bit packed representation using SIMD
pub fn pack_dna_sequence(src: &[u8], dst: &mut [u8]) -> usize {
    if has_avx2() {
        unsafe { pack_dna_sequence_avx2(src, dst) }
    } else if has_sse41() {
        unsafe { pack_dna_sequence_sse41(src, dst) }
    } else {
        pack_dna_sequence_scalar(src, dst)
    }
}

/// Unpack a 2-bit DNA sequence representation to ASCII using SIMD
pub fn unpack_dna_sequence(src: &[u8], dst: &mut [u8], len: usize) -> usize {
    if has_avx2() {
        unsafe { unpack_dna_sequence_avx2(src, dst, len) }
    } else if has_sse41() {
        unsafe { unpack_dna_sequence_sse41(src, dst, len) }
    } else {
        unpack_dna_sequence_scalar(src, dst, len)
    }
}

/// Scalar implementation for counting occurrences of a byte in a slice
fn count_byte_scalar(slice: &[u8], byte: u8) -> usize {
    slice.iter().filter(|&&b| b == byte).count()
}

/// Scalar implementation for finding a byte in a slice
fn find_byte_scalar(slice: &[u8], byte: u8) -> Option<usize> {
    slice.iter().position(|&b| b == byte)
}

/// Scalar implementation for comparing two slices
fn compare_slices_scalar(a: &[u8], b: &[u8]) -> bool {
    a == b
}

/// Scalar implementation for packing a DNA sequence to 2-bit representation
fn pack_dna_sequence_scalar(src: &[u8], dst: &mut [u8]) -> usize {
    let bytes_to_process = src.len();
    let bytes_required = (bytes_to_process + 3) / 4;
    
    // Make sure destination has enough space
    if dst.len() < bytes_required {
        return 0;
    }
    
    for i in 0..bytes_required {
        let mut packed_byte = 0u8;
        
        // Process 4 bases at once
        for j in 0..4 {
            let src_idx = i * 4 + j;
            if src_idx < bytes_to_process {
                let bits = match src[src_idx] {
                    b'A' | b'a' => 0b00,
                    b'C' | b'c' => 0b01,
                    b'G' | b'g' => 0b10,
                    b'T' | b't' | b'U' | b'u' => 0b11,
                    _ => 0b00, // Default to 'A' for invalid bases
                };
                packed_byte |= bits << (6 - j * 2);
            }
        }
        
        dst[i] = packed_byte;
    }
    
    bytes_required
}

/// Scalar implementation for unpacking a 2-bit DNA sequence to ASCII
fn unpack_dna_sequence_scalar(src: &[u8], dst: &mut [u8], len: usize) -> usize {
    let bases_to_unpack = len.min(dst.len());
    
    for i in 0..bases_to_unpack {
        let byte_idx = i / 4;
        let bit_offset = 6 - (i % 4) * 2;
        let bits = (src[byte_idx] >> bit_offset) & 0b11;
        
        dst[i] = match bits {
            0b00 => b'A',
            0b01 => b'C',
            0b10 => b'G',
            0b11 => b'T',
            _ => unreachable!(),
        };
    }
    
    bases_to_unpack
}

/// AVX2 implementation for counting occurrences of a byte in a slice
#[target_feature(enable = "avx2")]
unsafe fn count_byte_avx2(slice: &[u8], byte: u8) -> usize {
    let len = slice.len();
    let mut count = 0;

    if len >= 32 {
        // Broadcast byte to YMM register
        let broadcast = _mm256_set1_epi8(byte as i8);
        let mut i = 0;

        // Process 32 bytes at a time
        while i + 32 <= len {
            // Load 32 bytes
            let data = _mm256_loadu_si256(slice[i..].as_ptr() as *const __m256i);
            
            // Compare with broadcast byte
            let mask = _mm256_cmpeq_epi8(data, broadcast);
            
            // Get mask of matches
            let mask_bits = _mm256_movemask_epi8(mask) as u32;
            
            // Count set bits in mask
            count += mask_bits.count_ones() as usize;
            
            i += 32;
        }

        // Process remaining bytes with scalar method
        count += count_byte_scalar(&slice[i..], byte);
    } else {
        count = count_byte_scalar(slice, byte);
    }

    count
}

/// SSE4.1 implementation for counting occurrences of a byte in a slice
#[target_feature(enable = "sse4.1")]
unsafe fn count_byte_sse41(slice: &[u8], byte: u8) -> usize {
    let len = slice.len();
    let mut count = 0;

    if len >= 16 {
        // Broadcast byte to XMM register
        let broadcast = _mm_set1_epi8(byte as i8);
        let mut i = 0;

        // Process 16 bytes at a time
        while i + 16 <= len {
            // Load 16 bytes
            let data = _mm_loadu_si128(slice[i..].as_ptr() as *const __m128i);
            
            // Compare with broadcast byte
            let mask = _mm_cmpeq_epi8(data, broadcast);
            
            // Get mask of matches
            let mask_bits = _mm_movemask_epi8(mask) as u32;
            
            // Count set bits in mask
            count += mask_bits.count_ones() as usize;
            
            i += 16;
        }

        // Process remaining bytes with scalar method
        count += count_byte_scalar(&slice[i..], byte);
    } else {
        count = count_byte_scalar(slice, byte);
    }

    count
}

/// AVX2 implementation for finding a byte in a slice
#[target_feature(enable = "avx2")]
unsafe fn find_byte_avx2(slice: &[u8], byte: u8) -> Option<usize> {
    let len = slice.len();

    if len >= 32 {
        // Broadcast byte to YMM register
        let broadcast = _mm256_set1_epi8(byte as i8);
        let mut i = 0;

        // Process 32 bytes at a time
        while i + 32 <= len {
            // Load 32 bytes
            let data = _mm256_loadu_si256(slice[i..].as_ptr() as *const __m256i);
            
            // Compare with broadcast byte
            let mask = _mm256_cmpeq_epi8(data, broadcast);
            
            // Get mask of matches
            let mask_bits = _mm256_movemask_epi8(mask) as u32;
            
            // If there's a match, find its position
            if mask_bits != 0 {
                let pos = mask_bits.trailing_zeros() as usize;
                return Some(i + pos);
            }
            
            i += 32;
        }

        // Process remaining bytes with scalar method
        if let Some(pos) = find_byte_scalar(&slice[i..], byte) {
            return Some(i + pos);
        }
    } else {
        return find_byte_scalar(slice, byte);
    }

    None
}

/// SSE4.1 implementation for finding a byte in a slice
#[target_feature(enable = "sse4.1")]
unsafe fn find_byte_sse41(slice: &[u8], byte: u8) -> Option<usize> {
    let len = slice.len();

    if len >= 16 {
        // Broadcast byte to XMM register
        let broadcast = _mm_set1_epi8(byte as i8);
        let mut i = 0;

        // Process 16 bytes at a time
        while i + 16 <= len {
            // Load 16 bytes
            let data = _mm_loadu_si128(slice[i..].as_ptr() as *const __m128i);
            
            // Compare with broadcast byte
            let mask = _mm_cmpeq_epi8(data, broadcast);
            
            // Get mask of matches
            let mask_bits = _mm_movemask_epi8(mask) as u32;
            
            // If there's a match, find its position
            if mask_bits != 0 {
                let pos = mask_bits.trailing_zeros() as usize;
                return Some(i + pos);
            }
            
            i += 16;
        }

        // Process remaining bytes with scalar method
        if let Some(pos) = find_byte_scalar(&slice[i..], byte) {
            return Some(i + pos);
        }
    } else {
        return find_byte_scalar(slice, byte);
    }

    None
}

/// AVX2 implementation for comparing two slices
#[target_feature(enable = "avx2")]
unsafe fn compare_slices_avx2(a: &[u8], b: &[u8]) -> bool {
    let len = a.len();
    let mut i = 0;

    // Process 32 bytes at a time
    while i + 32 <= len {
        // Load 32 bytes from each slice
        let a_data = _mm256_loadu_si256(a[i..].as_ptr() as *const __m256i);
        let b_data = _mm256_loadu_si256(b[i..].as_ptr() as *const __m256i);
        
        // Compare the 32 bytes
        let mask = _mm256_cmpeq_epi8(a_data, b_data);
        
        // Get mask of matches (all bytes must match)
        let mask_bits = _mm256_movemask_epi8(mask) as u32;
        
        // If not all bits are set, slices are different
        if mask_bits != 0xFFFFFFFF {
            return false;
        }
        
        i += 32;
    }

    // Process remaining bytes with scalar method
    compare_slices_scalar(&a[i..], &b[i..])
}

/// SSE4.1 implementation for comparing two slices
#[target_feature(enable = "sse4.1")]
unsafe fn compare_slices_sse41(a: &[u8], b: &[u8]) -> bool {
    let len = a.len();
    let mut i = 0;

    // Process 16 bytes at a time
    while i + 16 <= len {
        // Load 16 bytes from each slice
        let a_data = _mm_loadu_si128(a[i..].as_ptr() as *const __m128i);
        let b_data = _mm_loadu_si128(b[i..].as_ptr() as *const __m128i);
        
        // Compare the 16 bytes
        let mask = _mm_cmpeq_epi8(a_data, b_data);
        
        // Get mask of matches (all bytes must match)
        let mask_bits = _mm_movemask_epi8(mask) as u32;
        
        // If not all bits are set, slices are different
        if mask_bits != 0xFFFF {
            return false;
        }
        
        i += 16;
    }

    // Process remaining bytes with scalar method
    compare_slices_scalar(&a[i..], &b[i..])
}

/// AVX2 implementation for packing a DNA sequence to 2-bit representation
#[target_feature(enable = "avx2")]
unsafe fn pack_dna_sequence_avx2(src: &[u8], dst: &mut [u8]) -> usize {
    // For simplicity, delegate to scalar implementation for now
    // In a real implementation, you would optimize this with AVX2 instructions
    pack_dna_sequence_scalar(src, dst)
}

/// SSE4.1 implementation for packing a DNA sequence to 2-bit representation
#[target_feature(enable = "sse4.1")]
unsafe fn pack_dna_sequence_sse41(src: &[u8], dst: &mut [u8]) -> usize {
    // For simplicity, delegate to scalar implementation for now
    // In a real implementation, you would optimize this with SSE4.1 instructions
    pack_dna_sequence_scalar(src, dst)
}

/// AVX2 implementation for unpacking a 2-bit DNA sequence to ASCII
#[target_feature(enable = "avx2")]
unsafe fn unpack_dna_sequence_avx2(src: &[u8], dst: &mut [u8], len: usize) -> usize {
    // For simplicity, delegate to scalar implementation for now
    // In a real implementation, you would optimize this with AVX2 instructions
    unpack_dna_sequence_scalar(src, dst, len)
}

/// SSE4.1 implementation for unpacking a 2-bit DNA sequence to ASCII
#[target_feature(enable = "sse4.1")]
unsafe fn unpack_dna_sequence_sse41(src: &[u8], dst: &mut [u8], len: usize) -> usize {
    // For simplicity, delegate to scalar implementation for now
    // In a real implementation, you would optimize this with SSE4.1 instructions
    unpack_dna_sequence_scalar(src, dst, len)
}

#[cfg(test)]
mod tests {
    use super::*;
    
    #[test]
    fn test_count_byte() {
        let data = b"ACGTACGTACGTACGTACGTACGTACGTACGT";
        
        // Count As
        assert_eq!(count_byte(data, b'A'), 8);
        
        // Count Cs
        assert_eq!(count_byte(data, b'C'), 8);
        
        // Count Ts
        assert_eq!(count_byte(data, b'T'), 8);
        
        // Count Gs
        assert_eq!(count_byte(data, b'G'), 8);
        
        // Count a byte not in the array
        assert_eq!(count_byte(data, b'N'), 0);
    }
    
    #[test]
    fn test_find_byte() {
        let data = b"ACGTACGTACGTACGTACGTACGTACGTACGT";
        
        // Find first A
        assert_eq!(find_byte(data, b'A'), Some(0));
        
        // Find first C
        assert_eq!(find_byte(data, b'C'), Some(1));
        
        // Find first G
        assert_eq!(find_byte(data, b'G'), Some(2));
        
        // Find first T
        assert_eq!(find_byte(data, b'T'), Some(3));
        
        // Find a byte not in the array
        assert_eq!(find_byte(data, b'N'), None);
    }
    
    #[test]
    fn test_compare_slices() {
        let a = b"ACGTACGTACGTACGTACGTACGTACGTACGT";
        let b = b"ACGTACGTACGTACGTACGTACGTACGTACGT";
        let c = b"ACGTACGTACGTACGTACGTACGTACGTACGA";
        
        // Equal slices
        assert!(compare_slices(a, b));
        
        // Different slices
        assert!(!compare_slices(a, c));
        
        // Different length slices
        assert!(!compare_slices(a, &c[0..30]));
    }
    
    #[test]
    fn test_pack_unpack_dna() {
        let dna = b"ACGTACGTACGTACGTACGTACGTACGTACGT";
        let mut packed = vec![0u8; (dna.len() + 3) / 4];
        let mut unpacked = vec![0u8; dna.len()];
        
        // Pack DNA sequence
        let packed_size = pack_dna_sequence(dna, &mut packed);
        assert_eq!(packed_size, (dna.len() + 3) / 4);
        
        // Unpack DNA sequence
        let unpacked_size = unpack_dna_sequence(&packed, &mut unpacked, dna.len());
        assert_eq!(unpacked_size, dna.len());
        
        // Check that unpacked sequence matches original
        assert_eq!(&unpacked, dna);
    }
}