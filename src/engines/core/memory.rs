//! Memory optimization primitives
//!
//! This module provides optimized memory allocation, mapping, and
//! management for biological sequence data.

use std::fs::File;
use std::path::Path;
use memmap2::{Mmap, MmapOptions};
use std::sync::atomic::{AtomicUsize, Ordering};

/// Alignment for memory allocations (in bytes)
/// Set to 64 for optimal cache line alignment on most processors
pub const MEMORY_ALIGNMENT: usize = 64;

/// Default page size for memory operations
pub const DEFAULT_PAGE_SIZE: usize = 4096;

// Track memory usage statistics
static TOTAL_ALLOCATED: AtomicUsize = AtomicUsize::new(0);
static PEAK_ALLOCATED: AtomicUsize = AtomicUsize::new(0);

/// Initialize the memory subsystem
pub fn initialize() {
    // Reset memory usage counters
    TOTAL_ALLOCATED.store(0, Ordering::SeqCst);
    PEAK_ALLOCATED.store(0, Ordering::SeqCst);
}

/// Memory mapping mode
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum MemoryMapMode {
    /// Read-only mapping
    ReadOnly,
    /// Read-write mapping
    ReadWrite,
    /// Copy-on-write mapping
    CopyOnWrite,
}

/// Memory-mapped file for efficient large sequence storage
pub struct MemoryMapped {
    mmap: Mmap,
    len: usize,
}

impl MemoryMapped {
    /// Create a new memory-mapped file
    pub fn new<P: AsRef<Path>>(path: P, mode: MemoryMapMode) -> std::io::Result<Self> {
        let file = File::open(path)?;
        let len = file.metadata()?.len() as usize;
        
        let mmap = match mode {
            MemoryMapMode::ReadOnly => unsafe { MmapOptions::new().map(&file)? },
            MemoryMapMode::ReadWrite => unsafe { MmapOptions::new().map_mut(&file)?.into() },
            MemoryMapMode::CopyOnWrite => unsafe { MmapOptions::new().map_copy(&file)? },
        };
        
        Ok(Self { mmap, len })
    }
    
    /// Get a reference to the underlying memory-mapped data
    pub fn as_slice(&self) -> &[u8] {
        &self.mmap[0..self.len]
    }
    
    /// Get the length of the memory-mapped data
    pub fn len(&self) -> usize {
        self.len
    }
    
    /// Check if the memory-mapped data is empty
    pub fn is_empty(&self) -> bool {
        self.len == 0
    }
    
    /// Get a slice of the memory-mapped data
    pub fn slice(&self, start: usize, end: usize) -> &[u8] {
        &self.mmap[start.min(self.len)..end.min(self.len)]
    }
}

/// Packed 2-bit encoding for DNA sequences
pub struct PackedDnaStorage {
    /// The packed sequence data (2 bits per base)
    data: Vec<u8>,
    /// The length of the sequence in bases
    len: usize,
}

impl PackedDnaStorage {
    /// Create a new packed DNA storage with the given capacity
    pub fn with_capacity(capacity: usize) -> Self {
        // Calculate number of bytes needed for 2-bit encoding
        let byte_capacity = (capacity + 3) / 4;
        Self {
            data: Vec::with_capacity(byte_capacity),
            len: 0,
        }
    }
    
    /// Pack a DNA sequence into the storage
    /// A=00, C=01, G=10, T=11
    pub fn pack(&mut self, sequence: &[u8]) {
        // Calculate required capacity
        let required_bytes = (sequence.len() + 3) / 4;
        self.data.clear();
        self.data.reserve(required_bytes);
        
        // Pack 4 bases per byte
        for chunk in sequence.chunks(4) {
            let mut byte = 0u8;
            for (i, &base) in chunk.iter().enumerate() {
                let bits = match base {
                    b'A' | b'a' => 0b00,
                    b'C' | b'c' => 0b01,
                    b'G' | b'g' => 0b10,
                    b'T' | b't' | b'U' | b'u' => 0b11,
                    _ => 0b00, // Default to 'A' for invalid bases
                };
                byte |= bits << (6 - i * 2);
            }
            self.data.push(byte);
        }
        
        self.len = sequence.len();
        
        // Update memory tracking
        update_memory_usage(self.data.capacity());
    }
    
    /// Unpack the DNA sequence into the provided buffer
    pub fn unpack(&self, buffer: &mut [u8]) -> usize {
        let unpack_len = self.len.min(buffer.len());
        
        for i in 0..unpack_len {
            let byte_idx = i / 4;
            let bit_offset = 6 - (i % 4) * 2;
            let bits = (self.data[byte_idx] >> bit_offset) & 0b11;
            
            buffer[i] = match bits {
                0b00 => b'A',
                0b01 => b'C',
                0b10 => b'G',
                0b11 => b'T',
                _ => unreachable!(),
            };
        }
        
        unpack_len
    }
    
    /// Get the length of the sequence in bases
    pub fn len(&self) -> usize {
        self.len
    }
    
    /// Check if the sequence is empty
    pub fn is_empty(&self) -> bool {
        self.len == 0
    }
    
    /// Get the memory usage in bytes
    pub fn memory_usage(&self) -> usize {
        self.data.capacity()
    }
}

/// Packed 4-bit encoding for protein sequences
pub struct PackedProteinStorage {
    /// The packed sequence data (4 bits per amino acid)
    data: Vec<u8>,
    /// The length of the sequence in amino acids
    len: usize,
}

impl PackedProteinStorage {
    /// Create a new packed protein storage with the given capacity
    pub fn with_capacity(capacity: usize) -> Self {
        // Calculate number of bytes needed for 4-bit encoding
        let byte_capacity = (capacity + 1) / 2;
        Self {
            data: Vec::with_capacity(byte_capacity),
            len: 0,
        }
    }
    
    /// Pack a protein sequence into the storage (4 bits per amino acid)
    pub fn pack(&mut self, sequence: &[u8]) {
        // Calculate required capacity
        let required_bytes = (sequence.len() + 1) / 2;
        self.data.clear();
        self.data.reserve(required_bytes);
        
        // Pack 2 amino acids per byte
        for chunk in sequence.chunks(2) {
            let mut byte = 0u8;
            
            // First amino acid in high nibble
            let first = encode_amino_acid(chunk[0]);
            byte |= first << 4;
            
            // Second amino acid in low nibble (if exists)
            if chunk.len() > 1 {
                let second = encode_amino_acid(chunk[1]);
                byte |= second;
            }
            
            self.data.push(byte);
        }
        
        self.len = sequence.len();
        
        // Update memory tracking
        update_memory_usage(self.data.capacity());
    }
    
    /// Unpack the protein sequence into the provided buffer
    pub fn unpack(&self, buffer: &mut [u8]) -> usize {
        let unpack_len = self.len.min(buffer.len());
        
        for i in 0..unpack_len {
            let byte_idx = i / 2;
            let is_high_nibble = i % 2 == 0;
            
            let nibble = if is_high_nibble {
                (self.data[byte_idx] >> 4) & 0x0F
            } else {
                self.data[byte_idx] & 0x0F
            };
            
            buffer[i] = decode_amino_acid(nibble);
        }
        
        unpack_len
    }
    
    /// Get the length of the sequence in amino acids
    pub fn len(&self) -> usize {
        self.len
    }
    
    /// Check if the sequence is empty
    pub fn is_empty(&self) -> bool {
        self.len == 0
    }
    
    /// Get the memory usage in bytes
    pub fn memory_usage(&self) -> usize {
        self.data.capacity()
    }
}

/// Update memory usage statistics
fn update_memory_usage(bytes: usize) {
    let current = TOTAL_ALLOCATED.fetch_add(bytes, Ordering::SeqCst) + bytes;
    let mut peak = PEAK_ALLOCATED.load(Ordering::SeqCst);
    
    while current > peak {
        match PEAK_ALLOCATED.compare_exchange_weak(
            peak, current, Ordering::SeqCst, Ordering::SeqCst
        ) {
            Ok(_) => break,
            Err(new_peak) => peak = new_peak,
        }
    }
}

/// Get current memory usage statistics
pub fn get_memory_stats() -> (usize, usize) {
    (
        TOTAL_ALLOCATED.load(Ordering::SeqCst),
        PEAK_ALLOCATED.load(Ordering::SeqCst)
    )
}

/// Encode an amino acid to a 4-bit representation
fn encode_amino_acid(aa: u8) -> u8 {
    match aa {
        b'A' | b'a' => 0x0,
        b'R' | b'r' => 0x1,
        b'N' | b'n' => 0x2,
        b'D' | b'd' => 0x3,
        b'C' | b'c' => 0x4,
        b'Q' | b'q' => 0x5,
        b'E' | b'e' => 0x6,
        b'G' | b'g' => 0x7,
        b'H' | b'h' => 0x8,
        b'I' | b'i' => 0x9,
        b'L' | b'l' => 0xA,
        b'K' | b'k' => 0xB,
        b'M' | b'm' => 0xC,
        b'F' | b'f' => 0xD,
        b'P' | b'p' => 0xE,
        _ => 0xF, // All other amino acids and stop codons
    }
}

/// Decode a 4-bit representation to an amino acid
fn decode_amino_acid(nibble: u8) -> u8 {
    match nibble {
        0x0 => b'A',
        0x1 => b'R',
        0x2 => b'N',
        0x3 => b'D',
        0x4 => b'C',
        0x5 => b'Q',
        0x6 => b'E',
        0x7 => b'G',
        0x8 => b'H',
        0x9 => b'I',
        0xA => b'L',
        0xB => b'K',
        0xC => b'M',
        0xD => b'F',
        0xE => b'P',
        _ => b'X', // Unknown or stop codon
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    
    #[test]
    fn test_dna_packing() {
        let dna = b"ACGTACGT";
        let mut packed = PackedDnaStorage::with_capacity(dna.len());
        packed.pack(dna);
        
        // Check length
        assert_eq!(packed.len(), dna.len());
        
        // Check unpacking
        let mut buffer = vec![0; dna.len()];
        let unpacked_len = packed.unpack(&mut buffer);
        assert_eq!(unpacked_len, dna.len());
        assert_eq!(&buffer, dna);
        
        // Check memory usage
        assert_eq!(packed.memory_usage(), 2); // 8 bases = 2 bytes
    }
    
    #[test]
    fn test_protein_packing() {
        let protein = b"ARNDCQEGHILKMFP";
        let mut packed = PackedProteinStorage::with_capacity(protein.len());
        packed.pack(protein);
        
        // Check length
        assert_eq!(packed.len(), protein.len());
        
        // Check unpacking
        let mut buffer = vec![0; protein.len()];
        let unpacked_len = packed.unpack(&mut buffer);
        assert_eq!(unpacked_len, protein.len());
        assert_eq!(&buffer, protein);
        
        // Check memory usage
        assert_eq!(packed.memory_usage(), 8); // 15 amino acids = 8 bytes
    }
    
    #[test]
    fn test_memory_stats() {
        // Reset counters
        initialize();
        
        // Create some allocations
        let mut dna = PackedDnaStorage::with_capacity(1000);
        dna.pack(b"ACGT"); // Small allocation
        
        let (total, peak) = get_memory_stats();
        assert!(total > 0);
        assert_eq!(total, peak);
    }
}