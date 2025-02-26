//! Storage and data format handling
//!
//! This module provides efficient storage and parsing for various
//! bioinformatics file formats.

pub mod formats;

use crate::engines::EngineResult;
use crate::engines::core::memory::MemoryMapped;
use std::path::Path;

/// Storage mode for sequence data
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum StorageMode {
    /// In-memory storage (full sequence loaded in RAM)
    InMemory,
    /// Memory-mapped storage (sequence accessed from disk via memory mapping)
    MemoryMapped,
    /// On-demand loading (sequence loaded in chunks as needed)
    OnDemand,
}

impl Default for StorageMode {
    fn default() -> Self {
        StorageMode::InMemory
    }
}

/// Trait for storable sequence data
pub trait StorableSequence: Send + Sync {
    /// Get the length of the sequence
    fn len(&self) -> usize;
    
    /// Check if the sequence is empty
    fn is_empty(&self) -> bool {
        self.len() == 0
    }
    
    /// Get a subsequence as a vector
    fn subsequence(&self, start: usize, end: usize) -> Vec<u8>;
    
    /// Get a slice of the sequence
    fn as_slice(&self) -> Option<&[u8]>;
    
    /// Get the storage mode
    fn storage_mode(&self) -> StorageMode;
    
    /// Get the memory usage of the sequence storage
    fn memory_usage(&self) -> usize;
}

/// In-memory sequence storage
#[derive(Debug, Clone)]
pub struct InMemoryStorage {
    /// The full sequence data
    data: Vec<u8>,
}

impl InMemoryStorage {
    /// Create a new in-memory storage
    pub fn new(data: Vec<u8>) -> Self {
        Self { data }
    }
    
    /// Create a new in-memory storage with the given capacity
    pub fn with_capacity(capacity: usize) -> Self {
        Self {
            data: Vec::with_capacity(capacity),
        }
    }
    
    /// Get a reference to the underlying data
    pub fn data(&self) -> &Vec<u8> {
        &self.data
    }
    
    /// Get a mutable reference to the underlying data
    pub fn data_mut(&mut self) -> &mut Vec<u8> {
        &mut self.data
    }
}

impl StorableSequence for InMemoryStorage {
    fn len(&self) -> usize {
        self.data.len()
    }
    
    fn subsequence(&self, start: usize, end: usize) -> Vec<u8> {
        let start = start.min(self.data.len());
        let end = end.min(self.data.len());
        self.data[start..end].to_vec()
    }
    
    fn as_slice(&self) -> Option<&[u8]> {
        Some(&self.data)
    }
    
    fn storage_mode(&self) -> StorageMode {
        StorageMode::InMemory
    }
    
    fn memory_usage(&self) -> usize {
        self.data.capacity()
    }
}

/// Memory-mapped sequence storage
#[derive(Debug)]
pub struct MemoryMappedStorage {
    /// The memory mapped file
    mmap: MemoryMapped,
    /// The length of the sequence
    length: usize,
}

impl MemoryMappedStorage {
    /// Create a new memory-mapped storage from a file
    pub fn new<P: AsRef<Path>>(path: P) -> EngineResult<Self> {
        let mmap = MemoryMapped::new(
            path,
            crate::engines::core::memory::MemoryMapMode::ReadOnly,
        )?;
        let length = mmap.len();
        
        Ok(Self { mmap, length })
    }
}

impl StorableSequence for MemoryMappedStorage {
    fn len(&self) -> usize {
        self.length
    }
    
    fn subsequence(&self, start: usize, end: usize) -> Vec<u8> {
        let start = start.min(self.length);
        let end = end.min(self.length);
        self.mmap.slice(start, end).to_vec()
    }
    
    fn as_slice(&self) -> Option<&[u8]> {
        Some(self.mmap.as_slice())
    }
    
    fn storage_mode(&self) -> StorageMode {
        StorageMode::MemoryMapped
    }
    
    fn memory_usage(&self) -> usize {
        // Only count metadata, not the mapped file
        std::mem::size_of::<Self>()
    }
}

/// Chunked on-demand sequence storage
pub struct OnDemandStorage {
    /// The path to the file
    path: String,
    /// The length of the sequence
    length: usize,
    /// The chunk size for loading
    chunk_size: usize,
    /// Currently loaded chunk
    current_chunk: Option<(usize, Vec<u8>)>,
}

impl OnDemandStorage {
    /// Create a new on-demand storage
    pub fn new<P: AsRef<Path>>(path: P, length: usize, chunk_size: usize) -> EngineResult<Self> {
        Ok(Self {
            path: path.as_ref().to_string_lossy().to_string(),
            length,
            chunk_size,
            current_chunk: None,
        })
    }
    
    /// Load a chunk containing the given position
    fn load_chunk(&mut self, position: usize) -> EngineResult<()> {
        // Check if the position is already in the current chunk
        if let Some((start, ref chunk)) = self.current_chunk {
            let end = start + chunk.len();
            if position >= start && position < end {
                return Ok(());
            }
        }
        
        // Calculate the chunk to load
        let chunk_start = (position / self.chunk_size) * self.chunk_size;
        let chunk_end = (chunk_start + self.chunk_size).min(self.length);
        
        // Load the chunk from the file
        let reader = crate::engines::core::io::FastReader::new(&self.path, None)?;
        let mut buffer = vec![0; chunk_end - chunk_start];
        let _ = reader.read_chunk(&mut buffer)?;
        
        // Store the loaded chunk
        self.current_chunk = Some((chunk_start, buffer));
        
        Ok(())
    }
}

impl StorableSequence for OnDemandStorage {
    fn len(&self) -> usize {
        self.length
    }
    
    fn subsequence(&self, start: usize, end: usize) -> Vec<u8> {
        let start = start.min(self.length);
        let end = end.min(self.length);
        let mut result = Vec::with_capacity(end - start);
        
        // Need a mutable reference to load chunks
        let mut storage = self.clone();
        
        // Load and copy each chunk that contains the requested subsequence
        let mut pos = start;
        while pos < end {
            let chunk_start = (pos / storage.chunk_size) * storage.chunk_size;
            let chunk_end = (chunk_start + storage.chunk_size).min(storage.length);
            
            if let Ok(()) = storage.load_chunk(pos) {
                if let Some((chunk_pos, ref chunk)) = storage.current_chunk {
                    let offset = pos - chunk_pos;
                    let copy_end = (end - chunk_pos).min(chunk.len());
                    result.extend_from_slice(&chunk[offset..copy_end]);
                    pos += copy_end - offset;
                }
            } else {
                // Error loading chunk, fill with placeholder value
                let remaining = end - pos;
                result.extend(vec![b'N'; remaining]);
                break;
            }
        }
        
        result
    }
    
    fn as_slice(&self) -> Option<&[u8]> {
        // On-demand storage doesn't provide direct slice access
        None
    }
    
    fn storage_mode(&self) -> StorageMode {
        StorageMode::OnDemand
    }
    
    fn memory_usage(&self) -> usize {
        let chunk_size = match &self.current_chunk {
            Some((_, chunk)) => chunk.capacity(),
            None => 0,
        };
        
        // Count metadata and currently loaded chunk
        std::mem::size_of::<Self>() + chunk_size
    }
}

impl Clone for OnDemandStorage {
    fn clone(&self) -> Self {
        Self {
            path: self.path.clone(),
            length: self.length,
            chunk_size: self.chunk_size,
            current_chunk: self.current_chunk.clone(),
        }
    }
}

/// Factory for creating appropriate storage based on sequence size and preferences
pub struct StorageFactory;

impl StorageFactory {
    /// Create the most appropriate storage for a sequence with the given parameters
    pub fn create_storage(
        data: Option<Vec<u8>>,
        path: Option<&Path>,
        length: Option<usize>,
        preferred_mode: Option<StorageMode>,
    ) -> EngineResult<Box<dyn StorableSequence>> {
        // Use the preferred mode if specified
        let mode = preferred_mode.unwrap_or_else(|| {
            if let Some(len) = length {
                // Use memory mapping for sequences > 100MB
                if len > 100 * 1024 * 1024 {
                    StorageMode::MemoryMapped
                // Use on-demand for sequences > 10MB
                } else if len > 10 * 1024 * 1024 {
                    StorageMode::OnDemand
                } else {
                    StorageMode::InMemory
                }
            } else {
                StorageMode::InMemory
            }
        });
        
        match mode {
            StorageMode::InMemory => {
                if let Some(data) = data {
                    Ok(Box::new(InMemoryStorage::new(data)))
                } else if let (Some(p), Some(len)) = (path, length) {
                    // Read the file into memory
                    let mut reader = crate::engines::core::io::FastReader::new(p, None)?;
                    let data = reader.read_all()?;
                    Ok(Box::new(InMemoryStorage::new(data)))
                } else {
                    Err(crate::engines::EngineError::InvalidSequenceData(
                        "Cannot create in-memory storage without data or path".to_string(),
                    ))
                }
            },
            StorageMode::MemoryMapped => {
                if let Some(p) = path {
                    let storage = MemoryMappedStorage::new(p)?;
                    Ok(Box::new(storage))
                } else {
                    Err(crate::engines::EngineError::InvalidSequenceData(
                        "Cannot create memory-mapped storage without a path".to_string(),
                    ))
                }
            },
            StorageMode::OnDemand => {
                if let (Some(p), Some(len)) = (path, length) {
                    let chunk_size = 1024 * 1024; // 1MB chunks
                    let storage = OnDemandStorage::new(p, len, chunk_size)?;
                    Ok(Box::new(storage))
                } else {
                    Err(crate::engines::EngineError::InvalidSequenceData(
                        "Cannot create on-demand storage without path and length".to_string(),
                    ))
                }
            },
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;
    use tempfile::tempdir;
    
    #[test]
    fn test_in_memory_storage() {
        let data = b"ACGTACGTACGT".to_vec();
        let storage = InMemoryStorage::new(data.clone());
        
        // Check length
        assert_eq!(storage.len(), data.len());
        
        // Check subsequence
        let sub = storage.subsequence(4, 8);
        assert_eq!(sub, b"ACGT");
        
        // Check slice access
        assert_eq!(storage.as_slice().unwrap(), data.as_slice());
        
        // Check storage mode
        assert_eq!(storage.storage_mode(), StorageMode::InMemory);
    }
    
    #[test]
    fn test_memory_mapped_storage() -> std::io::Result<()> {
        // Create a temporary file
        let dir = tempdir()?;
        let file_path = dir.path().join("test.seq");
        
        // Write test data
        let data = b"ACGTACGTACGT";
        {
            let mut file = std::fs::File::create(&file_path)?;
            file.write_all(data)?;
        }
        
        // Create memory-mapped storage
        let storage = MemoryMappedStorage::new(&file_path).unwrap();
        
        // Check length
        assert_eq!(storage.len(), data.len());
        
        // Check subsequence
        let sub = storage.subsequence(4, 8);
        assert_eq!(sub, b"ACGT");
        
        // Check slice access
        assert_eq!(storage.as_slice().unwrap(), data);
        
        // Check storage mode
        assert_eq!(storage.storage_mode(), StorageMode::MemoryMapped);
        
        Ok(())
    }
    
    #[test]
    fn test_storage_factory() -> std::io::Result<()> {
        // Create a temporary file
        let dir = tempdir()?;
        let file_path = dir.path().join("test.seq");
        
        // Write test data
        let data = b"ACGTACGTACGT";
        {
            let mut file = std::fs::File::create(&file_path)?;
            file.write_all(data)?;
        }
        
        // Test in-memory storage creation
        let storage = StorageFactory::create_storage(
            Some(data.to_vec()),
            None,
            None,
            Some(StorageMode::InMemory),
        ).unwrap();
        
        assert_eq!(storage.storage_mode(), StorageMode::InMemory);
        assert_eq!(storage.len(), data.len());
        
        // Test memory-mapped storage creation
        let storage = StorageFactory::create_storage(
            None,
            Some(&file_path),
            Some(data.len()),
            Some(StorageMode::MemoryMapped),
        ).unwrap();
        
        assert_eq!(storage.storage_mode(), StorageMode::MemoryMapped);
        assert_eq!(storage.len(), data.len());
        
        Ok(())
    }
}