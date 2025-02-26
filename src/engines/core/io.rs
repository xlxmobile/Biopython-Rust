//! High-performance I/O operations
//!
//! This module provides optimized I/O operations for biological sequence data,
//! focusing on efficient reading and writing of large sequence files.

use std::fs::{File, OpenOptions};
use std::io::{self, Read, Write, BufReader, BufWriter, SeekFrom, Seek};
use std::path::Path;
use memmap2::{Mmap, MmapOptions};
use std::sync::atomic::{AtomicUsize, Ordering};
use crate::engines::core::memory::MemoryMapped;

// Default buffer sizes
const DEFAULT_READ_BUFFER_SIZE: usize = 1024 * 1024; // 1MB
const DEFAULT_WRITE_BUFFER_SIZE: usize = 1024 * 1024; // 1MB

// Tracking I/O statistics
static TOTAL_BYTES_READ: AtomicUsize = AtomicUsize::new(0);
static TOTAL_BYTES_WRITTEN: AtomicUsize = AtomicUsize::new(0);

/// Initialize the I/O subsystem
pub fn initialize() {
    // Reset I/O counters
    TOTAL_BYTES_READ.store(0, Ordering::SeqCst);
    TOTAL_BYTES_WRITTEN.store(0, Ordering::SeqCst);
}

/// High-performance buffered file reader
pub struct FastReader {
    reader: BufReader<File>,
    path: String,
    buffer_size: usize,
}

impl FastReader {
    /// Create a new fast reader for the given file path
    pub fn new<P: AsRef<Path>>(path: P, buffer_size: Option<usize>) -> io::Result<Self> {
        let file = File::open(path.as_ref())?;
        let buf_size = buffer_size.unwrap_or(DEFAULT_READ_BUFFER_SIZE);
        let reader = BufReader::with_capacity(buf_size, file);
        
        Ok(Self {
            reader,
            path: path.as_ref().to_string_lossy().to_string(),
            buffer_size: buf_size,
        })
    }
    
    /// Read the entire file into a vector
    pub fn read_all(&mut self) -> io::Result<Vec<u8>> {
        let mut buffer = Vec::new();
        let bytes_read = self.reader.read_to_end(&mut buffer)?;
        
        // Update read statistics
        TOTAL_BYTES_READ.fetch_add(bytes_read, Ordering::SeqCst);
        
        Ok(buffer)
    }
    
    /// Read a chunk of data from the file
    pub fn read_chunk(&mut self, buffer: &mut [u8]) -> io::Result<usize> {
        let bytes_read = self.reader.read(buffer)?;
        
        // Update read statistics
        TOTAL_BYTES_READ.fetch_add(bytes_read, Ordering::SeqCst);
        
        Ok(bytes_read)
    }
    
    /// Read the file line by line (efficiency optimized)
    pub fn read_lines(&mut self) -> Lines<'_> {
        Lines {
            reader: &mut self.reader,
            buffer: String::new(),
        }
    }
    
    /// Get the path of the file being read
    pub fn path(&self) -> &str {
        &self.path
    }
    
    /// Get the buffer size being used
    pub fn buffer_size(&self) -> usize {
        self.buffer_size
    }
    
    /// Reset the reader to the beginning of the file
    pub fn reset(&mut self) -> io::Result<()> {
        self.reader.seek(SeekFrom::Start(0))?;
        Ok(())
    }
}

/// Iterator over lines in a file
pub struct Lines<'a> {
    reader: &'a mut BufReader<File>,
    buffer: String,
}

impl<'a> Iterator for Lines<'a> {
    type Item = io::Result<String>;
    
    fn next(&mut self) -> Option<Self::Item> {
        self.buffer.clear();
        match self.reader.read_line(&mut self.buffer) {
            Ok(0) => None, // EOF
            Ok(bytes) => {
                // Update read statistics
                TOTAL_BYTES_READ.fetch_add(bytes, Ordering::SeqCst);
                
                // Trim the trailing newline
                if self.buffer.ends_with('\n') {
                    self.buffer.pop();
                    if self.buffer.ends_with('\r') {
                        self.buffer.pop();
                    }
                }
                
                Some(Ok(self.buffer.clone()))
            }
            Err(e) => Some(Err(e)),
        }
    }
}

/// High-performance buffered file writer
pub struct FastWriter {
    writer: BufWriter<File>,
    path: String,
    buffer_size: usize,
}

impl FastWriter {
    /// Create a new fast writer for the given file path
    pub fn new<P: AsRef<Path>>(path: P, buffer_size: Option<usize>) -> io::Result<Self> {
        let file = OpenOptions::new()
            .write(true)
            .create(true)
            .truncate(true)
            .open(path.as_ref())?;
        
        let buf_size = buffer_size.unwrap_or(DEFAULT_WRITE_BUFFER_SIZE);
        let writer = BufWriter::with_capacity(buf_size, file);
        
        Ok(Self {
            writer,
            path: path.as_ref().to_string_lossy().to_string(),
            buffer_size: buf_size,
        })
    }
    
    /// Append to an existing file instead of overwriting
    pub fn append<P: AsRef<Path>>(path: P, buffer_size: Option<usize>) -> io::Result<Self> {
        let file = OpenOptions::new()
            .write(true)
            .create(true)
            .append(true)
            .open(path.as_ref())?;
        
        let buf_size = buffer_size.unwrap_or(DEFAULT_WRITE_BUFFER_SIZE);
        let writer = BufWriter::with_capacity(buf_size, file);
        
        Ok(Self {
            writer,
            path: path.as_ref().to_string_lossy().to_string(),
            buffer_size: buf_size,
        })
    }
    
    /// Write data to the file
    pub fn write(&mut self, data: &[u8]) -> io::Result<usize> {
        let bytes_written = self.writer.write(data)?;
        
        // Update write statistics
        TOTAL_BYTES_WRITTEN.fetch_add(bytes_written, Ordering::SeqCst);
        
        Ok(bytes_written)
    }
    
    /// Write a line to the file (appends a newline)
    pub fn write_line(&mut self, line: &str) -> io::Result<usize> {
        let bytes_written = self.writer.write(line.as_bytes())?;
        let newline_written = self.writer.write(b"\n")?;
        
        // Update write statistics
        TOTAL_BYTES_WRITTEN.fetch_add(bytes_written + newline_written, Ordering::SeqCst);
        
        Ok(bytes_written + newline_written)
    }
    
    /// Flush any buffered data to disk
    pub fn flush(&mut self) -> io::Result<()> {
        self.writer.flush()
    }
    
    /// Get the path of the file being written
    pub fn path(&self) -> &str {
        &self.path
    }
    
    /// Get the buffer size being used
    pub fn buffer_size(&self) -> usize {
        self.buffer_size
    }
}

/// Memory-mapped sequence file reader for efficient processing of large files
pub struct MemoryMappedReader {
    mmap: Mmap,
    path: String,
    position: usize,
}

impl MemoryMappedReader {
    /// Create a new memory-mapped reader for the given file path
    pub fn new<P: AsRef<Path>>(path: P) -> io::Result<Self> {
        let file = File::open(path.as_ref())?;
        let mmap = unsafe { MmapOptions::new().map(&file)? };
        
        Ok(Self {
            mmap,
            path: path.as_ref().to_string_lossy().to_string(),
            position: 0,
        })
    }
    
    /// Get a slice of the entire memory-mapped file
    pub fn as_slice(&self) -> &[u8] {
        &self.mmap[..]
    }
    
    /// Get the length of the memory-mapped file
    pub fn len(&self) -> usize {
        self.mmap.len()
    }
    
    /// Check if the memory-mapped file is empty
    pub fn is_empty(&self) -> bool {
        self.mmap.is_empty()
    }
    
    /// Get the path of the file being read
    pub fn path(&self) -> &str {
        &self.path
    }
    
    /// Get a slice of the memory-mapped file from the current position
    pub fn current_slice(&self, length: usize) -> &[u8] {
        let end = (self.position + length).min(self.mmap.len());
        &self.mmap[self.position..end]
    }
    
    /// Seek to a position in the memory-mapped file
    pub fn seek(&mut self, pos: usize) {
        self.position = pos.min(self.mmap.len());
    }
    
    /// Advance the current position by the given amount
    pub fn advance(&mut self, amount: usize) {
        self.position = (self.position + amount).min(self.mmap.len());
    }
    
    /// Get the current position in the memory-mapped file
    pub fn position(&self) -> usize {
        self.position
    }
    
    /// Check if we've reached the end of the file
    pub fn is_eof(&self) -> bool {
        self.position >= self.mmap.len()
    }
}

/// Split a file into chunks for parallel processing
pub fn split_file_into_chunks<P: AsRef<Path>>(
    path: P,
    chunk_size: usize,
) -> io::Result<Vec<(usize, usize)>> {
    let file = File::open(path.as_ref())?;
    let file_size = file.metadata()?.len() as usize;
    
    let mut chunks = Vec::new();
    let mut start = 0;
    
    while start < file_size {
        let end = (start + chunk_size).min(file_size);
        chunks.push((start, end));
        start = end;
    }
    
    Ok(chunks)
}

/// Process a file in parallel using memory-mapped I/O
pub fn process_file_parallel<P, F, R>(
    path: P,
    chunk_size: usize,
    processor: F,
) -> io::Result<Vec<R>>
where
    P: AsRef<Path>,
    F: Fn(&[u8]) -> R + Send + Sync + 'static,
    R: Send + 'static,
{
    // Create memory map
    let mmap = MemoryMapped::new(path, crate::engines::core::memory::MemoryMapMode::ReadOnly)?;
    let data = mmap.as_slice();
    
    // Split into chunks
    let chunk_bounds = crate::engines::core::parallel::chunk_slice(data, Some(chunk_size));
    
    // Process chunks in parallel
    let processor = &processor;
    let results = crate::engines::core::parallel::execute(|pool| {
        pool.install(|| {
            chunk_bounds
                .par_iter()
                .map(|chunk| processor(chunk))
                .collect::<Vec<R>>()
        })
    });
    
    Ok(results)
}

/// Get the current I/O statistics
pub fn get_io_stats() -> (usize, usize) {
    (
        TOTAL_BYTES_READ.load(Ordering::SeqCst),
        TOTAL_BYTES_WRITTEN.load(Ordering::SeqCst),
    )
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;
    use tempfile::tempdir;
    
    #[test]
    fn test_fast_reader_writer() -> io::Result<()> {
        // Create a temporary directory
        let dir = tempdir()?;
        let file_path = dir.path().join("test.txt");
        
        // Write test data
        let test_data = b"Line 1\nLine 2\nLine 3";
        {
            let mut writer = FastWriter::new(&file_path, None)?;
            writer.write(test_data)?;
            writer.flush()?;
        }
        
        // Read test data
        let mut reader = FastReader::new(&file_path, None)?;
        let data = reader.read_all()?;
        assert_eq!(data, test_data);
        
        // Test reading lines
        reader.reset()?;
        let lines: Result<Vec<String>, _> = reader.read_lines().collect();
        let lines = lines?;
        assert_eq!(lines, vec!["Line 1", "Line 2", "Line 3"]);
        
        // Check I/O statistics
        let (bytes_read, bytes_written) = get_io_stats();
        assert_eq!(bytes_written, test_data.len());
        assert_eq!(bytes_read, test_data.len() * 2); // read_all + read_lines
        
        Ok(())
    }
    
    #[test]
    fn test_memory_mapped_reader() -> io::Result<()> {
        // Create a temporary directory
        let dir = tempdir()?;
        let file_path = dir.path().join("test_mmap.txt");
        
        // Write test data
        let test_data = b"Memory mapped test data";
        {
            let mut file = File::create(&file_path)?;
            file.write_all(test_data)?;
        }
        
        // Test memory mapping
        let reader = MemoryMappedReader::new(&file_path)?;
        assert_eq!(reader.len(), test_data.len());
        assert_eq!(reader.as_slice(), test_data);
        
        // Test slicing
        let mut slice_reader = MemoryMappedReader::new(&file_path)?;
        slice_reader.seek(7); // Start at "mapped"
        assert_eq!(slice_reader.current_slice(6), b"mapped");
        
        // Test position tracking
        assert_eq!(slice_reader.position(), 7);
        slice_reader.advance(6);
        assert_eq!(slice_reader.position(), 13);
        
        Ok(())
    }
    
    #[test]
    fn test_split_file_into_chunks() -> io::Result<()> {
        // Create a temporary directory
        let dir = tempdir()?;
        let file_path = dir.path().join("chunk_test.txt");
        
        // Create a file with 1000 bytes
        {
            let mut file = File::create(&file_path)?;
            let data = vec![b'A'; 1000];
            file.write_all(&data)?;
        }
        
        // Split into chunks of 300 bytes
        let chunks = split_file_into_chunks(&file_path, 300)?;
        
        // Should have 4 chunks (300, 300, 300, 100)
        assert_eq!(chunks.len(), 4);
        assert_eq!(chunks[0], (0, 300));
        assert_eq!(chunks[1], (300, 600));
        assert_eq!(chunks[2], (600, 900));
        assert_eq!(chunks[3], (900, 1000));
        
        Ok(())
    }
}