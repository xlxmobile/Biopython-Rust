//! Parallel processing primitives
//!
//! This module provides a high-performance parallel processing framework
//! for biological sequence data, using work-stealing and adaptive chunking.

use std::sync::Once;
use rayon::{ThreadPool, ThreadPoolBuilder};
use std::sync::Mutex;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::Arc;
use parking_lot::RwLock;
use rayon::prelude::*;

// Initialize once
static INIT: Once = Once::new();

// Global thread pool for parallel operations
static mut GLOBAL_POOL: Option<ThreadPool> = None;

// Default chunk size for adaptive chunking
const DEFAULT_CHUNK_SIZE: usize = 1024 * 1024; // 1MB

// Default minimum chunks per thread
const MIN_CHUNKS_PER_THREAD: usize = 4;

/// Initialize the thread pool for parallel processing
pub fn initialize_thread_pool() {
    INIT.call_once(|| {
        // Create a thread pool with the number of CPUs
        let num_threads = default_num_threads();
        
        let pool = ThreadPoolBuilder::new()
            .num_threads(num_threads)
            .thread_name(|idx| format!("bioseq-worker-{}", idx))
            .build()
            .expect("Failed to create thread pool");
        
        // Store in global variable
        unsafe {
            GLOBAL_POOL = Some(pool);
        }
        
        log::info!("Initialized thread pool with {} threads", num_threads);
    });
}

/// Get the default number of threads to use
pub fn default_num_threads() -> usize {
    num_cpus::get()
}

/// Get a reference to the global thread pool
pub fn global_pool() -> &'static ThreadPool {
    unsafe {
        GLOBAL_POOL.as_ref().expect("Thread pool not initialized")
    }
}

/// Execute a closure in parallel with the global thread pool
pub fn execute<F, R>(f: F) -> R
where
    F: FnOnce(&ThreadPool) -> R + Send,
    R: Send,
{
    let pool = global_pool();
    f(pool)
}

/// Calculate optimal chunk size for parallel processing
pub fn calculate_chunk_size(total_size: usize, min_chunk_size: Option<usize>) -> usize {
    let num_threads = default_num_threads();
    let min_size = min_chunk_size.unwrap_or(1024); // Minimum 1KB
    
    // Calculate optimal chunk size
    let chunks_per_thread = MIN_CHUNKS_PER_THREAD;
    let total_chunks = num_threads * chunks_per_thread;
    
    // Ensure chunk size is at least min_size
    let chunk_size = (total_size / total_chunks).max(min_size);
    
    // Round to the nearest multiple of 1KB for better memory alignment
    let alignment = 1024;
    ((chunk_size + alignment - 1) / alignment) * alignment
}

/// Work-stealing scheduler for balanced parallel execution
pub struct WorkStealingScheduler<T> {
    /// Work items to process
    work_items: Mutex<Vec<T>>,
    /// Number of work items initially submitted
    total_items: usize,
    /// Number of completed work items
    completed: AtomicUsize,
}

impl<T: Send + 'static> WorkStealingScheduler<T> {
    /// Create a new work-stealing scheduler with the given work items
    pub fn new(work_items: Vec<T>) -> Self {
        let total_items = work_items.len();
        Self {
            work_items: Mutex::new(work_items),
            total_items,
            completed: AtomicUsize::new(0),
        }
    }
    
    /// Execute the work items in parallel using the given function
    pub fn execute<F>(&self, f: F)
    where
        F: Fn(T) + Send + Sync + Clone + 'static,
    {
        let pool = global_pool();
        
        pool.install(|| {
            rayon::scope(|s| {
                // Start workers equal to the number of threads
                for _ in 0..pool.current_num_threads() {
                    let f_clone = f.clone();
                    s.spawn(move |_| {
                        // Worker loop: grab work items and process them
                        loop {
                            // Try to get work
                            let work_item = {
                                let mut guard = self.work_items.lock().unwrap();
                                if guard.is_empty() {
                                    break;
                                }
                                guard.pop()
                            };
                            
                            if let Some(item) = work_item {
                                // Process the work item
                                f_clone(item);
                                
                                // Update completed count
                                self.completed.fetch_add(1, Ordering::SeqCst);
                            } else {
                                break;
                            }
                        }
                    });
                }
            });
        });
    }
    
    /// Get the progress of the execution (0.0-1.0)
    pub fn progress(&self) -> f64 {
        if self.total_items == 0 {
            return 1.0;
        }
        
        let completed = self.completed.load(Ordering::SeqCst);
        completed as f64 / self.total_items as f64
    }
    
    /// Check if all work items have been processed
    pub fn is_completed(&self) -> bool {
        let completed = self.completed.load(Ordering::SeqCst);
        completed == self.total_items
    }
}

/// Parallel sequence chunk processor
pub struct ParallelChunkProcessor<T> {
    /// Data divided into chunks
    chunks: Vec<T>,
    /// Results of processing
    results: Arc<RwLock<Vec<usize>>>,
}

impl<T: Send + Sync + 'static> ParallelChunkProcessor<T> {
    /// Create a new parallel chunk processor with the given chunks
    pub fn new(chunks: Vec<T>) -> Self {
        let num_chunks = chunks.len();
        Self {
            chunks,
            results: Arc::new(RwLock::new(Vec::with_capacity(num_chunks))),
        }
    }
    
    /// Process the chunks in parallel using the given function
    pub fn process<F, R>(&self, f: F) -> Vec<R>
    where
        F: Fn(&T) -> R + Send + Sync + Clone + 'static,
        R: Send + 'static,
    {
        let pool = global_pool();
        let results = Arc::new(Mutex::new(Vec::with_capacity(self.chunks.len())));
        
        pool.install(|| {
            let chunks_ref = &self.chunks;
            let results_ref = Arc::clone(&results);
            
            // Process chunks in parallel
            chunks_ref
                .par_iter()
                .map(|chunk| f(chunk))
                .collect_into_vec(&mut *results_ref.lock().unwrap());
        });
        
        // Return results
        let guard = results.lock().unwrap();
        guard.clone()
    }
}

/// Adaptive parallel execution based on workload
pub fn adaptive_parallel_execute<T, F, R>(items: Vec<T>, f: F) -> Vec<R>
where
    T: Send + Sync + 'static,
    F: Fn(&T) -> R + Send + Sync + Clone + 'static,
    R: Send + Default + 'static,
{
    let num_items = items.len();
    
    // For very small workloads, use sequential processing
    if num_items <= 8 {
        return items.iter().map(|item| f(item)).collect();
    }
    
    // For small to medium workloads, use rayon's par_iter directly
    if num_items <= 1000 {
        let pool = global_pool();
        let results = Arc::new(Mutex::new(Vec::with_capacity(num_items)));
        
        pool.install(|| {
            items.par_iter()
                .map(|item| f(item))
                .collect_into_vec(&mut *results.lock().unwrap());
        });
        
        let guard = results.lock().unwrap();
        return guard.clone();
    }
    
    // For large workloads, use the work-stealing scheduler
    let scheduler = WorkStealingScheduler::new(items);
    let results = Arc::new(Mutex::new(Vec::with_capacity(num_items)));
    let results_ref = Arc::clone(&results);
    
    scheduler.execute(move |item| {
        let result = f(&item);
        let mut guard = results_ref.lock().unwrap();
        guard.push(result);
    });
    
    let guard = results.lock().unwrap();
    guard.clone()
}

/// Chunk a slice into optimally sized chunks for parallel processing
pub fn chunk_slice<T>(slice: &[T], min_chunk_size: Option<usize>) -> Vec<&[T]> {
    let len = slice.len();
    let chunk_size = calculate_chunk_size(len, min_chunk_size);
    
    // Create chunks
    let mut chunks = Vec::new();
    let mut start = 0;
    
    while start < len {
        let end = (start + chunk_size).min(len);
        chunks.push(&slice[start..end]);
        start = end;
    }
    
    chunks
}

/// Split a task into parallel subtasks and join the results
pub fn parallel_split_join<S, R, FS, FR, J>(
    split_func: FS,
    process_func: FR,
    join_func: J,
) -> R
where
    S: Send + Sync + 'static,
    R: Send + 'static,
    FS: FnOnce() -> Vec<S> + Send + 'static,
    FR: Fn(S) -> R + Send + Sync + 'static,
    J: FnOnce(Vec<R>) -> R + Send + 'static,
{
    let pool = global_pool();
    
    pool.install(|| {
        // Split the task
        let subtasks = split_func();
        
        // Process subtasks in parallel
        let results: Vec<R> = subtasks
            .into_par_iter()
            .map(|subtask| process_func(subtask))
            .collect();
        
        // Join the results
        join_func(results)
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    
    #[test]
    fn test_chunk_size_calculation() {
        // Test with various total sizes
        let size_1mb = 1024 * 1024;
        let size_1gb = 1024 * 1024 * 1024;
        
        let chunk_1mb = calculate_chunk_size(size_1mb, None);
        let chunk_1gb = calculate_chunk_size(size_1gb, None);
        
        // Chunks should be properly sized
        assert!(chunk_1mb > 0);
        assert!(chunk_1gb > chunk_1mb);
        
        // Should be a multiple of 1KB
        assert_eq!(chunk_1mb % 1024, 0);
        assert_eq!(chunk_1gb % 1024, 0);
    }
    
    #[test]
    fn test_work_stealing_scheduler() {
        // Initialize thread pool
        initialize_thread_pool();
        
        // Create work items
        let items: Vec<usize> = (0..1000).collect();
        let scheduler = WorkStealingScheduler::new(items);
        
        // Run the work items
        let sum = Arc::new(AtomicUsize::new(0));
        let sum_ref = Arc::clone(&sum);
        
        scheduler.execute(move |item| {
            sum_ref.fetch_add(item, Ordering::SeqCst);
        });
        
        // Check that all work was completed
        assert!(scheduler.is_completed());
        
        // Check the result (sum of numbers 0-999)
        let expected_sum = (0..1000).sum();
        assert_eq!(sum.load(Ordering::SeqCst), expected_sum);
    }
    
    #[test]
    fn test_parallel_chunk_processing() {
        // Initialize thread pool
        initialize_thread_pool();
        
        // Create data chunks
        let data: Vec<Vec<usize>> = vec![
            vec![1, 2, 3],
            vec![4, 5, 6],
            vec![7, 8, 9],
            vec![10, 11, 12],
        ];
        
        let processor = ParallelChunkProcessor::new(data);
        
        // Process chunks (sum each chunk)
        let results = processor.process(|chunk| chunk.iter().sum::<usize>());
        
        // Check results
        assert_eq!(results, vec![6, 15, 24, 33]);
    }
}