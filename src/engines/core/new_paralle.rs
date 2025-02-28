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
    
    // For large workloads, use our advanced work-stealing scheduler with adaptive chunking
    
    // Create chunks with adaptive sizing
    let config = ChunkingConfig::new()
        .with_strategy(ChunkingStrategy::FullyAdaptive);
    
    // Create data chunks
    let chunks: Vec<_> = items
        .chunks(config.calculate_chunk_size(num_items))
        .enumerate()
        .map(|(i, chunk_items)| {
            let is_first = i == 0;
            let is_last = (i + 1) * config.calculate_chunk_size(num_items) >= num_items;
            DataChunk::new(chunk_items, i * chunk_items.len(), (i + 1) * chunk_items.len(), is_first, is_last)
        })
        .collect();
    
    // Process using the parallel chunk processor
    let processor = ParallelChunkProcessor::new(chunks);
    
    // Use the process_with_scheduler method to get the benefits of our prioritized scheduler
    processor.process_with_scheduler(|chunk| {
        chunk.data.iter().map(|item| f(item)).collect::<Vec<R>>()
    })
    .into_iter()
    .flatten()
    .collect()
}