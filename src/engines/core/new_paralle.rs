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
        
        // For very small number of subtasks, use simpler approach
        if subtasks.len() <= 4 {
            // Process subtasks in parallel
            let results: Vec<R> = subtasks
                .into_par_iter()
                .map(|subtask| process_func(subtask))
                .collect();
            
            // Join the results
            join_func(results)
        } else {
            // Use our advanced work-stealing scheduler for better load balancing
            
            // Prioritize tasks - first and last tasks often contain edge cases
            let mut prioritized_tasks = Vec::with_capacity(subtasks.len());
            
            for (i, subtask) in subtasks.into_iter().enumerate() {
                let priority = if i == 0 || i == prioritized_tasks.len() - 1 {
                    // First and last tasks get high priority
                    Priority::High
                } else {
                    // Other tasks get normal priority
                    Priority::Normal
                };
                
                prioritized_tasks.push(PrioritizedWork::new(subtask, priority));
            }
            
            // Create scheduler
            let scheduler = WorkStealingScheduler::with_prioritized_items(prioritized_tasks);
            
            // Process tasks and collect results
            let results = scheduler.execute_with_results(|subtask| process_func(subtask));
            
            // Join the results
            join_func(results)
        }
    })
}