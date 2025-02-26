use pyo3::prelude::*;

pub mod memory;
pub mod parallel;
pub mod io;
pub mod simd;

/// 核心优化引擎
#[pymodule]
pub fn core_module(_py: Python, m: &PyModule) -> PyResult<()> {
    memory::register_module(_py, m)?;
    // 其他模块将在后续阶段添加
    
    Ok(())
}