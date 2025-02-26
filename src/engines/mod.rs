use pyo3::prelude::*;

pub mod core;
pub mod compute;
pub mod storage;

/// Rust引擎模块包含高性能计算组件
#[pymodule]
pub fn engine_module(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add_submodule(core::core_module(_py)?)?;
    // 其他模块将在后续阶段添加
    
    Ok(())
}