use pyo3::prelude::*;

pub mod fasta;

/// 文件I/O模块
#[pymodule]
pub fn io_module(_py: Python, m: &PyModule) -> PyResult<()> {
    // 这将在下一个开发阶段实现
    
    Ok(())
}