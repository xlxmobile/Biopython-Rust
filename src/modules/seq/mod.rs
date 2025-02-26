use pyo3::prelude::*;

pub mod sequence;
pub mod alphabet;

/// 序列处理模块
#[pymodule]
pub fn seq_module(_py: Python, m: &PyModule) -> PyResult<()> {
    sequence::register_module(_py, m)?;
    alphabet::register_module(_py, m)?;
    
    Ok(())
}