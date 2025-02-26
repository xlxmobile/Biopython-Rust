use pyo3::prelude::*;

pub mod seq;
pub mod io;

/// 生物信息学功能模块
#[pymodule]
pub fn module_module(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add_submodule(seq::seq_module(_py)?)?;
    m.add_submodule(io::io_module(_py)?)?;
    
    Ok(())
}