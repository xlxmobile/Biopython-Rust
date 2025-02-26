use pyo3::prelude::*;
use pyo3::wrap_pymodule;

pub mod engines;
pub mod modules;

/// 这个模块提供了Biopython的Rust实现
/// 专注于高性能序列操作
#[pymodule]
fn biopython_rust(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add_wrapped(wrap_pymodule!(engines::engine_module))?;
    m.add_wrapped(wrap_pymodule!(modules::module_module))?;
    
    m.add("__version__", env!("CARGO_PKG_VERSION"))?;
    
    Ok(())
}