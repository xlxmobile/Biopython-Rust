use pyo3::prelude::*;

/// 序列字母表类型
#[pyclass]
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AlphabetType {
    #[pyo3(name = "DNA")]
    DNA,
    #[pyo3(name = "RNA")]
    RNA,
    #[pyo3(name = "PROTEIN")]
    Protein,
    #[pyo3(name = "GENERIC")]
    Generic,
}

#[pyfunction]
fn get_alphabet_letters(alphabet_type: AlphabetType) -> Vec<char> {
    match alphabet_type {
        AlphabetType::DNA => vec!['A', 'C', 'G', 'T', 'N'],
        AlphabetType::RNA => vec!['A', 'C', 'G', 'U', 'N'],
        AlphabetType::Protein => vec![
            'A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L',
            'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y',
            'B', 'Z', 'X', '*'
        ],
        AlphabetType::Generic => vec!['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 
                                      'N', 'O', 'P', 'Q', 'R', 'S', 'T', 'U', 'V', 'W', 'X', 'Y', 'Z'],
    }
}

#[pyfunction]
fn is_valid_sequence(sequence: &str, alphabet_type: AlphabetType) -> bool {
    let letters = get_alphabet_letters(alphabet_type);
    
    for c in sequence.chars() {
        let upper_c = c.to_ascii_uppercase();
        if !letters.contains(&upper_c) {
            return false;
        }
    }
    
    true
}

pub fn register_module(_py: Python<'_>, m: &PyModule) -> PyResult<()> {
    let submod = PyModule::new(_py, "alphabet")?;
    
    submod.add_class::<AlphabetType>()?;
    submod.add_function(wrap_pyfunction!(get_alphabet_letters, submod)?)?;
    submod.add_function(wrap_pyfunction!(is_valid_sequence, submod)?)?;
    
    m.add_submodule(submod)?;
    
    Ok(())
}