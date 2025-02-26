use pyo3::prelude::*;
use bitvec::prelude::*;
use std::borrow::Cow;
use std::fs::File;
use std::io::{self, Read, Seek};
use std::path::Path;
use std::sync::Arc;
use memmap2::{Mmap, MmapOptions};
use thiserror::Error;

/// 内存模块错误类型
#[derive(Error, Debug)]
pub enum MemoryError {
    #[error("I/O error: {0}")]
    IoError(#[from] io::Error),
    
    #[error("Invalid sequence data: {0}")]
    InvalidSequence(String),
    
    #[error("Unsupported operation: {0}")]
    UnsupportedOperation(String),
    
    #[error("Index out of bounds")]
    IndexOutOfBounds,
}

/// 将内存错误转换为Python异常
impl From<MemoryError> for PyErr {
    fn from(err: MemoryError) -> PyErr {
        match err {
            MemoryError::IoError(e) => {
                pyo3::exceptions::PyIOError::new_err(e.to_string())
            }
            MemoryError::InvalidSequence(msg) => {
                pyo3::exceptions::PyValueError::new_err(msg)
            }
            MemoryError::UnsupportedOperation(msg) => {
                pyo3::exceptions::PyNotImplementedError::new_err(msg)
            }
            MemoryError::IndexOutOfBounds => {
                pyo3::exceptions::PyIndexError::new_err("Index out of bounds")
            }
        }
    }
}

// DNA/RNA 碱基的 2-bit 编码
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[repr(u8)]
pub enum NucleotideCode {
    A = 0b00,
    C = 0b01,
    G = 0b10,
    T = 0b11,
}

impl NucleotideCode {
    pub fn from_char(c: char) -> Result<Self, MemoryError> {
        match c.to_ascii_uppercase() {
            'A' => Ok(NucleotideCode::A),
            'C' => Ok(NucleotideCode::C),
            'G' => Ok(NucleotideCode::G),
            'T' | 'U' => Ok(NucleotideCode::T),
            _ => Err(MemoryError::InvalidSequence(format!("Invalid nucleotide: {}", c)))
        }
    }
    
    pub fn to_char(self, is_rna: bool) -> char {
        match self {
            NucleotideCode::A => 'A',
            NucleotideCode::C => 'C',
            NucleotideCode::G => 'G',
            NucleotideCode::T => if is_rna { 'U' } else { 'T' },
        }
    }
}

// 蛋白质字母表编码 (5-bit)
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[repr(u8)]
pub enum AminoAcidCode {
    A = 0,   // Alanine
    R = 1,   // Arginine
    N = 2,   // Asparagine
    D = 3,   // Aspartic acid
    C = 4,   // Cysteine
    Q = 5,   // Glutamine
    E = 6,   // Glutamic acid
    G = 7,   // Glycine
    H = 8,   // Histidine
    I = 9,   // Isoleucine
    L = 10,  // Leucine
    K = 11,  // Lysine
    M = 12,  // Methionine
    F = 13,  // Phenylalanine
    P = 14,  // Proline
    S = 15,  // Serine
    T = 16,  // Threonine
    W = 17,  // Tryptophan
    Y = 18,  // Tyrosine
    V = 19,  // Valine
    B = 20,  // Asparagine or aspartic acid
    Z = 21,  // Glutamine or glutamic acid
    X = 22,  // Any amino acid
    Stop = 23, // Stop codon '*'
}

impl AminoAcidCode {
    pub fn from_char(c: char) -> Result<Self, MemoryError> {
        match c.to_ascii_uppercase() {
            'A' => Ok(AminoAcidCode::A),
            'R' => Ok(AminoAcidCode::R),
            'N' => Ok(AminoAcidCode::N),
            'D' => Ok(AminoAcidCode::D),
            'C' => Ok(AminoAcidCode::C),
            'Q' => Ok(AminoAcidCode::Q),
            'E' => Ok(AminoAcidCode::E),
            'G' => Ok(AminoAcidCode::G),
            'H' => Ok(AminoAcidCode::H),
            'I' => Ok(AminoAcidCode::I),
            'L' => Ok(AminoAcidCode::L),
            'K' => Ok(AminoAcidCode::K),
            'M' => Ok(AminoAcidCode::M),
            'F' => Ok(AminoAcidCode::F),
            'P' => Ok(AminoAcidCode::P),
            'S' => Ok(AminoAcidCode::S),
            'T' => Ok(AminoAcidCode::T),
            'W' => Ok(AminoAcidCode::W),
            'Y' => Ok(AminoAcidCode::Y),
            'V' => Ok(AminoAcidCode::V),
            'B' => Ok(AminoAcidCode::B),
            'Z' => Ok(AminoAcidCode::Z),
            'X' => Ok(AminoAcidCode::X),
            '*' => Ok(AminoAcidCode::Stop),
            _ => Err(MemoryError::InvalidSequence(format!("Invalid amino acid: {}", c)))
        }
    }
    
    pub fn to_char(self) -> char {
        match self {
            AminoAcidCode::A => 'A',
            AminoAcidCode::R => 'R',
            AminoAcidCode::N => 'N',
            AminoAcidCode::D => 'D',
            AminoAcidCode::C => 'C',
            AminoAcidCode::Q => 'Q',
            AminoAcidCode::E => 'E',
            AminoAcidCode::G => 'G',
            AminoAcidCode::H => 'H',
            AminoAcidCode::I => 'I',
            AminoAcidCode::L => 'L',
            AminoAcidCode::K => 'K',
            AminoAcidCode::M => 'M',
            AminoAcidCode::F => 'F',
            AminoAcidCode::P => 'P',
            AminoAcidCode::S => 'S',
            AminoAcidCode::T => 'T',
            AminoAcidCode::W => 'W',
            AminoAcidCode::Y => 'Y',
            AminoAcidCode::V => 'V',
            AminoAcidCode::B => 'B',
            AminoAcidCode::Z => 'Z',
            AminoAcidCode::X => 'X',
            AminoAcidCode::Stop => '*',
        }
    }
}

/// 序列类型枚举
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SequenceType {
    DNA,
    RNA,
    Protein,
    Generic,
}

/// 序列存储特性，为不同的存储策略定义公共接口
pub trait SequenceStorage {
    fn len(&self) -> usize;
    fn is_empty(&self) -> bool {
        self.len() == 0
    }
    fn get_char(&self, index: usize) -> Result<char, MemoryError>;
    fn to_string(&self) -> String;
    fn slice(&self, start: usize, end: usize) -> Result<Box<dyn SequenceStorage>, MemoryError>;
    fn get_type(&self) -> SequenceType;
}

/// 紧凑型核酸序列存储 (2-bit per base)
#[pyclass]
#[derive(Debug, Clone)]
pub struct CompactDnaStorage {
    data: BitVec<u8, Msb0>,
    length: usize,
    seq_type: SequenceType,
}

impl CompactDnaStorage {
    pub fn new(sequence: &str, seq_type: SequenceType) -> Result<Self, MemoryError> {
        if seq_type != SequenceType::DNA && seq_type != SequenceType::RNA {
            return Err(MemoryError::InvalidSequence(
                "CompactDnaStorage only supports DNA or RNA sequences".to_string()));
        }
        
        let length = sequence.len();
        // 每个碱基需要2位，需要预先分配足够的空间
        let mut data = BitVec::<u8, Msb0>::with_capacity(length * 2);
        
        for c in sequence.chars() {
            let code = NucleotideCode::from_char(c)?;
            // 添加2位，表示一个碱基
            data.push(code as u8 & 0b10 != 0);
            data.push(code as u8 & 0b01 != 0);
        }
        
        Ok(CompactDnaStorage {
            data,
            length,
            seq_type,
        })
    }
}

impl SequenceStorage for CompactDnaStorage {
    fn len(&self) -> usize {
        self.length
    }
    
    fn get_char(&self, index: usize) -> Result<char, MemoryError> {
        if index >= self.length {
            return Err(MemoryError::IndexOutOfBounds);
        }
        
        let bit_index = index * 2;
        let msb = self.data[bit_index];
        let lsb = self.data[bit_index + 1];
        
        let code_value = (msb as u8) << 1 | (lsb as u8);
        let code = match code_value {
            0b00 => NucleotideCode::A,
            0b01 => NucleotideCode::C,
            0b10 => NucleotideCode::G,
            0b11 => NucleotideCode::T,
            _ => unreachable!(),
        };
        
        Ok(code.to_char(self.seq_type == SequenceType::RNA))
    }
    
    fn to_string(&self) -> String {
        let mut result = String::with_capacity(self.length);
        for i in 0..self.length {
            if let Ok(c) = self.get_char(i) {
                result.push(c);
            }
        }
        result
    }
    
    fn slice(&self, start: usize, end: usize) -> Result<Box<dyn SequenceStorage>, MemoryError> {
        if start >= self.length || end > self.length || start > end {
            return Err(MemoryError::IndexOutOfBounds);
        }
        
        let slice_length = end - start;
        let mut result = BitVec::<u8, Msb0>::with_capacity(slice_length * 2);
        
        for i in start..end {
            let bit_index = i * 2;
            result.push(self.data[bit_index]);
            result.push(self.data[bit_index + 1]);
        }
        
        Ok(Box::new(CompactDnaStorage {
            data: result,
            length: slice_length,
            seq_type: self.seq_type,
        }))
    }
    
    fn get_type(&self) -> SequenceType {
        self.seq_type
    }
}

/// 紧凑型蛋白质序列存储 (5-bit per amino acid)
#[pyclass]
#[derive(Debug, Clone)]
pub struct CompactProteinStorage {
    // 每个氨基酸用5位表示，打包到u8数组中
    data: Vec<u8>,
    // 存储序列的实际长度（氨基酸数量）
    length: usize,
}

impl CompactProteinStorage {
    pub fn new(sequence: &str) -> Result<Self, MemoryError> {
        let length = sequence.len();
        
        // 计算需要的字节数：每个氨基酸5位，8位一个字节
        // 需要向上取整以确保足够空间
        let bytes_needed = (length * 5 + 7) / 8;
        let mut data = vec![0u8; bytes_needed];
        
        for (i, c) in sequence.chars().enumerate() {
            let code = AminoAcidCode::from_char(c)? as u8;
            
            // 计算氨基酸在bit序列中的位置
            let bit_pos = i * 5;
            let byte_idx = bit_pos / 8;
            let bit_offset = bit_pos % 8;
            
            // 处理跨字节边界情况
            if bit_offset <= 3 {
                // 氨基酸编码完全在一个字节内
                data[byte_idx] |= code << (3 - bit_offset);
            } else {
                // 氨基酸编码跨越两个字节
                let bits_in_first = 8 - bit_offset;
                let bits_in_second = 5 - bits_in_first;
                
                // 第一个字节
                data[byte_idx] |= code >> bits_in_second;
                
                // 第二个字节（如果在范围内）
                if byte_idx + 1 < bytes_needed {
                    data[byte_idx + 1] |= code << (8 - bits_in_second);
                }
            }
        }
        
        Ok(CompactProteinStorage {
            data,
            length,
        })
    }
    
    fn get_code(&self, index: usize) -> Result<AminoAcidCode, MemoryError> {
        if index >= self.length {
            return Err(MemoryError::IndexOutOfBounds);
        }
        
        let bit_pos = index * 5;
        let byte_idx = bit_pos / 8;
        let bit_offset = bit_pos % 8;
        
        let mut code;
        
        if bit_offset <= 3 {
            // 氨基酸完全在一个字节内
            code = (self.data[byte_idx] >> (3 - bit_offset)) & 0x1F;
        } else {
            // 氨基酸跨越两个字节
            let bits_in_first = 8 - bit_offset;
            let bits_in_second = 5 - bits_in_first;
            
            code = (self.data[byte_idx] & ((1 << bits_in_first) - 1)) << bits_in_second;
            
            // 如果第二个字节在范围内
            if byte_idx + 1 < self.data.len() {
                code |= self.data[byte_idx + 1] >> (8 - bits_in_second);
            }
        }
        
        match code {
            0 => Ok(AminoAcidCode::A),
            1 => Ok(AminoAcidCode::R),
            2 => Ok(AminoAcidCode::N),
            3 => Ok(AminoAcidCode::D),
            4 => Ok(AminoAcidCode::C),
            5 => Ok(AminoAcidCode::Q),
            6 => Ok(AminoAcidCode::E),
            7 => Ok(AminoAcidCode::G),
            8 => Ok(AminoAcidCode::H),
            9 => Ok(AminoAcidCode::I),
            10 => Ok(AminoAcidCode::L),
            11 => Ok(AminoAcidCode::K),
            12 => Ok(AminoAcidCode::M),
            13 => Ok(AminoAcidCode::F),
            14 => Ok(AminoAcidCode::P),
            15 => Ok(AminoAcidCode::S),
            16 => Ok(AminoAcidCode::T),
            17 => Ok(AminoAcidCode::W),
            18 => Ok(AminoAcidCode::Y),
            19 => Ok(AminoAcidCode::V),
            20 => Ok(AminoAcidCode::B),
            21 => Ok(AminoAcidCode::Z),
            22 => Ok(AminoAcidCode::X),
            23 => Ok(AminoAcidCode::Stop),
            _ => Err(MemoryError::InvalidSequence(format!("Invalid amino acid code: {}", code))),
        }
    }
}

impl SequenceStorage for CompactProteinStorage {
    fn len(&self) -> usize {
        self.length
    }
    
    fn get_char(&self, index: usize) -> Result<char, MemoryError> {
        let code = self.get_code(index)?;
        Ok(code.to_char())
    }
    
    fn to_string(&self) -> String {
        let mut result = String::with_capacity(self.length);
        for i in 0..self.length {
            if let Ok(c) = self.get_char(i) {
                result.push(c);
            }
        }
        result
    }
    
    fn slice(&self, start: usize, end: usize) -> Result<Box<dyn SequenceStorage>, MemoryError> {
        if start >= self.length || end > self.length || start > end {
            return Err(MemoryError::IndexOutOfBounds);
        }
        
        // 创建一个新的序列字符串，然后从中构建压缩存储
        let mut slice_str = String::with_capacity(end - start);
        for i in start..end {
            slice_str.push(self.get_char(i)?);
        }
        
        Ok(Box::new(CompactProteinStorage::new(&slice_str)?))
    }
    
    fn get_type(&self) -> SequenceType {
        SequenceType::Protein
    }
}

/// 内存映射序列存储，用于处理大型序列文件
#[pyclass]
#[derive(Debug)]
pub struct MmapSequenceStorage {
    mmap: Arc<Mmap>,
    offset: usize,
    length: usize,
    seq_type: SequenceType,
}

impl MmapSequenceStorage {
    pub fn from_file<P: AsRef<Path>>(path: P, seq_type: SequenceType) -> Result<Self, MemoryError> {
        let file = File::open(path)?;
        let mmap = unsafe { MmapOptions::new().map(&file)? };
        
        Ok(MmapSequenceStorage {
            mmap: Arc::new(mmap),
            offset: 0,
            length: mmap.len(),
            seq_type,
        })
    }
    
    pub fn new(mmap: Arc<Mmap>, offset: usize, length: usize, seq_type: SequenceType) -> Self {
        MmapSequenceStorage {
            mmap,
            offset,
            length,
            seq_type,
        }
    }
}

impl SequenceStorage for MmapSequenceStorage {
    fn len(&self) -> usize {
        self.length
    }
    
    fn get_char(&self, index: usize) -> Result<char, MemoryError> {
        if index >= self.length {
            return Err(MemoryError::IndexOutOfBounds);
        }
        
        let byte = self.mmap[self.offset + index] as char;
        
        // 验证读取的数据是有效的序列字符
        match self.seq_type {
            SequenceType::DNA => {
                match byte.to_ascii_uppercase() {
                    'A' | 'C' | 'G' | 'T' | 'N' => Ok(byte),
                    _ => Err(MemoryError::InvalidSequence(format!("Invalid DNA base: {}", byte))),
                }
            },
            SequenceType::RNA => {
                match byte.to_ascii_uppercase() {
                    'A' | 'C' | 'G' | 'U' | 'N' => Ok(byte),
                    _ => Err(MemoryError::InvalidSequence(format!("Invalid RNA base: {}", byte))),
                }
            },
            SequenceType::Protein => {
                // 简单验证是否是有效的氨基酸字符
                let valid_aa = "ACDEFGHIKLMNPQRSTVWYBZX*";
                if valid_aa.contains(byte.to_ascii_uppercase()) {
                    Ok(byte)
                } else {
                    Err(MemoryError::InvalidSequence(format!("Invalid amino acid: {}", byte)))
                }
            },
            SequenceType::Generic => Ok(byte),
        }
    }
    
    fn to_string(&self) -> String {
        let slice = &self.mmap[self.offset..self.offset + self.length];
        
        // 假设数据是有效的UTF-8数据
        match std::str::from_utf8(slice) {
            Ok(s) => s.to_string(),
            Err(_) => {
                // 如果不是有效的UTF-8，则一个字符一个字符地构建
                let mut result = String::with_capacity(self.length);
                for i in 0..self.length {
                    if let Ok(c) = self.get_char(i) {
                        result.push(c);
                    } else {
                        result.push('?'); // 用问号代替无效字符
                    }
                }
                result
            }
        }
    }
    
    fn slice(&self, start: usize, end: usize) -> Result<Box<dyn SequenceStorage>, MemoryError> {
        if start >= self.length || end > self.length || start > end {
            return Err(MemoryError::IndexOutOfBounds);
        }
        
        Ok(Box::new(MmapSequenceStorage::new(
            self.mmap.clone(),
            self.offset + start,
            end - start,
            self.seq_type
        )))
    }
    
    fn get_type(&self) -> SequenceType {
        self.seq_type
    }
}

/// 标准字符串序列存储
#[pyclass]
#[derive(Debug, Clone)]
pub struct StringSequenceStorage {
    data: String,
    seq_type: SequenceType,
}

impl StringSequenceStorage {
    pub fn new(sequence: &str, seq_type: SequenceType) -> Self {
        StringSequenceStorage {
            data: sequence.to_string(),
            seq_type,
        }
    }
}

impl SequenceStorage for StringSequenceStorage {
    fn len(&self) -> usize {
        self.data.len()
    }
    
    fn get_char(&self, index: usize) -> Result<char, MemoryError> {
        self.data.chars().nth(index).ok_or(MemoryError::IndexOutOfBounds)
    }
    
    fn to_string(&self) -> String {
        self.data.clone()
    }
    
    fn slice(&self, start: usize, end: usize) -> Result<Box<dyn SequenceStorage>, MemoryError> {
        if start >= self.data.len() || end > self.data.len() || start > end {
            return Err(MemoryError::IndexOutOfBounds);
        }
        
        Ok(Box::new(StringSequenceStorage {
            data: self.data[start..end].to_string(),
            seq_type: self.seq_type,
        }))
    }
    
    fn get_type(&self) -> SequenceType {
        self.seq_type
    }
}

/// 存储工厂函数，用于根据序列类型和长度选择最优的存储方式
#[pyfunction]
pub fn create_optimal_storage(sequence: &str, seq_type: &str) -> Result<PyObject, MemoryError> {
    let seq_type = match seq_type.to_lowercase().as_str() {
        "dna" => SequenceType::DNA,
        "rna" => SequenceType::RNA,
        "protein" => SequenceType::Protein,
        _ => SequenceType::Generic,
    };
    
    let length = sequence.len();
    
    Python::with_gil(|py| {
        // 根据序列类型和长度选择最优存储
        let storage: Box<dyn SequenceStorage> = match seq_type {
            SequenceType::DNA | SequenceType::RNA if length > 1000 => {
                // 对于较长的核酸序列，使用紧凑存储
                Box::new(CompactDnaStorage::new(sequence, seq_type)?)
            },
            SequenceType::Protein if length > 1000 => {
                // 对于较长的蛋白质序列，使用紧凑存储
                Box::new(CompactProteinStorage::new(sequence)?)
            },
            _ => {
                // 对于较短的序列，使用标准字符串存储
                Box::new(StringSequenceStorage::new(sequence, seq_type))
            }
        };
        
        // 创建PyObject
        let py_storage = PySequenceStorage::new(storage);
        Ok(Py::new(py, py_storage)?.into_py(py))
    })
}

#[pyfunction]
pub fn mmap_sequence_file(path: &str, seq_type: &str) -> Result<PyObject, MemoryError> {
    let seq_type = match seq_type.to_lowercase().as_str() {
        "dna" => SequenceType::DNA,
        "rna" => SequenceType::RNA,
        "protein" => SequenceType::Protein,
        _ => SequenceType::Generic,
    };
    
    let storage = MmapSequenceStorage::from_file(path, seq_type)?;
    
    Python::with_gil(|py| {
        let py_storage = PySequenceStorage::new(Box::new(storage));
        Ok(Py::new(py, py_storage)?.into_py(py))
    })
}

/// Python绑定的序列存储类
#[pyclass]
#[derive(Debug)]
pub struct PySequenceStorage {
    storage: Box<dyn SequenceStorage>,
}

#[pymethods]
impl PySequenceStorage {
    #[new]
    fn new(storage: Box<dyn SequenceStorage>) -> Self {
        PySequenceStorage { storage }
    }
    
    #[getter]
    fn get_length(&self) -> usize {
        self.storage.len()
    }
    
    fn get(&self, index: usize) -> PyResult<String> {
        match self.storage.get_char(index) {
            Ok(c) => Ok(c.to_string()),
            Err(e) => Err(e.into()),
        }
    }
    
    fn slice(&self, start: usize, end: Option<usize>) -> PyResult<Self> {
        let end = end.unwrap_or(self.storage.len());
        match self.storage.slice(start, end) {
            Ok(slice) => Ok(PySequenceStorage { storage: slice }),
            Err(e) => Err(e.into()),
        }
    }
    
    fn to_string(&self) -> String {
        self.storage.to_string()
    }
    
    fn __str__(&self) -> String {
        self.storage.to_string()
    }
    
    fn __len__(&self) -> usize {
        self.storage.len()
    }
    
    fn __getitem__(&self, idx: isize) -> PyResult<String> {
        let len = self.storage.len() as isize;
        let idx = if idx < 0 { len + idx } else { idx };
        
        if idx < 0 || idx >= len {
            return Err(PyErr::new::<pyo3::exceptions::PyIndexError, _>("Index out of bounds"));
        }
        
        self.get(idx as usize)
    }
}

/// 为Python模块注册函数和类
pub fn register_module(_py: Python<'_>, m: &PyModule) -> PyResult<()> {
    let submod = PyModule::new(_py, "memory")?;
    
    submod.add_class::<PySequenceStorage>()?;
    submod.add_function(wrap_pyfunction!(create_optimal_storage, submod)?)?;
    submod.add_function(wrap_pyfunction!(mmap_sequence_file, submod)?)?;
    
    m.add_submodule(submod)?;
    
    Ok(())
}