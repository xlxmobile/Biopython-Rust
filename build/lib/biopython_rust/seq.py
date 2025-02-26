from biopython_rust.engines.core.memory import PySequenceStorage, create_optimal_storage, mmap_sequence_file

class Sequence:
    """高性能序列类，使用Rust实现的内存优化存储"""
    
    def __init__(self, sequence, seq_type="generic", id="", name="", description=""):
        """初始化序列对象
        
        参数:
            sequence (str): 序列字符串
            seq_type (str): 序列类型 ('dna', 'rna', 'protein', 'generic')
            id (str): 序列ID
            name (str): 序列名称
            description (str): 序列描述
        """
        self.storage = create_optimal_storage(sequence, seq_type)
        self.id = id
        self.name = name
        self.description = description
        self.seq_type = seq_type
        
    def __len__(self):
        """返回序列长度"""
        return self.storage.get_length()
    
    def __getitem__(self, index):
        """通过索引访问序列"""
        if isinstance(index, slice):
            start = index.start or 0
            stop = index.stop or len(self)
            return self.slice(start, stop)
        return self.storage[index]
    
    def __str__(self):
        """返回序列字符串表示"""
        return self.storage.to_string()
    
    def slice(self, start, end=None):
        """获取序列片段"""
        if end is None:
            end = len(self)
            
        new_storage = self.storage.slice(start, end)
        result = Sequence("", self.seq_type, self.id, self.name, self.description)
        result.storage = new_storage
        return result
    
    @classmethod
    def from_file(cls, path, seq_type="dna", id="", name="", description=""):
        """从文件加载序列，使用内存映射优化大文件处理"""
        storage = mmap_sequence_file(path, seq_type)
        seq = cls("", seq_type, id, name, description)
        seq.storage = storage
        return seq

class DNASequence(Sequence):
    """DNA序列特化类"""
    
    def __init__(self, sequence, id="", name="", description=""):
        super().__init__(sequence, "dna", id, name, description)

class RNASequence(Sequence):
    """RNA序列特化类"""
    
    def __init__(self, sequence, id="", name="", description=""):
        super().__init__(sequence, "rna", id, name, description)

class ProteinSequence(Sequence):
    """蛋白质序列特化类"""
    
    def __init__(self, sequence, id="", name="", description=""):
        super().__init__(sequence, "protein", id, name, description)