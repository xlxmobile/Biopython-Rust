import unittest
from biopython.seq import Sequence, DNASequence, RNASequence, ProteinSequence

class TestSequence(unittest.TestCase):
    def test_basic_sequence(self):
        """测试基本序列操作"""
        seq = Sequence("ACGTACGT", "dna")
        self.assertEqual(len(seq), 8)
        self.assertEqual(str(seq), "ACGTACGT")
        self.assertEqual(seq[0], "A")
        self.assertEqual(seq[4], "A")
        
    def test_dna_sequence(self):
        """测试DNA序列创建"""
        dna = DNASequence("ACGTACGT")
        self.assertEqual(len(dna), 8)
        self.assertEqual(str(dna), "ACGTACGT")
        
    def test_rna_sequence(self):
        """测试RNA序列创建"""
        rna = RNASequence("ACGUACGU")
        self.assertEqual(len(rna), 8)
        self.assertEqual(str(rna), "ACGUACGU")
        
    def test_protein_sequence(self):
        """测试蛋白质序列创建"""
        protein = ProteinSequence("MKVWDIFQEY")
        self.assertEqual(len(protein), 10)
        self.assertEqual(str(protein), "MKVWDIFQEY")
        
    def test_slice(self):
        """测试序列切片"""
        seq = Sequence("ACGTACGT", "dna")
        slice1 = seq.slice(2, 6)
        self.assertEqual(str(slice1), "GTAC")
        slice2 = seq[2:6]
        self.assertEqual(str(slice2), "GTAC")

if __name__ == "__main__":
    unittest.main()