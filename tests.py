import unittest
import solve

class SolveTests(unittest.TestCase):
    def test_nucleotide_count(self):
        seq = 'AGCTTTTCATTCTGACTGCAACGGGCAATATGTCTCTGTGTGGATTAAAAAAAGAGTGTCTGATAGCAGC'
        self.assertEqual(solve.nucleotide_count(seq), (20, 12, 17, 21))

if __name__ == '__main__':
    unittest.main()
