"""
Test that input provided by Rosalind gives correct answer.

Input and output are either taken from the example problem on the
Rosalind problem page (to make the tests run as fast as possible) or
from the actual question generated by Rosalind to pass to the next level.
"""

import solve
import unittest


class SolveTests(unittest.TestCase):

    # p.1
    def test_nucleotide_count(self):
        seq = 'AGCTTTTCATTCTGACTGCAACGGGCAATATGTCTCTGTGTGGATTAAAAAAAGAGTGTCTGATAGCAGC'
        self.assertEqual(solve.nucleotide_count(seq), (20, 12, 17, 21))

    # p.2
    def test_dna_to_rna(self):
        dna = 'TCTCACGAGCCCGCGATATATGGCAATCACCCCTGGCGTCACGGGGGGCAGTACATTCGCTGAGACGACCGTAGAGCTATACTAAAAGTCGTAGCTTTGTTTAGATCCTTAGAAACCACGTCACTGAAGCCGCCTTCGGCGCGTGGACGTCGGTCAAGTCATTATATGGATAGTGCGTTGGTGTAAAACGCCGCAACACTAATTAGTATGATAGAGTTAACACAGACGCACCGACCTCTTTATCCAAGGCAGGCATTCAGACGAAAGCGTGCCAGAAAGCGGAAGTCTCGGAGCGCCTCGTTCTGCATTTTTTGGCGAACGTCGAGCTCCGGAGTCCGGGGTTTCATCTACGGAGCCTGTCTCATGGGTGTAAGCCTGTTCAGTAATCTCCGACTGGCGACCTCCCGTTTAAACGGCACTTATGAATGTTCCTGCTCATTAGCCTTTGCATAAGCTTAATTTACGCACGAGTACACCTTGTAAGCATGTCTGCACTTCCCAATGTCCCCGCTCCGGATTGGAAGGGAGGTTCTATCTTTACGAGACTGGTGTGCCGCATTTCCTGTCCAGCTTGAGGGACCTCCGCGGGGGCGCCCCTATTTAATGTGAGGACTTATAAATCGCAGGGACTCACCAAGTACCATTACAAACACTACACGCTCGGCCATATTGCGTCGCCGGACACTACCGGCAGACCGAAATTCCTTTGTTGCCGCACCAGGGTTCTGAGGTTTATATACCGTTCGCCCTTGAACAACCCGGGACCAATGGGTAGTGACGATGAGGAGTTGCGTAAGGGCAAGTTATGCTTTTCTGAAGTATAAACTCAACTAGGTCGACAGATGCAAGGTGTGCGAAAAAAACCAACTAGTTATTTGTCCCCTAGAGTGGTGTTGGAAAGCCGAGGGGCACATATTCATGTAAGGGGGGACCGCTTACTAGCCCGATTGGTCCCTCACCCGGTGTAAAACGGCT'
        rna = 'UCUCACGAGCCCGCGAUAUAUGGCAAUCACCCCUGGCGUCACGGGGGGCAGUACAUUCGCUGAGACGACCGUAGAGCUAUACUAAAAGUCGUAGCUUUGUUUAGAUCCUUAGAAACCACGUCACUGAAGCCGCCUUCGGCGCGUGGACGUCGGUCAAGUCAUUAUAUGGAUAGUGCGUUGGUGUAAAACGCCGCAACACUAAUUAGUAUGAUAGAGUUAACACAGACGCACCGACCUCUUUAUCCAAGGCAGGCAUUCAGACGAAAGCGUGCCAGAAAGCGGAAGUCUCGGAGCGCCUCGUUCUGCAUUUUUUGGCGAACGUCGAGCUCCGGAGUCCGGGGUUUCAUCUACGGAGCCUGUCUCAUGGGUGUAAGCCUGUUCAGUAAUCUCCGACUGGCGACCUCCCGUUUAAACGGCACUUAUGAAUGUUCCUGCUCAUUAGCCUUUGCAUAAGCUUAAUUUACGCACGAGUACACCUUGUAAGCAUGUCUGCACUUCCCAAUGUCCCCGCUCCGGAUUGGAAGGGAGGUUCUAUCUUUACGAGACUGGUGUGCCGCAUUUCCUGUCCAGCUUGAGGGACCUCCGCGGGGGCGCCCCUAUUUAAUGUGAGGACUUAUAAAUCGCAGGGACUCACCAAGUACCAUUACAAACACUACACGCUCGGCCAUAUUGCGUCGCCGGACACUACCGGCAGACCGAAAUUCCUUUGUUGCCGCACCAGGGUUCUGAGGUUUAUAUACCGUUCGCCCUUGAACAACCCGGGACCAAUGGGUAGUGACGAUGAGGAGUUGCGUAAGGGCAAGUUAUGCUUUUCUGAAGUAUAAACUCAACUAGGUCGACAGAUGCAAGGUGUGCGAAAAAAACCAACUAGUUAUUUGUCCCCUAGAGUGGUGUUGGAAAGCCGAGGGGCACAUAUUCAUGUAAGGGGGGACCGCUUACUAGCCCGAUUGGUCCCUCACCCGGUGUAAAACGGCU'
        self.assertEqual(solve.dna_to_rna(dna), rna)

    # p.3
    def test_reverse_complement_dna(self):
        dna = 'ACGTTACAGGTAATGGCGTCCCAGGCGATCTCGGACGTAATTTTCATTAACACCTCGTGGGCAAACCTGTAGACGGCCATTAGACCGTTATTTCATCTAACTGGTGGGTTACCGAGACGCGAATCCTAGATTAACGGTCTTGGTCGTCCGGTTAGGATGGCATGCCAGTCCGAATAACCCTTTTTATAGGCCAAGATATTCCGTTTAAACATGTTATACTTAGACGCGTTGCAGTGGTTTAGAGGAAGTTACTCGCTCCCAATCTCCATCCTTGCGCGGCCGGGCTGACGCGGGGCGAGCACCCCAAGCTAGCAAAGCGCCTAGTACCTGTTAGGTCAAAGAAAGAGGCAACCACGACATCACCCTTTAACTGTTATCCCAGATCCCTAAGCGACTTCTCGATTTCAAAATTATAGACGGATGATGCAATAATTGGTCAGGGCTTCGGACAATCAAAAATGTGTGTCGTCTCCTATACTGAGACATTAGTTAATAGGCCTAGCGGTCATTTTGAGACGGAGTTCCCGAAAATCCGCTTGTATCACCATGAGGTCCTTCACGCTCCATAGTTCATAATTTTCCCCGTATCCATCGTCCCCGTGCCATAAGCTCACTACTGAGTCTCCCACGACAAGGGACTCACCAGGGGACTGAGATGGCATTCGAACTAGCACTGCTGTCGAATTCCCAACACAGGGTGGGGCGCCGCTAGAAAACACGATCGTCTCCTCAGACGGTACAACGATAGGTTTGACTCTGCACGTACTACACCCGATAACTAACTGATGCTGTTGTTCACTGGCTCCCCGTAACGATTGTCCTTGCAATTGTTTACTGGCCGGCTGTAGGTTGACGGTCGAAACACAAGGTTAAAACTCGCGTTGTGGTGCCGCCATATCACATGTCCAGAGCCGAC'
        reverse = 'GTCGGCTCTGGACATGTGATATGGCGGCACCACAACGCGAGTTTTAACCTTGTGTTTCGACCGTCAACCTACAGCCGGCCAGTAAACAATTGCAAGGACAATCGTTACGGGGAGCCAGTGAACAACAGCATCAGTTAGTTATCGGGTGTAGTACGTGCAGAGTCAAACCTATCGTTGTACCGTCTGAGGAGACGATCGTGTTTTCTAGCGGCGCCCCACCCTGTGTTGGGAATTCGACAGCAGTGCTAGTTCGAATGCCATCTCAGTCCCCTGGTGAGTCCCTTGTCGTGGGAGACTCAGTAGTGAGCTTATGGCACGGGGACGATGGATACGGGGAAAATTATGAACTATGGAGCGTGAAGGACCTCATGGTGATACAAGCGGATTTTCGGGAACTCCGTCTCAAAATGACCGCTAGGCCTATTAACTAATGTCTCAGTATAGGAGACGACACACATTTTTGATTGTCCGAAGCCCTGACCAATTATTGCATCATCCGTCTATAATTTTGAAATCGAGAAGTCGCTTAGGGATCTGGGATAACAGTTAAAGGGTGATGTCGTGGTTGCCTCTTTCTTTGACCTAACAGGTACTAGGCGCTTTGCTAGCTTGGGGTGCTCGCCCCGCGTCAGCCCGGCCGCGCAAGGATGGAGATTGGGAGCGAGTAACTTCCTCTAAACCACTGCAACGCGTCTAAGTATAACATGTTTAAACGGAATATCTTGGCCTATAAAAAGGGTTATTCGGACTGGCATGCCATCCTAACCGGACGACCAAGACCGTTAATCTAGGATTCGCGTCTCGGTAACCCACCAGTTAGATGAAATAACGGTCTAATGGCCGTCTACAGGTTTGCCCACGAGGTGTTAATGAAAATTACGTCCGAGATCGCCTGGGACGCCATTACCTGTAACGT'
        self.assertEqual(solve.reverse_complement_dna(dna), reverse)

    # p.4
    def test_fibonacci_rabbits(self):
        months, litter_size = 32, 3
        mature_rabbits = 108412748857
        self.assertEqual(solve.fibonacci_rabbits(months, litter_size), mature_rabbits)

    # p.5
    def test_find_gc_content(self):
        seq = '>Rosalind_0783CCGTGCCCCCCACAAAGGTTTGAACCCGGCCACTTGCCACGCAGGTCCCAATACGGAATCTGCATTGGAAAGGCACTGGAGCGATGACACTAAATTGTAAATGGCTCACGCATCTAGTTAACAGAACCTTGGCGTGATAGGTTTATGTCGCCCCTTTTGTGCGGTGCGGAGTAGTTGAGCCCACAACCGCGTTTCTGAGTCCGCCTAAGTGGCATTGGTTAGACCGATCAGACTGCAGTGAGTACGAGGTATCGTGGATTATGTCCGTATCGTCACAGCTTGAATGCAAATTAGATTACCTATCCCCTCAGACCAAAGTGATTCGGCTATGCAACAATTGTCCGGACCAGCTGACAAACTCGGCCCCACTTGAACGAAATGGGGAACTGTTTCGCGCAACTTGCTCACATAGAAGACACGCCCGTGCAAGCGCGAAAGCTCTATGATCCTAAGCAAATACCAAACGTTCTGATTGTAGTGGCCGATACCAATGGATGGGCATCACTCATGTTGACCGCACAAAATTTGTTTTCTTTAATTGGACATCTACTCTTAAACCGCCGGAGCTGCTCTTCTTTCAGTGCTTTGACGGTGCGCGTGTGTCTGAGTCCAGATACGTGAACTCTGTATGAGGAAAAAATAAAGCCGCAGTAGGTCGAGTTCCGACATAGGCTTCCCCGTATTCAGTGCCGGACTCCCCGAGCAGGATGGGGGAACGCCTTTCCTAAGCTAAACCTATTTTCCACTCGGTTTCCTCACCATGCTGCTCGGTATCCGTCGTCAAACTAGGTTAATGAATAGTTCCAGTGACTGGTGCTCCCAGTTCGATGTAGCTTATGGTTGGGCGGAGTTGAAACTCAGAATTTTTCGGG>Rosalind_6965GCGGAAAACGGTACCAGTCTTTCTATCTGGCACGTTCACCAGAGCCGTTTTGGCATGTACAGATATCAAAAAGCTGCCCGACACCCGCGGGGGTTGCAAGACCCCCAACAATTTGTACCCCTCGGAGCGATTGGCCGGTAGGACCCTCGTTCTGAGGTATTGTACACCGCCCTTATCGGCCGGCTTATCGGTATCGCTTCACAATTAACTCGGAAATGAGATATTCTAATTTAAACAGAATACATTATGAGTTGCGGACGTGTTTTCAAACAGTTCTACGTGTTCCCTATGTCGGTCTGCCAGTATAGTGGGATGCACAAGAGTCGGCCGGACGGTTTCCATGTGAATCTCTAAGGGAGACCTGAACGGACAACCGGGCATGCTAAATGTGAAAGTCAAAACCATGTAGCGTGCCCGTAGCCCGTACTATGCGTCGATTTCGGTCGCTGACGGGGTTGATTAATCCCATATAGGATGGTGCGAGAGGCGCTCGCGCCGTCCTATTCACGAGCGTATCGGTGTTCTCGCCTGCCCTCCTGCATCTCAGCCCGTGTATAAATCATGTCGGGGGCACTCAAAGTACTGGACGAACTAATTCGCGGGTCCCACACCTAACGTATAATCCCGGGTCAATTGTTCACTTGTCGCCCACACGCGATGTGCAGATCACGTGTCCTTCTCCTTTAACATTGGAACTGAACAGCGGTGTGTCTGGTTGTGTCGTTACGGTTTATCATGAAATTTAAGATTGGAGCTCCGGGCCTGTAATGGGCAAGCGACTGGGACTTTTCAAGAGATGATCCCACATGAGTGGTTCTTTTGAAGACGTGCGAGCCGCTGATTCCCATGCGTCGCAGCGCGCCTGGCAAGGAGCCGATAATCAGTTGTAGCACTGTGAGCATTGACTGGCACCCCGTGTGGAATGGTGTTCCAATACGGTTCAGCAACT>Rosalind_2659TACAGTACGTCCTATAAATATGACCTCAGCGCATTTTTACGACTCGGTTGTATGTCACGAAGATGGCGTGATTACACTCGCATGAAGATTACAATGGAGTTGGTTAGGCATGTGACTTATGGGTCCCGAGGCTGGCATCTAGTACATCATCTCAAAGCCTCCGCTGCCTGTACAACAAGTCGGCAGGTCCTTGTCTCTGAATTCGGCGAAGCTAGGGTTGCTTAAAGGTGGTAATACCTGATCACCTCACCACGCCATTCCCTCAAATGTGATTGCTACGGGCCCGCATAATTGGTAGCTCATTTCAATTTAACCTCGAGTCAGTGAACATTGAATATATGGCATTGAGGCCAGCAGGTCTGTCACTCGAGCACATCCCGAGTTAATCATCATTGGGTGGGTGGTTTCGTGACATTTATGGGCCAGATTTCTGAAATCAATGATCTCACTTTTGCATGACCTAGTCGTTGCGCTTCCTCGGAAGGTCTGTGGTACAGCGCCGCCTAGTTGTTGAGTAGCCGGAACTCTTCGCCACTTCCGTGCCCCGTTTGATGGACCCGCACGTGATTCTATGAGTCATTCGATAGCGCATCTAATCACCTCGGCGAGTGACGGGTGTAGTCTCTATGTCCGAGAGTAAGCTCTCGCAAACAAGGAGTTTGGCGTAGGTCTAAGACGCTGCCTGTTGCATCGCAAACCCTACGACCACCTGCTACTGACTGAGAGCCAATCGGAGTTACTGGGATTAAGTGACCTCTTGATAGGCTTATTATTTGGTCGCGACTTATAAACACGCAGGTTGGCCCATCCGGAGCGACGCTGCACGTGAGACATTGCGAAGGGCATGGTGACCTGTAGACTACCATCGAGTCTCCTGGTGCTCACGCCGGCCCCCCGGCGATCACTGTCTCTGTCCCCATTTGTTACTCCACATTACGCGATAGTATTCCTCGCGAGGACAGAGCATGACCTTAGATCCGAGGATTGACATGTG>Rosalind_8693AGCGATCTGACCATCGTAATGCACCCGGGCATCAATTGTTGCTGCACCACATTTACGATAAGTCCACGGCCAGTGGATAATGTGTTTGCAAACGCCGGGTTCCAGAGTTAGACATCTGAAAACCAATAGTTTCTTTCGTCACGGATTAGTAGGGTCATAGGTACTAGCCCAGTAAACCGGGGTGCTTAAACACTCGTCGAAACCGCCGGTCGCAGACTTCGTCTTGGCCCCAAAAATCGGATTACACGCTTGTACGAGTGCAGCCCAAAAAATATCCCTCCATGGCAGGATGGTGCGGATTAACTGGCGATCGAACATCTTGTCTCCTCTCTATTTTGCATAACGAAGCAAGACCGGAGGACACAAATCAGTTTCGCGCAGTTGCCTCCACAACGCAGCTCTGCTGTCACTTTATCCGCCTACGCACATGCGGAGAGTGACATCGATCGACAGCACCATTTGGCAAGGATAACCAAGTTAGAATAGTTTGATCTGATGTACTTCGGTAGATCTGCGTCCTGAAGCCCAAGGCAGGTTGGAGGTTATAGAGAAAGTCGAAAATACTACTAGGGCTGTCGATACAGTCTACTGACTACCGCTATATGGCTGACTAGAGGCGTAGCAGGAGTGGGGACGGCCATCTTTCAGGTACGTAACCTTCCTATAGCATAGCCGCCACGTAAAACGTCATCACTCTGCTGCCTTAAACACATTTGAGACACAGGTATGCATGGGTAGTACATCCGGGAGCAGATCTGATAGCGTTGACGGCAGTTGCACACCTTCAATCCGCTAAGAGTAAAAATAGGTTTGCCGTAAGTGACCCTCGACTCCATCATGTATTTCAAACGAGGGGAGTTTACGCCCCGGC>Rosalind_7571ATATTTACGCTGATCTTTACTAGGAAACCTTGACTCCACGCACTACGGCCATTTGAGGCGTAATTTAGTCAAGCTCTGCGTAAAACGCTTGACTGTAGCTTTTAGGCGGCAGACACGCTCCTCGGACCCGGACGCTGAGCGGGTTCCCGTACCTTGGGACCCAAATTAGTGAGAGCGCATATGAGCTGTAGGGATACTCGCCATCAGGGCGAGACTGTGCAGCAGCTAACCATTACTTCGGCGCGATTGGCTGCGAATGATAATCACCAGCTAAATTTAGTATTTTAATGGACTACAACTACTCCACATGGCTGAACGGCATGTCAAACTGACTGGAACCGCTAAGACGACTTGCGGGACATGGTTGTTTAATTGCAACAAGCGGAGCATAGGCGACGCTAGGTTGCCAACTTACTTTCACGTCTAAACATTGCCAACTCGAGAGATCGCCAAATAGTATTATGACTATGACGTGTGTTCGGTAGAGTCTAGACGCCTGTAATTGGGACTGGATACCGTCCAATACGCGACGGGCATAGTGTTGTGAGGGTTCGTGCAGTGCTGCAATTGTCCAGCACACCAGGTTAATGCCAACACGTTAATCATGAAGAAGGTACTTTCCTCTTAGTCCCACATGCCCTACGCGATTCAACCACTCCACCCCCGCTAGTACTGTCGACTTCCTTTACACTATAGGATCGGGCTAGAGCCACATTTTGCATTTACTGGAAGATGGCGGCTATATAGGAGTGCCTGACAGGGGCGAAAGACGGTCGGCCGTTCTGCTTGTCTAACGTTAGTAACAGGCGCCGTAGTTAGTCAACACCAAGTAGTCAAATGAAATAGTTCACTAGTGGTGGCAGGAACAGTGCACTGTTATACGCCTGACCCTTGATCATTCGGATGATTTCATGATTGTGCATCGATCTCGAGTAGGGACAACTGGGGCATACAAGAACGTTGCCCAGCCGGTCACATGCAGA>Rosalind_6946AAGTTCGTTCTATCGTTTCAAGACGGATCGTGTGTACCTACATGGTCCTCCCGCGGCGCGAGGATGTAGTTTCGCAGAGTACCCTACCAACGCGCACCAACACTTGTCGGAATGTCTTCATAGGACCGCCAGCTTGTCTTGAAGGGCTAAACCTTCGCGAACACTAGGACCCGCAATGTCTGGGTAAAGCTGGCCCATCTCTTGTAAGGTTTAAGTTTCGTGTCAGGATGAGACAAACACCGAAATCAAGCGCACAGTAGCTCGTCACAAAGGACTCCTCATAACGGAAATTAAGTCGCAAAAGAATCAGTGATAAAACTAGTCATCCGAATTTCCGACGTCCTTAACCTTGGTCCTAGGACACATCGGTTAGTGTCAGCTTCCCGATAAAGCGGATCGTGAATTCGTGACTGAGAGACGATACTAGCAAGGCTCTCCCCAAAGGGTTCGCCCATCCAATAGGCTGGCAGAAGCTGGTGAGCGAATAGCAACAGTCTGTATGCCAGCGCGTTTGCGATGTCCCGCGGGTCAGATTCGTCCTGGCCAATGCGGCCCCCGCCTGTGCGTCCCGTATTACCTACTCCTTACCTCCTAGGGGTGTAATTGGAGCACCCCTCGACTATGGCAATCCAACCATGACGTGTTAAGCGTTCCTATTCTTAAGAAACGGCAGTCTGGGTTCCTGATGGAAGTATGACAGTGCGACTGGCGCTCCATCCCGTGTTAAAAAGCCCGCACGCGTTGGTTCTCCGTTCGAAAGTTCTATTAGGGGCTCCGGCTAGACGCGATTTAGCCGTAGGACCAGTCTGTCCACGAATCTTCACGCAGGGCTCTGTACTCATAACCCACGTGTGAACCACTCTGTTATGAATATTAAC>Rosalind_5156GTGATGCGCTCGACTCGCCGATCAGTGGCAGCCGAGAAACAAGCTTCCCCTTCCAGTCGTCTAAGGTCACAGGAAACCGTGTAACTGTCGGCGCTCCTAAGTACCGCATGTGATTTCACCGCGAGGGGCACAAGAACCCTTCGGAGAGGACATCCCCACATTCAGCTTCAGTGGGGAAGGATTTCTAAGGGTGAGGTGTAGGCATAATACTAATGAAGCTACTGCGCACGTCGCGGCGATCTATCTCTAAAGATTCCTTGCATAAGACGCTCGTTGCGCTACGACCCTATAAGGGTGTAAATGCTAGAGTCGTCACAAGCTGTAACCGCCGGGCTCGCTAGCGTCTCATGACATGCCTAGCGGTGCTAACCCATTTGGTGATCTATACGAATCCTCAACGGTCCCACCATAAGTAAAAGCCCACGTAGGAGAAGACTCTTTTGAGTCGTACTTGCAATCAGTTGGTGCGCAATTAATCACTACTTTGACAAGGGACGTGGTGGGCTAACCCACAACTTCGATGTCTTATGATAAAAGTAAGCATTCGCGAAACGCGTACAGAAGAGCTCGTTCTGTGACGAAGTTCCGGTCCATCCGATTGATGCTACCTGTCTTATCTAGTCAGCCGCCGTCTTAATTTCTCCCTCTACATATCTGAAAACCCTAGGAGTACGACATCCAAGAAAATAGCCGCGCAATTCCTAGCTCTTACTCGCATTCGAGGACATTGGTTCCGACGCCATACGAAAATCAATTTAGGATCATATATTATATATGATCTATTGCTCCGTCCCGATATACCATTAACCAGGTTCCAAACCGCACGTAAGGTGCTTGTGGCACCAATAATTAAACTACTACGATTCGCTATCTACCCGGTAACTAGTTCATTGTATTTACTGATCAATAGCCGATTATTCGGGGTAAGTACACAGGGGCGCACAACGGCTGC'
        gc = 'Rosalind_6965 51.633298208640674'
        self.assertEqual(solve.find_gc_content(seq), gc)

    # p.6
    def test_find_hamming_distance(self):
        first, second = 'GAGCCTACTAACGGGAT', 'CATCGTAATGACGGCCT'
        score = 7
        self.assertEqual(solve.find_hamming_distance(first, second), score)

    # p.7
    def test_proba_dominant_allele(self):
        homozygous_dominant, heterozygous, homozygous_recessive = 2, 2, 2
        proba = 0.78333
        self.assertAlmostEqual(solve.proba_dominant_allele(
            homozygous_dominant, heterozygous, homozygous_recessive), proba, 5)

    # p.8
    def test_rna_to_protein(self):
        rna = 'AUGGCCAUGGCGCCCAGAACUGAGAUCAAUAGUACCCGUAUUAACGGGUGA'
        protein = 'MAMAPRTEINSTRING'
        self.assertEqual(solve.rna_to_protein(rna), protein)

    # p.9
    def test_find_motif(self):
        seq, subseq = 'GATATATGCATATACTT', 'ATAT'
        index_list = [2, 4, 10]
        self.assertEqual(solve.find_motif(seq, subseq), index_list)

    # p.10
    def test_find_profile(self):
        sequences = '>Rosalind_1ATCCAGCT>Rosalind_2GGGCAACT>Rosalind_3ATGGATCT>Rosalind_4AAGCAACC>Rosalind_5TTGGAACT>Rosalind_6ATGCCATT>Rosalind_7ATGGCACT'
        profile = (
            'ATGCAACT',
            {
                'A': '5 1 0 0 5 5 0 0',
                'T': '1 5 0 0 0 1 1 6',
                'C': '0 0 1 4 2 0 6 1',
                'G': '1 1 6 3 0 1 0 0'
            }
        )
        self.assertEqual(solve.find_profile(sequences, test=True), profile)

    # p.11
    def test_mortal_rabbits(self):
        months, life_span = 6, 3
        number_rabbits = 4
        self.assertEqual(solve.mortal_rabbits(months, life_span), number_rabbits)

    # p.12
    def test_overlap_graphs(self):
        sequences = '>Rosalind_0498AAATAAA>Rosalind_2391AAATTTT>Rosalind_2323TTTTCCC>Rosalind_0442AAATCCC>Rosalind_5013GGGTGGG'
        overlaps = [
            'Rosalind_0498 Rosalind_2391',
            'Rosalind_2391 Rosalind_2323',
            'Rosalind_0498 Rosalind_0442'
        ]
        self.assertEqual(solve.overlap_graphs(sequences), overlaps)

    # p.13
    def test_expected_offspring(self):
        genotype_population = '1 0 0 1 0 1'
        expected_value = 3.5
        self.assertEqual(solve.expected_offspring(
            genotype_population), expected_value)

    # p.14
    def test_shared_motif(self):
        sequences = '>Rosalind_1GATTACA>Rosalind_2TAGACCA>Rosalind_3ATACA'
        longest_motif = 'TA'
        self.assertEqual(solve.shared_motif(
            sequences, test=True), longest_motif)

    # p.15
    def test_calc_proba_heterozygous(self):
        generation, number_with_trait = 2, 1
        proba_heterozygous = 0.684
        self.assertAlmostEqual(solve.calc_proba_heterozygous(
            generation, number_with_trait), proba_heterozygous, 3)

    # p.16
    def test_find_protein_motif(self):
        dataset = 'A2Z669,B5ZC00,P07204_TRBM_HUMAN,P20840_SAG1_YEAST'
        protein_dict = {
            'B5ZC00': [85, 118, 142, 306, 395],
            'P20840_SAG1_YEAST': [79, 109, 135, 248, 306, 348, 364, 402, 485, 501, 614],
            'P07204_TRBM_HUMAN': [47, 115, 116, 382, 409]
        }
        self.assertEqual(solve.find_protein_motif(dataset), protein_dict)

    # p.17
    def test_protein_to_mrna(self):
        sequence = 'MA'
        num_possible = 12
        self.assertEqual(solve.protein_to_mrna(sequence), num_possible)

    # p. 18
    def test_open_reading_frames(self):
        sequence = 'AGCCATGTAGCTAACTCAGGTTACATGGGGATGACCCCGCGACTTGGATTAGAGTCTCTTTTGGAATAAGCCTGAATGATCCGAGTAGCATCTCAG'
        possible_proteins = {
            'MLLGSFRLIPKETLIQVAGSSPCNLS',
            'MGMTPRLGLESLLE',
            'M',
            'MTPRLGLESLLE'
        }
        self.assertEqual(solve.open_reading_frames(
            sequence), possible_proteins)

    # p.19
    def test_n_permutations(self):
        number = 3
        permutations = (
            6,
            [(1, 2, 3),
             (1, 3, 2),
             (2, 1, 3),
             (2, 3, 1),
             (3, 1, 2),
             (3, 2, 1)]
        )
        self.assertEqual(solve.n_permutations(number), permutations)

    # p.20
    def test_protein_mass(self):
        sequence = 'SKADYEK'
        mass = 821.392
        self.assertAlmostEqual(solve.protein_mass(sequence), mass, 3)

    # p.21
    def test_reverse_palindrome(self):
        sequence = 'TCAATGCATGCGGGTCTATATGCAT'
        kmer_list = [(5, 4), (7, 4), (17, 4), (18, 4),
                     (21, 4), (4, 6), (6, 6), (20, 6)]
        self.assertEqual(solve.reverse_palindrome(sequence), kmer_list)
    # p.22

    def test_rna_splicing(self):
        sequences = '>Rosalind_10ATGGTCTACATAGCTGACAAACAGCACGTAGCAATCGGTCGAATCTCGAGAGGCATATGGTCACATGATCGGTCGAGCGTGTTTCAAAGTTTGCGCCTAG>Rosalind_12ATCGGTCGAA>Rosalind_15ATCGGTCGAGCGTGT'
        protein = 'MVYIADKQHVASREAYGHMFKVCA'
        self.assertEqual(solve.rna_splicing(sequences, test=True), protein)

    # p.23
    def test_lexicographic_permutations(self):
        alphabet, string_length = 'T A G C', 2
        permutations = ['AA', 'AC', 'AG', 'AT', 'CA', 'CC', 'CG',
                        'CT', 'GA', 'GC', 'GG', 'GT', 'TA', 'TC', 'TG', 'TT']
        self.assertEqual(solve.lexicographic_permutations(
            alphabet, string_length), permutations)

    # p.24
    def test_longest_subsequence(self):
        number, sequence = 5, '5 1 4 2 3'
        longest_subsequences1 = (['1', '2', '3'], ['5', '4', '2'])
        longest_subsequences2 = (['1', '2', '3'], ['5', '4', '3'])
        result = solve.longest_subsequence(number, sequence)
        self.assertTrue(
            result == longest_subsequences1 or result == longest_subsequences2)

    # p.25
    def test_shortest_superstring(self):
        sequences = '>Rosalind_56ATTAGACCTG>Rosalind_57CCTGCCGGAA>Rosalind_58AGACCTGCCG>Rosalind_59GCCGGAATAC'
        superstring = 'ATTAGACCTGCCGGAATAC'
        self.assertEqual(solve.shortest_superstring(
            sequences, test=True), superstring)
        self.assertEqual(solve.shortest_superstring(
            sequences, brute=False, test=True), superstring)

    # p.26
    def test_matching_graph(self):
        sequence = 'AGCUAGUCAU'
        number_matching = 12
        self.assertEqual(solve.matching_graph(sequence), number_matching)

    # p.27
    def test_partial_permutations(self):
        number, subset = 21, 7
        partial = 51200
        self.assertEqual(solve.partial_permutations(number, subset), partial)


if __name__ == '__main__':
    unittest.main()
