#algo
import random
import numpy as np

class BlastSequence(object):
    '''
    Generate random DNA sequences
    '''
    def __init__(self):
        self.sequence = self.generate_sequence()

    def __repr__(self):
        print('{:>6} {:<25} {:<5}'.format('number', 'sequence', 'length'))
        return ("\n").join(['{:>6} {:<25} {:<5}'.format(num, str(seq), len(seq)) for num, seq in enumerate(self.sequence)])

    def generate_sequence(self):
        number_of_sequence = 2 #int(random.random()*10) + 1
        sequences = []
        letters = ['A', 'C', 'G', 'T']
        for number in range(number_of_sequence):
            sequence_length = 15 #int(random.random()*25)
            if sequence_length < 10:
                sequence_length = sequence_length + 10
            generated_sequence = []
            for i in range(sequence_length):
                generated_sequence.append(np.random.choice(letters))
            generated_sequence = ('').join(generated_sequence)
            sequences.append(generated_sequence)
        return sequences

    #Calculate Hamming distance between two sequences
    def calc_hamming(self):
        score = 0
        for i, y in enumerate(self.sequence[0]):
            first_seq_letter = self.sequence[0][i]
            second_seq_letter = self.sequence[1][i]
            letter_is_equal = first_seq_letter == second_seq_letter
            if letter_is_equal:
                comparison_str = "is equal to"
            else:
                score += 1
                comparison_str = "is not equal to"
            print('{} {} {}'.format(first_seq_letter, comparison_str, second_seq_letter))
        print("Hamming distance: {}".format(score))
        print("{:.1f}% similar".format((1-(score/len(self.sequence[0])))*100))
        




