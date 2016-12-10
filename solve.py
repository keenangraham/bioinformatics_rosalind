#solve bioinformatics problems on Rosalind (http://rosalind.info/)
#each function returns an answer to specified question given proper input

def nucleotide_count(sequence):
    #p.1 - print ACGT count separated by spaces
    #alt. soln.: return sequence.count('A'), sequence.count('C'), etc. 
    seq = sequence.upper()
    seq_dict = {}
    for letter in seq:
        if letter not in seq_dict.keys():
            seq_dict[letter] = 0
        seq_dict[letter] += 1
    print(seq_dict['A'], seq_dict['C'], seq_dict['G'], seq_dict['T'])

def dna_to_rna(sequence):
    #p.2 - convert DNA to RNA by replacing T with U
    print(sequence.replace('T', 'U'))

def reverse_complement_dna(sequence):
    #p.3 - output reverse complement DNA
    # A:T, C:G
    complement_map = {'A':'T', 'C':'G', 'T':'A', 'G':'C'}
    print(''.join(list(map(lambda x: complement_map[x], sequence)))[::-1])

def fibonacci_rabbits(months, litter_size):
    #p.4 - output no. of pairs of mature rabbits after given no. of months, litter_size
    mature_rabbits = 1
    baby_rabbits = 1
    for i in range(months+1):
        if i == 0: 
            print('month: {}, mature: {}, baby: {}'.format(i, mature_rabbits, baby_rabbits))
            continue
        if i == 1:
            mature_rabbits = mature_rabbits + baby_rabbits
            baby_rabbits = 0
            print('month: {}, mature: {}, baby: {}'.format(i, mature_rabbits, baby_rabbits))
            continue
        mature_rabbits, baby_rabbits = mature_rabbits+baby_rabbits,mature_rabbits*litter_size        
        print('month: {}, mature: {}, baby: {}'.format(i, mature_rabbits, baby_rabbits))
    print('total pairs: {}'.format(int(mature_rabbits/2)))

def find_gc_content(fasta):
    #p.5 - parses >Rosalind_xxxx FASTA file, returns name and GC percent of sequence with max GC content
    sequences = fasta.split('>Rosalind_')
    sequences.pop(0)
    sequence_dict = {}
    for seq in sequences:
        id = seq[:4]
        sequence = seq[4:]
        seq_length = len(sequence)
        gc_content = ((sequence.count('G') + sequence.count('C'))/seq_length)*100
        sequence_dict[id] = gc_content
    max_gc = max(sequence_dict, key=sequence_dict.get)
    print('Rosalind_{} has the highest GC content at {}%'.format(max_gc, sequence_dict[max_gc]))

def find_hamming_distance(first_sequence, second_sequence):
    #p.6 - output Hamming distance between first and second sequence
    #number of nucleotides that differ
    score = 0
    for i, y in enumerate(first_sequence):
        if first_sequence[i] != second_sequence[i]:
            score += 1
    print(score)

def proba_dominant_allele(homozygous_dominant, heterozygous, homozygous_recessive):
    #p.7 - given proportions of population homozygous dominant/recessive and heterozygous
    #return probability that two random mating partners will produce offspring with dominant allele
    population_size = homozygous_dominant + heterozygous + homozygous_recessive
    default_values = {'dom': homozygous_dominant, 'het': heterozygous, 'rec': homozygous_recessive}
    
    def reset_dict():
        mating_dict = {
        'population_size': population_size,
        'dom': default_values['dom'],
        'het': default_values['het'],
        'rec': default_values['rec'],
        'previous_choice': None
        }
        return mating_dict

    def update_population_probability(mating_dict):
        population_size = mating_dict['population_size']
        mating_dict['prob_dom'] = mating_dict['dom']/population_size
        mating_dict['prob_het'] = mating_dict['het']/population_size
        mating_dict['prob_rec'] = mating_dict['rec']/population_size
        return mating_dict

    def change_population_size(mating_dict):
        removed_individual = mating_dict['previous_choice']
        mating_dict['population_size'] -= 1
        mating_dict[removed_individual] -= 1
        return mating_dict
        
    def calc_mating_probability(mating_dict):
        if mating_dict['previous_choice'] is None:
            mating_dict = update_population_probability(mating_dict)
            return mating_dict
        else:
            new_mating_dict = change_population_size(mating_dict)
            mating_dict = update_population_probability(new_mating_dict)
            return mating_dict

    #calculate probability of dominant allele given all combinations of mating partners
    #dominant-dominant, dominant-heterozygous, dominant-recessive
    #heterozygous-heterozygous, heterozygous-recessive
    #recessive-recessive

    dom_allele_map = {
        'dom_dom': 1.0,
        'dom_het': 1.0,
        'het_dom': 1.0,
        'dom_rec': 1.0,
        'het_het': 0.75,
        'het_rec': 0.5,
        'rec_het': 0.5,
        'rec_rec': 0,
        'rec_dom': 1.0
       }

    #calculate probability of every combination of mating partners
    #dom-dom, dom-het, dom-rec
    #het-dom, het-het, het-rec
    #rec-dom, rec-het, rec-rec
    choice_list = ['dom', 'het', 'rec']
    proba_dict = {}
    for first_choice in choice_list:
       mating_dict = reset_dict()
       mating_dict = calc_mating_probability(mating_dict)
       first_proba = mating_dict['prob_'+first_choice]
       mating_dict['previous_choice'] = first_choice
       mating_dict = calc_mating_probability(mating_dict)
       for second_choice in choice_list:
           second_proba = mating_dict['prob_'+second_choice]
           proba_dict[first_choice + '_' + second_choice] = first_proba * second_proba
           print(first_choice, first_proba, second_choice, second_proba)
    sum = 0
    for key in proba_dict.keys():
        sum = sum + (proba_dict[key]*dom_allele_map[key])
    print(sum)

def rna_to_protein(sequence):
    #p.7 - print protein sequence given DNA
    #generate codon table (borrowed code)
    bases = ['U', 'C', 'A', 'G']
    codons = [a+b+c for a in bases for b in bases for c in bases]
    amino_acids = 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'
    codon_table = dict(zip(codons, amino_acids))
    rna_chunked = [sequence[i:i+3] for i in range(0, len(sequence), 3)]
    protein = "".join(list(map(lambda x: codon_table[x], rna_chunked)))
    print(protein)
 
def find_motif(sequence, motif):
    #p.8 - find index of motif given sequence
    #alt. soln. if sequence[i:].startswith(motif)
    motif_length = len(motif)
    index_list = []
    for i, y in enumerate(sequence):
        if sequence[i:i+motif_length] == motif:
            index_list.append(i+1)
    print(index_list)
 
def find_profile(sequences):
    #p.9 - generate profile and consensus string given many sequences
    import numpy as np
    sequence_list = sequences.split('>Rosalind_')
    sequence_list = [seq[4:] for seq in sequence_list if seq != '']
    sequence_length = len(sequence_list[0])
    sequence_list = np.array([[i for i in seq] for seq in sequence_list])
    seq_dict = {}
    profile_list = ['A', 'C', 'G', 'T']
    consensus_seq = ''
    for i in range(sequence_length):
        unique, counts = np.unique(sequence_list[:,i], return_counts=True)
        count_dict = dict(zip(unique, counts))
        seq_dict[i] = count_dict
        max_key = max(count_dict, key=count_dict.get)
        consensus_seq = consensus_seq + max_key
       
    print(consensus_seq)
    print('')
    for letter in profile_list:
        letter_row = []
        for item in seq_dict.items():
            try:
                letter_row.append(item[1][letter])
            except:
                letter_row.append(0)
        print('{}: {}'.format(letter, ' '.join(str(x) for x in letter_row)))
        print('')

def fibonacci(max_value):
    #example of Fibonacci sequence using yield statement
    def fibonacci_generator(max_value):
        a, b = 0, 1
        i = 0
        while i < max_value:
            yield b
            a, b = b, a+b
            i += 1
    print(*fibonacci_generator(max_value), sep='\n')

def mortal_rabbits(months, life_span):
    #p.10 - print number of rabbits after given number of months
    #and given lifespan
    life_span -= 1
    baby_rabbit_pairs = 1
    mature_rabbit_pairs = 0
    baby_rabbit_list = [] 
    for i in range(months-1):
        baby_rabbit_list.append(baby_rabbit_pairs)
        baby_rabbit_pairs, mature_rabbit_pairs = mature_rabbit_pairs,\
            mature_rabbit_pairs + baby_rabbit_pairs\
            - (baby_rabbit_list[i-life_span] if i >= life_span else 0)
    print(baby_rabbit_pairs + mature_rabbit_pairs)

def overlap_graphs(sequences):
    #p.11 - print list of sequences where last three letters of string one
    #match first three letters of string two
    sequences = sequences.split(">Rosalind_")
    sequences = [seq for seq in sequences if seq != '']
    for seq_one in sequences:
        seq_id_one = seq_one[:4]
        for seq_two in sequences:
            seq_id_two = seq_two[:4]
            if seq_id_one != seq_id_two:
                if seq_one[4:].startswith(seq_two[-3:]):
                    print("Rosalind_{} Rosalind_{}".format(seq_id_two, seq_id_one))

def expected_offspring(genotype_population):
    #p.12 - print expected number of offspring with dominant allele
    #given number of mating couples with certain alleles
    #AA-AA, AA-Aa, AA-aa, Aa-Aa, Aa-aa, aa-aa
    
    #hardcode probability of dominant allele given mating couples
    dominant_probability = [1, 1, 1, 0.75, 0.5, 0]
    genotype_population = [int(gen) for gen in genotype_population.split()]
    zip_proba_genotype = zip(dominant_probability, genotype_population)
    expected_value = 0
    for z in zip_proba_genotype:
        a, b = z
        expected_value = expected_value + 2*(a*b)
    print(expected_value)

def shared_motif(sequences):
    #p.13 - print longest shared motif contained in all sequences
    sequences = sequences.split(">Rosalind_")
    sequences = [seq[4:] for seq in sequences if seq != '']
    longest_motif = ''
    all_contain_motif = True
    for index, letter in enumerate(sequences[0]):
        motif_length = 1
        all_contain_motif = True
        while all_contain_motif:
    #could instead start with longest sequence and shorten
    #for more efficiency
            motif = sequences[0][index:index+motif_length]
            for seq in sequences:
                if motif not in seq:
                    all_contain_motif = False
                    break
                if all_contain_motif:
                    if len(motif) > len(longest_motif):
                        longest_motif = motif
                        motif_length += 1
                        if motif_length > len(sequences[0]):
                            break
    print(longest_motif)

def calc_proba_heterozygous(generation, number_with_trait):
    #p.14 - print probability of seeing number_with_trait individuals
    #heterozygous for two independent traits in specified generation
    #assuming each descendant mates with a heterozygous individual
    #and has two offspring
 
    #key insight: any individual mating with a heterozygous individual
    #will have 1/4 chance of offspring with heterozygous alleles
    #probability of no heterozygous for each descendant is 1-(1/4)
    
    #when k = 2, n= 1: 1 - ((0.75**2)**((2**k)*(1/2)))
    #alt. soln. - use scipy.stats.binom.pmf(n, k, p=0.25) for n in range(n, 2**k)
    #alt. soln. 2:  1 - scipy.stats.binom.cdf(n-1, 2**k, 0.25)
    import math
    proba_heterozygous = 0
    for num in range(number_with_trait, (2**generation)+1):
        number_offspring = 2**generation
        number_no_trait = number_offspring-num
        number_of_combinations = math.factorial(number_offspring)/(math.factorial(num)*(math.factorial((number_no_trait))))
        proba_of_number_with_trait = (3/4)**(number_no_trait) * (1/4)**(num)
        proba_heterozygous = proba_heterozygous + (proba_of_number_with_trait * number_of_combinations)
    print(proba_heterozygous)
 
def find_protein_motif(dataset):
    #p.15 - return proteins in given dataset with
    #N-glycosylation motif: N{P}[ST]{P}
    
    #pull sequences from UniProt
    import regex as re
    import requests
    ids = [i for i in dataset.split(',')]
    protein_dict = {}
    for i in ids:
        url = 'http://www.uniprot.org/uniprot/{}.fasta'.format(i)
        r = requests.get(url)
        seq = r.text.split("SV=")[1][2:]
        seq = ("").join(seq.split("\n"))
        pattern = re.compile('N[^P][ST][^P]')
        positions = [index.start() + 1 for index in pattern.finditer(seq, overlapped=True)]
        protein_dict[i] = positions
        if protein_dict[i]:
            print(i)
            print(*protein_dict[i], sep=' ')

def protein_to_mrna(sequence):
    #p.16 - print number of possible mRNA sequences that
    #would produce given protein string, modulo 1,000,000

    #generate codon table (borrowed code)
    bases = ['U', 'C', 'A', 'G']
    codons = [a+b+c for a in bases for b in bases for c in bases]
    amino_acids = 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'
    codon_table = dict(zip(codons, amino_acids))
    possible_sequences = 1
    for seq in sequence: 
        possible_sequences *= len([i for i in codon_table.values() if i == seq])
    #multiply by three because three possible stop codons
    print((3 * possible_sequences) % 1000000)
 
def open_reading_frames(sequence):
    #p. 17 - given DNA sequence return all possible proteins
    #made from open reading frames of strand and
    #reverse strand
    import regex as re
    def dna_to_rna(sequence):
        return sequence.replace('T', 'U')
    def reverse_complement_dna(sequence):
        complement_map = {'A':'U', 'C':'G', 'U':'A', 'G':'C'}
        return ('').join(list(map(lambda x: complement_map[x], sequence)))[::-1]
    def rna_to_protein(sequence):
        bases = ['U', 'C', 'A', 'G']
        codons = [a+b+c for a in bases for b in bases for c in bases]
        amino_acids = 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'
        codon_table = dict(zip(codons, amino_acids))
        rna_chunked = [sequence[i:i+3] for i in range(0, len(sequence), 3)]
        protein = "".join(list(map(lambda x: codon_table[x], rna_chunked)))
        return protein
    def reading_frame(sequence):
        seq = ''
        for codon in [sequence[i:i+3] for i in range(0, len(sequence), 3)]:
            if codon in ["UAA", "UAG", "UGA"]:
                seq += codon
                return seq
            else:
                seq += codon
        return 0
    def start_index(sequence):
        pattern = re.compile('AUG')
        return [index.start() for index in pattern.finditer(sequence, overlapped=True)]
    def valid_sequences(sequence, start_index):
        seq_list = []
        for i in start_index:
            seq = reading_frame(sequence[i:])
            if seq and seq not in seq_list:
                seq_list.append(rna_to_protein(seq).replace('*', ''))
            return seq_list                 

    pattern = re.compile('AUG')
    sequence = dna_to_rna(sequence)
    sequence_reverse = reverse_complement_dna(sequence)
    start_forward = start_index(sequence)
    start_reverse = start_index(sequence_reverse)
    forward_seqs = valid_sequences(sequence, start_forward)
    reverse_seqs = valid_sequences(sequence_reverse, start_reverse)
    print(*set(forward_seqs + reverse_seqs), sep='\n')

def n_permutations(number):
    #p.18 - print number of possible
    #permutations from 1...number
    #as well as all of the permutations
    import itertools
    permutations = list(itertools.permutations([i for i in range(1, number+1)]))
    print(len(permutations))
    for perm in permutations:
        print(*perm, sep=" ")
 
def protein_mass(sequence):
    #p.19 print monoisotopic mass of protein

    #monoisotopic mass table
    mass_table = {
        "A":71.03711,
        "C":103.00919,
        "D":115.02694,
        "E":129.04259,
        "F":147.06841,
        "G":57.02146,
        "H":137.05891,
        "I":113.08406,
        "K":128.09496,
        "L":113.08406,
        "M":131.04049,
        "N":114.04293,
        "P":97.05276,
        "Q":128.05858,
        "R":156.10111,
        "S":87.03203,
        "T":101.04768,
        "V":99.06841,
        "W":186.07931,
        "Y":163.06333
        }
    mass = 0
    for letter in sequence:
        mass += mass_table[letter]
    print(mass)
 
def reverse_palindrome(sequence):
    #p.20 - print position and length of every
    #palindrome with length between 4 and 12
    def reverse_complement_dna(sequence):
        complement_map = {'A':'T', 'C':'G', 'T':'A', 'G':'C'}
        return ''.join(list(map(lambda x: complement_map[x], sequence)))[::-1]
    kmer_dict = {}
    counter = 0
    for length in range(4,13):
        for index, letter in enumerate(sequence):
            seq = sequence[index:index + length]
            if len(seq) == length:
                kmer_dict[counter] = {"seq":seq, "position":index+1, "length":length, "rc":reverse_complement_dna(seq)}
                counter += 1
                for key, value in kmer_dict.items():
                    if value['seq'] == value["rc"]:
                        print(value["position"], value["length"])

def rna_splicing(sequences):
    #p.21 - return protein string resulting from DNA sequence after
    #introns (substrings provided) have been removed
    def dna_to_rna(sequence):
        return sequence.replace('T', 'U')
    def rna_to_protein(sequence):
        bases = ['U', 'C', 'A', 'G']
        codons = [a+b+c for a in bases for b in bases for c in bases]
        amino_acids = 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'
        codon_table = dict(zip(codons, amino_acids))
        rna_chunked = [sequence[i:i+3] for i in range(0, len(sequence), 3)]
        protein = "".join(list(map(lambda x: codon_table[x], rna_chunked)))
        return protein
    sequences = [seq[4:] for seq in sequences.split(">Rosalind_") if seq != '']
    seq = sequences[0]
    for intron in sequences[1:]:
        seq = seq.replace(intron, '')
    print(rna_to_protein(dna_to_rna(seq)).replace('*', ''))


                        



