import random

def hamming_distance(motif_a,motif_b):
    count=0
    for i,j in zip(range(len(motif_a)),range(len(motif_b))):
        if motif_a[i] == motif_b[j]:
            continue
        else:
            count+=1
    dist=len(motif_a)-count 
    return dist

def generate_genomic_sequence(length):
    bases = ['A', 'T', 'C', 'G']
    sequence = ''.join(random.choice(bases) for _ in range(length))
    return sequence

def create_genomic_dataset(num_sequences,min_length,max_length):
    dataset = []
    for _ in range(num_sequences):
        sequence_length = random.randint(min_length, max_length)
        sequence = generate_genomic_sequence(sequence_length)
        dataset.append(sequence)
    return dataset

def write_fasta_file(file_path, sequences):
    with open(file_path, 'w') as fasta_file:
        for i, sequence in enumerate(sequences, start=1):
            fasta_file.write(f">Sequence_{i}\n{sequence}\n")

def calculate_hamming_distance(motif1, motif2):
    return sum(base1 != base2 for base1, base2 in zip(motif1, motif2))

def mutate_genomic_motif(prototype_motif, mutation_rate):

    motif_length = len(prototype_motif)
    num_mutations = int(mutation_rate * motif_length)

    mutated_motif = list(prototype_motif)

    for _ in range(num_mutations):
        # Mutate the nucleotide
        index_to_mutate = random.randint(0, motif_length - 1)
        mutated_motif[index_to_mutate] = random.choice(['A','T','C','G'])
    
    mutated_motif = ''.join(mutated_motif)
    return mutated_motif

def genomic_sequence_with_motif(dataset,motif):
    seq_motif,motif_valid=[],[]
    for i,sequence in enumerate(dataset,start=1):
        random_pos=random.randint(0,len(sequence)-1)
        mutated_motif=mutate_genomic_motif(motif,mutation_rate=0.2)
        similarity = hamming_distance(motif,mutated_motif)*100/len(motif)
        motif_valid.append(sequence[:random_pos-1] + ' ' + mutated_motif + ' ' + sequence[random_pos:] + ' ' + str(similarity))
        seq_motif.append(sequence[:random_pos-1] + mutated_motif + sequence[random_pos:])
    return seq_motif, motif_valid


# Random Sequences:
num_sequences = random.randint(5000,10000)
fasta_file_path = 'generated_sequences.fasta'
min_sequence_length=50
max_sequence_length=100

# Random Motif
motif_length = random.randint(8,13)

genomic_dataset = create_genomic_dataset(num_sequences, min_sequence_length, max_sequence_length)
random_motif = generate_genomic_sequence(motif_length)
sequence_with_motif, motif_validation = genomic_sequence_with_motif(genomic_dataset,random_motif)
write_fasta_file(fasta_file_path, sequence_with_motif)
with open('genomic_sequences.fasta', 'w') as fasta_file:
    fasta_file.writelines(random_motif+'\n'+'\n'.join(motif_validation))

print(f"FASTA file '{fasta_file_path}' has been created with {num_sequences} sequences.")
print(f"The random motif for your dataset: {random_motif} with length {motif_length}")