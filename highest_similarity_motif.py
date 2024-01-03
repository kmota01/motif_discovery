from fasta_dictionary import fasta_dictionary
import time

start = time.time()

input = 'generated_sequences.fasta'
my_dict=fasta_dictionary(input)
genomic_dataset=list(my_dict.values())

def calculate_similarity(sequence_a, sequence_b, start_a, start_b, length):
    count = 0
    for i in range(length):
        if sequence_a[start_a + i] == sequence_b[start_b + i]:
            count += 1
    return count

def high_similarity_motifs(sequence_a,sequence_b, min_len, max_len):

    motif_length = min_len
    max_similarity = 0
    best_motifs = []

    while min_len <= motif_length <= max_len:
        
        for i in range(len(sequence_a) - motif_length + 1):
            for j in range(len(sequence_b) - motif_length + 1):
                similarity = calculate_similarity(sequence_a, sequence_b, i, j, motif_length)
                current_motif = sequence_a[i:i+motif_length]

                if similarity >= motif_length*0.5:

                    best_motifs.append([current_motif, similarity, i, j])
                        
        motif_length += 1

    return best_motifs

def is_motif_present_in_all_sequences(dataset, motif, threshold):

    motif_length = len(motif)
    sum,average_sim = 0,0

    for sequence in dataset:
        sequence_length = len(sequence)

        # Iterate over all possible k-mers in the sequence
        for i in range(sequence_length - motif_length + 1):
            subsequence = sequence[i:i + motif_length]

            similarity = calculate_similarity(subsequence, motif, 0, 0, motif_length) * 100 / motif_length
            
            

            if similarity >= threshold:
                sum += similarity
                break 

        else:
            return False,0
          
    average_sim = sum/len(dataset)
    return True, average_sim
    # If all sequences passed the check, return True and the average similarity


high_similarity_motif,motif_dictionary={},{}
filtered_motifs={}

#calculate high similarity motifs between all the pairwise sequences of the genomic dataset
for i in range(51-1):
    list_name = f"high_similarity_motif_{i}"

    high_similarity_motif = high_similarity_motifs(genomic_dataset[i],genomic_dataset[i+1],min_len=8,max_len=13)

    #remove duplicates and sort
    high_similarity_motif = [high_similarity_motif[i][0] for i in range(len(high_similarity_motif)-1)]
    high_similarity_motif = sorted(list(set(high_similarity_motif)))

    
    motif_dictionary[list_name] = []
    motif_dictionary[list_name].append(high_similarity_motif)


for i in range(50-1):
    for motif1 in list(motif_dictionary.values())[i][0]:
        result,sim = is_motif_present_in_all_sequences(genomic_dataset,motif1,threshold=80)
        if result:
                if not motif1 in filtered_motifs:
                    filtered_motifs[motif1]=round(sim,2)
                    



sorted_dict = {key: filtered_motifs[key] for key in sorted(filtered_motifs, key=filtered_motifs.get, reverse=True)}
print(sorted_dict)
max_value = max(filtered_motifs.values())
max_keys = [key for key, value in filtered_motifs.items() if value == max_value]
print(max_keys)


end=time.time()

print(f'Elapsed time: {end-start}')

