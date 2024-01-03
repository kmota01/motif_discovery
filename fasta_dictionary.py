#python3

def fasta_dictionary(input_file):

    dict = {}
    headers, sequences = '', ''

    with open(input_file,'r') as infile: 
        lines = infile.read().splitlines()
        for f in lines:
            if f.startswith(">"):
                count = 0
                headers = f
            else:
                count += 1
                if count == 1:
                    sequences = f
                else:
                    sequences = sequences + f
            dict[headers] = sequences
        return dict