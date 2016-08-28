fasta = open('dna.example.fasta','r')

# function no longer neccessary, since it can be computed from the sequence info
# def find_records(fastafile):
#     records = 0
#     for line in fastafile:
#         if line[:1] == '>':
#             print('record found');
#             records += 1
#     return records
#
# records = find_records(fasta)
# print("Total Number of records:", records)

def get_sequences(fastafile):
    sequences = []
    tempsequence = {'identifier': '', 'sequence': '', 'length': ''}
    for line in fastafile:
        if line[:1] != '>':
            tempsequence['sequence'] += line.rstrip()
        else:
            if tempsequence['identifier'] != '':
                # print("Appending Sequence: ", tempsequence)
                sequences.append(tempsequence)
                tempsequence['identifier'] = line.split('|')[0].split('>')[1]
                tempsequence['sequence'] = ''
            else: #first line
                # print("First Line Found!")
                tempsequence['identifier'] = line.split('|')[0].split('>')[1]
                tempsequence['sequence'] = ''
    sequences.append(tempsequence) #at the end, append the last sequence
    for sequence in sequences:
        sequence['length'] = len(sequence['sequence'])
    return sequences

def sequence_lengths(sequences):
    longest_list = []
    shortest_list = []
    longest = 0
    shortest = sequences[0]['length']
    for sequence in sequences:
        if sequence['length'] == longest:
            longest_list.append(sequence)
        elif sequence['length'] > longest:
            longest_list = [sequence]
            longest = sequence['length']
        if sequence['length'] == shortest:
            shortest_list.append(sequence)
        elif sequence['length'] < shortest:
            shortest_list = [sequence]
            shortest= sequence['length']
    return {'longest': longest, 'longest_list': longest_list, 'shortest': shortest, 'shortest_list': shortest_list}


sequences_list = get_sequences(fasta)
print("*** Sequences List ***")
print(sequences_list)
print("\n**********************\n")
print("Number of Records: ", len(sequences_list))
length_list = sequence_lengths(sequences_list)
print("*** Sequence Lengths ***")
print("Longest Sequence: ", length_list['longest'])
print("Shortest Sequnce: ", length_list['shortest'])
print("Longest List: ", length_list['longest_list'])
print("\n**********************\n")
print("Shortest List: ", length_list['shortest_list'])
