fasta = open('dna.example.fasta')

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
    tempsequence = ''
    identifier = ''
    for line in fastafile:
        if line[:1] != '>':
            tempsequence += line.rstrip()
        else:
            if identifier != '':
                print("Appending Sequence: ", tempsequence)
                sequences.append({'sequence': tempsequence, 'identifier': identifier})
                identifier = line.split('|')[0].split('>')[1]
                tempsequence = ''
            else: #first line
                # print("First Line Found!")
                identifier = line.split('|')[0].split('>')[1]
                tempsequence = ''

    sequences.append({'sequence': tempsequence, 'identifier': identifier}) #at the end, append the last sequence
    print("Sequences within get_sequences: ", sequences)
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

def split_by_n( seq, n ):
    """A generator to divide a sequence into chunks of n units."""
    while seq:
        yield seq[:n]
        seq = seq[n:]

def find_frames(sequences):
    orf_list = {'frame_1_list': [], 'frame_2_list': [], 'frame_3_list': []}
    frame_1 = []
    frame_2 = []
    frame_3 = []
    for sequence in sequences:
        #to do - find all ORFs
        frame_1 = list(split_by_n(sequence['sequence'], 3))
        orf_list['frame_1_list'].append({'identifier': sequence['identifier'], 'sequence': frame_1})
        frame_2 = [sequence['sequence'][:1]]
        frame_2 += (list(split_by_n(sequence['sequence'][1:], 3)))
        orf_list['frame_2_list'].append({'identifier': sequence['identifier'], 'sequence': frame_2})
        frame_3 = [sequence['sequence'][:2]]
        frame_3 += (list(split_by_n(sequence['sequence'][2:], 3)))
        orf_list['frame_3_list'].append({'identifier': sequence['identifier'], 'sequence': frame_3})
    return orf_list

def find_ORFs(framelist):
    orf_list = []
    start_index = 0
    end_index = 0
    for i, frame in enumerate(framelist):
        if frame.lower() == 'atg':
            start_index = i
        elif (frame.lower() == 'taa' or frame.lower() == 'tag' or frame.lower() == 'tga') and start_index != 0:
            end_index = i
            orf_list.append({'orf': framelist[start_index:end_index], 'index': start_index})
            end_index = 0
            start_index = 0
    return orf_list

def get_ORFS(framedict):
    frame_1_list = []
    frame_2_list = []
    frame_3_list = []
    for frame in framedict['frame_1_list']:
        frame_1_list.append({'identifier': frame['identifier'], 'orfs': find_ORFs(frame['sequence'])})
    for frame in framedict['frame_2_list']:
        frame_2_list.append({'identifier': frame['identifier'], 'orfs': find_ORFs(frame['sequence'])})
    for frame in framedict['frame_3_list']:
        frame_3_list.append({'identifier': frame['identifier'], 'orfs': find_ORFs(frame['sequence'])})
    return {'frame_1_list': frame_1_list, 'frame_2_list': frame_2_list, 'frame_3_list': frame_3_list}

def longest_ORF(orflist):
    longest_value = 0
    longest_index = 0
    longest_identifier = ''
    for index, orf_entry in enumerate(orflist):
        for orf in orf_entry['orfs']:
            if len(orf['orf']) > longest_value:
                longest_value = len(orf['orf'])
                longest_index = orf['index']
                longest_identifier = orf_entry['identifier']
    result = {'longest_value': longest_value * 3, 'index': longest_index+1, 'identifier': longest_identifier}
    print("Longest ORF: ", result)
    return result

def longest_by_ID(orflist, ID):
    longest_value = 0
    longest_index = 0
    for index, orf_entry in enumerate(orflist):
        if orf_entry['identifier'].lower() == ID.lower():
            for orf in orf_entry['orfs']:
                if len(orf['orf']) > longest_value:
                    longest_value = len(orf['orf'])
                    longest_index = orf['index']
                    longest_identifier = orf_entry['identifier']
    if longest_value == 0 or longest_index == 0:
        result = "NO Results Found"
    else:
        result = {'longest_value': longest_value * 3, 'index': longest_index+1, 'identifier': longest_identifier}
    print("Longest ORF: ", result)
    return result

sequences_list = get_sequences(fasta)
print("*** Sequences List ***")
print(sequences_list)
print("\n**********************\n")
print("Number of Records: ", len(sequences_list))
length_list = sequence_lengths(sequences_list)
print("*** Sequence Lengths ***")
print("Longest Sequence: ", length_list['longest'])
print("Shortest Sequence: ", length_list['shortest'])
print("Longest List: ", length_list['longest_list'])
print("\n**********************\n")
print("Shortest List: ", length_list['shortest_list'])
print("Find ORFs")
frame_list = find_frames(sequences_list)

orf_list = get_ORFS(frame_list)
# print("Frame 1 list: ", orf_list['frame_1_list'])
print('\n****** Longest Frame 1 *******\n')
longest_ORF(orf_list['frame_1_list'])
# print("Frame 2 list: ", orf_list['frame_2_list'])
print('\n****** Longest Frame 2 *******\n')
longest_ORF(orf_list['frame_2_list'])
# print("Frame 3 list: ", orf_list['frame_3_list'])
print('\n****** Longest Frame 3 *******\n')
longest_ORF(orf_list['frame_3_list'])

print('\n****** Longest Frame 3 by ID *******\n')
longest_by_ID(orf_list['frame_3_list'], 'g')
# for ORF in orf_list['frame_1_list']:
#     print("ORF: ", ORF['orfs'])
# for ORF in orf_list['frame_2_list']:
#     print("ORF: ", ORF['orfs'])
# for ORF in orf_list['frame_3_list']:
#     print("ORF: ", ORF['orfs'])
# print("Max Frame 1: ", max(orf_list['frame_1_list']))
# print("Max Frame 2: ", max(orf_list['frame_2_list']))
# print("Max Frame 3: ", max(orf_list['frame_3_list']))
