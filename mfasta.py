fasta = open('dna.example.fasta','r')

def find_records(fasta):
    records = 0
    for line in fasta:
        if line[:1] == '>':
            print('record found');
            records += 1
    return records

records = find_records(fasta);
print("Total Number of records:", records)
