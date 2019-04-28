from Bio import SeqIO

# Nucleotides: Genering seeds
def seed_list_of_query_generating(query_seq,w=11):
    seeds = {}
    num = len(query_seq) - w + 1
    if num >= 0:
        for i in range(num):
            if not query_seq[i:i+w] in seeds:
                seeds[query_seq[i:i+w]] = []
            seeds[query_seq[i:i+w]].append(i)
    return seeds

# Nucleotides: Preparing database
def database_preparing(fasta,w=11):
    database = {}
    for db_record in SeqIO.parse(fasta, "fasta"):
        db_seq = db_record.seq
        num = len(db_seq) - w + 1
        if num >= 0:
            for i in range(num):
                if not db_seq[i:i+w] in database:
                    database[db_seq[i:i+w]] = []
                database[db_seq[i:i+w]].append('|'.join(db_record.id,str(i)))
    return database   


test = 'AKDJAMSHJWHZKSHJFLSFJDDLSLMCSSKNDJF'
seed_list_of_query_generating(test)

