from Bio import SeqIO
import numpy as np
import fileinput
from scoring import score_nt_seq
from extend import extending_words,extending_words_positions

# Genering seeds
def seed_list_of_query_generating(query_seq,w=11):
    one_seed = {}
    num = len(query_seq) - w + 1
    if num >= 0:
        for i in range(num):
            if not query_seq[i:i+w] in one_seed:
                one_seed[query_seq[i:i+w]] = []
            one_seed[query_seq[i:i+w]].append(i)
    return one_seed

def query_seed_preparing(query_file,w=11):
    seeds = {}
    for line in fileinput.input(query_file):
        one_query = line.rstrip('\n')
        one_seed = seed_list_of_query_generating(one_query,w)
        seeds.update({one_query:one_seed})
    return seeds

# Preparing database
def database_seed_preparing(fasta,w=11):
    database = {}
    genes = {}
    for db_record in SeqIO.parse(fasta, "fasta"):
        db_seq = str(db_record.seq)
        db_id = str(db_record.id).split(':')[1]
        genes[db_id] = db_seq
        num = len(db_seq) - w + 1
        if num >= 0:
            for i in range(num):
                if not db_seq[i:i+w] in database.keys():
                    database[db_seq[i:i+w]] = {}
                    if not db_id in database[db_seq[i:i+w]].keys():
                        database[db_seq[i:i+w]][db_id] = [i]
                    else:
                        database[db_seq[i:i+w]][db_id].append(i)          
    return database,genes

def one_query_scan_and_scoring(one_word_dict,db_dict,threshold=11):
    scores = {}
    for word in one_word_dict.keys():
        for db_record in db_dict.keys():
            score = score_nt_seq(word,db_record)
            if score >= threshold:
                if not word in scores.keys():
                    scores[word] = {}
                if not db_record in scores[word].keys():
                    scores[word][db_record] = []
                scores[word][db_record] = score
    return scores

def merge_scan_and_scoring(word_dict,db_dict,threshold=11):
    scores = {}
    for query in word_dict.keys():
        one_score = one_query_scan_and_scoring(word_dict[query],db_dict,threshold)
        scores.update({query:one_score})
    return scores
    
def build_hit_matrix(scores,query_dict,db_dict,genes):
    
    hit_matrix = {}
    for (query,one_score) in scores.items():
        hit_matrix[query] = {}
        for (gene,seq) in genes.items():
            matrix = np.zeros((len(query),len(seq)))
            for word in one_score.keys():
                for db_record in one_score[word]:
                    for i in range(len(word)):
                        for query_pos in query_dict[query][word]:
                            if gene in db_dict[db_record]:
                                for db_pos in db_dict[db_record][gene]:
                                    if word[i] == db_record[i]:
                                        matrix[query_pos+1,db_pos+i] = 1
            if matrix.any():
                hit_matrix[query][gene] = matrix

    return hit_matrix



# one_query = 'GACAGCGACGCCGCGAGCCAGAAGATGGAGCCGCGGGCGCCGTGGATAGAGCAGGAGGGGCCGGAGTATTGGGACCAGGAGACACGGAATATGTTGGCCCACTCACAGACTGACCGAGCGAACCTGGGGACCCTGCGCTACTACTACAACCAGAGCGAGGACGGTTCTCACACCATCCAGATAATGTATGGCTGCGACGTGGGGCCGGACGGGCGCTTCCTCCGCGTACCGGCAGG'
# test_database = 'A_nuc.fasta'
# test_query = 'query.txt'
# seeds = query_seed_preparing(test_query)
# db,genes = database_seed_preparing(test_database)
# scores  = merge_scan_and_scoring(seeds,db)
# segment_hits = {}

# for query in seeds.keys():
#     one_query_dict = seeds[one_query]
    


# one_scores = scores[one_query]
# #print(scores)
# word = 'AGAGCAGGAGG'
# word2 = 'GAGCAGGAGGG'
# one_query_dict = seeds[one_query]
# positions,db_positions = extending_words(word, one_query, one_scores, one_query_dict, db)
# print(positions)
# print(db_positions)
#hit_matrix = build_hit_matrix(scores,seeds,db,genes)

# for query in scores.keys():
#     print(query+':')
#     for word in scores[query].keys():
#         print('\t'+word+':')
#         for pair_w in scores[query][word].keys():
#             print('\t\t'+pair_w+'\t'+str(scores[query][word][pair_w]))

# for query in seeds.keys():
#     print(query+':')
#     for word in seeds[query].keys():
#         print('\t'+word+':')
#         for pos in seeds[query][word]:
#             print('\t\t'+str(pos))

# for word in db.keys():
#     print(word+':')
#     for gene in db[word].keys():
#         print('\t'+gene+':')
#         for pos in db[word][gene]:
#             print('\t\t'+str(pos))

# for gene in genes.keys():
#     print(gene+':')
#     print(genes[gene])