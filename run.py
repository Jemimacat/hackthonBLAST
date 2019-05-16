import os,sys,subprocess
import numpy as np
import pkg_resources
from seeding import query_seed_preparing,database_seed_preparing,merge_scan_and_scoring
from extend import extending_words,extending_words_positions

#one_query = 'GACAGCGACGCCGCGAGCCAGAAGATGGAGCCGCGGGCGCCGTGGATAGAGCAGGAGGGGCCGGAGTATTGGGACCAGGAGACACGGAATATGTTGGCCCACTCACAGACTGACCGAGCGAACCTGGGGACCCTGCGCTACTACTACAACCAGAGCGAGGACGGTTCTCACACCATCCAGATAATGTATGGCTGCGACGTGGGGCCGGACGGGCGCTTCCTCCGCGTACCGGCAGG'
test_database = 'A_nuc.fasta'
test_query = 'query.txt'
seeds = query_seed_preparing(test_query)
db,genes = database_seed_preparing(test_database)
scores  = merge_scan_and_scoring(seeds,db)
segment_hits = {}

for query in seeds.keys():
    one_query_dict = seeds[query]
    one_scores = scores[query]
    for word in one_scores.keys():
        positions,db_positions = extending_words_positions(word, query, one_scores, one_query_dict, db, max_gap=5, step=1)
        # print(db_positions)
        one_segment_hits = extending_words(word, query, positions, db_positions, genes)
        segment_hits.update({query:one_segment_hits})

## output result

result1 = open('result_mod1.txt','w+')
#result2 = open('result_mod2.txt','w+')

for query in segment_hits.keys():
    one_segment_hits = segment_hits[query]
    result1.writelines('query: '+query+'\n')
    count = 1
    for score in sorted(one_segment_hits.keys(),reverse=True):
        for consensus in sorted(one_segment_hits[score].keys(),reverse=True):
            for q_length in sorted(one_segment_hits[score][consensus].keys()):
                if q_length < 0:
                    continue
                for gene in sorted(one_segment_hits[score][consensus][q_length].keys()):
                    for q_segment in sorted(one_segment_hits[score][consensus][q_length][gene].keys()):
                        for db_segment in sorted(one_segment_hits[score][consensus][q_length][gene][q_segment].keys()):
                            for items in sorted(one_segment_hits[score][consensus][q_length][gene][q_segment][db_segment].keys()):
                                (layer,q_start,q_end,db_start,db_end) = items
                                result1.writelines(str(count) + '\tGene:' + gene + '\tLength:' + str(q_length) + '\tCon:' + str(consensus) + '\tScore:' + str(score) + '\tQuery:' + str(q_start) + ',' + str(q_end) + '\tRef:' + str(db_start) + ',' + str(db_end) +'\n')
                                result1.writelines('Ref: '+db_segment +'\n')
                                result1.writelines('     '+layer +'\n')
                                result1.writelines('Que: '+q_segment +'\n')
                                result1.writelines('\n')
                                count += 1
result1.close()

# if __name__ == '__main__':
#     sys.exit(main())