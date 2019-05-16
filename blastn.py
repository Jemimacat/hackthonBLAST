import os, sys, subprocess, copy, re
import argparse
import fileinput
import time
from multiprocessing import Pool, Manager
from seeding import query_seed_preparing,database_seed_preparing,merge_scan_and_scoring
from extend import extending_words,extending_words_positions


class BLASTn():

    def __init__(self): 
        self._parser = argparse.ArgumentParser()
    
    def parameter(self):

        print(time.asctime( time.localtime(time.time())) + " Reading parameters...")
        self._parser.add_argument('-db','--database',help='Input database in .fasta format.')
        self._parser.add_argument('-i','--in_query',help='Input query sequence.')
        self._parser.add_argument('-o','--outdir',help='Output directory (path).',type=str,default='.')
        self._parser.add_argument('-w','--seed_len',help='Seed length',type=int,default=11)
        self._parser.add_argument('-T','--threshold',help='Score threshold',type=int,default=11)
        self._parser.add_argument('-nt','--proc_num',help='Max number of processes (int)',type=int,default=4)

        self.args = self._parser.parse_args()

        if (self.args.database == None or self.args.in_query == None):
            raise Exception('Error: query sequence and database must be provided!')
    
    def init_database_and_query(self):

        self.db, self.genes = database_seed_preparing(self.args.databse)
        self.seeds = query_seed_preparing(self.args.in_query)

    def scaning_and_extending(self,max_gap=5, step=1):
        self.segment_hits = {}
        self.scores  = merge_scan_and_scoring(self.seeds,self.db)
        for query in self.seeds.keys():
            one_query_dict = self.seeds[query]
            one_scores = self.scores[query]
            for word in one_scores.keys():
                positions,db_positions = extending_words_positions(word, query, one_scores, one_query_dict, self.db, max_gap=max_gap, step=step)
                one_segment_hits = extending_words(word, query, positions, db_positions, self.genes, max_gap=max_gap)
                self.segment_hits.update({query:one_segment_hits})

    def output_result(self):
        file = self.args.outdir + '/result.txt'
        result = open(file,'w+')
        for query in self.segment_hits.keys():
            one_segment_hits = self.segment_hits[query]
            result.writelines('query: '+query+'\n')
            count = 1
            for score in sorted(one_segment_hits.keys(),reverse=True):
                for consensus in sorted(one_segment_hits[score].keys(),reverse=True):
                    for q_length in sorted(one_segment_hits[score][consensus].keys()):
                        for gene in sorted(one_segment_hits[score][consensus][q_length].keys()):
                            for q_segment in sorted(one_segment_hits[score][consensus][q_length][gene].keys()):
                                for db_segment in sorted(one_segment_hits[score][consensus][q_length][gene][q_segment].keys()):
                                    for items in sorted(one_segment_hits[score][consensus][q_length][gene][q_segment][db_segment].keys()):
                                        (layer,q_start,q_end,db_start,db_end) = items
                                        result.writelines(str(count) + '\tGene:' + gene + '\tLength:' + str(q_length) + '\tCon:' + str(consensus) + '\tScore:' + str(score) + '\tQuery:' + str(q_start) + ',' + str(q_end) + '\tRef:' + str(db_start) + ',' + str(db_end) +'\n')
                                        result.writelines('Ref: '+db_segment +'\n')
                                        result.writelines('     '+layer +'\n')
                                        result.writelines('Que: '+q_segment +'\n')
                                        result.writelines('\n')
                                        count += 1
        result.close()



