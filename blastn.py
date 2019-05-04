import os, sys, subprocess, copy, re
import argparse
import fileinput
import time
from multiprocessing import Pool, Manager
from seeding import *


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

        self.db = database_preparing(self.args.databse)
        self.
