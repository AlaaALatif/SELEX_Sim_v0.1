import time
import random
import linecache
from itertools import izip, imap, product
import operator
from collections import OrderedDict
from scipy import stats 

class Aptamers:

# Generate all possible sequences 
    def aptamerGenerator(self, alphabetSet, seqLength):
        initialLibrary = product(alphabetSet, repeat=seqLength)
        
        return initialLibrary





# choose a random subset of sequences to be aptamers from the initial pool
# input params are the number of aptamers to choose and the initial pool 
# sequence file
# returns total number of seqs in initial pool and set of aptamers 
    def optimumSeqs(self, num_optm_seqs, init_pool_file):                                 
        with open(init_pool_file) as f: #open initial library
            for totalseqs, l in enumerate(f): #calculate total num of seqs
                pass
            totalseqs+=1
        f.close()
        print(totalseqs)
        optimseqs_idx={}
        optimseqs={}
        i=1
        while(i <= num_optm_seqs):
            optimseq_idx=random.randint(0,totalseqs) #random seq index
            if optimseq_idx not in optimseqs_idx:
                optimseqs_idx.setdefault(optimseq_idx, [])                          
                optimseq = str(open(init_pool_file).readlines()[optimseq_idx])
                optimseqs.setdefault(optimseq, [])
            i+=1
        print(optimseqs)
        return optimseqs, totalseqs
