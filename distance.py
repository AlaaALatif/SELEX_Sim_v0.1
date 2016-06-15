import time
import random
import linecache
from itertools import izip, imap
import operator
from collections import OrderedDict
from scipy import stats

class distance

def hamming_func(str1, str2):
   assert len(str1) == len(str2)
   ne = operator.ne
   return sum(imap(ne, str1, str2))


def seqs_hamming(optimseqs, pool_file):
    seqs={}
    seqsfile = open(pool_file)
    seq=seqsfile.readline()
    while(seq):
        if seq not in seqs:
            seqs.setdefault(seq, []).append(1)
            for optimseq in optimseqs:
               hammdist = hamming_func(optimseq, seq)
               seqs.setdefault(seq, []).append(hammdist)
        else:
            seqs[seq][0]+=1
        seq=seqsfile.readline()
    seqsfile.close()
    print("distance calculations passed")
    return seqs
