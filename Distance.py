import time
import random
import linecache
from itertools import izip, imap, islice
import operator
from collections import OrderedDict
from scipy import stats

def hamming_func(str1, str2):
   assert len(str1) == len(str2)
   ne = operator.ne
   return sum(imap(ne, str1, str2))


class Distance:

    def hamming_func(self, str1, str2):
        assert len(str1) == len(str2)
        ne = operator.ne
        return sum(imap(ne, str1, str2))


    def biasedHamming_initLib(self, seqLen, optimseqs, pool_file, scale, partition):
        seqs=np.zeros((scale, 3))
        seqIdx = 0
        with open(pool_file) as p:
            while True:
                next_n_seqs = list(islice(p, partition))
                if not next_n_seqs:
                    break
                for i, seq in enumerate(next_n_seqs):
                    seqs[seqIdx][0] += 1 #increment count
                    seqs[seqIdx][1] = hamming_func(optimseqs, seq[:-1])
                    for nt in seq[:-1]:
                        if(nt == 'C') or (nt == 'T'):
                            seqs[seqIdx][2] += 1 #increment no. of pyrimidines
                    seqs[seqIdx][2] = 0.1*(2*seqs[seqIdx][2] - seqLen)/seqLen #compute bias
                    seqIdx += 1
        p.close()

        return seqs

##TEST AREA
"""
d = Distance()
seqs = d.biasedHamming_initLib(20, 'AAAAAAAAAAAAAAAAAAAA', 'random_initLib_1KM', 100000000, 1000000)
"""


    def seqsHamming(self, optimseqs, pool_file):
        seqs=np.zeros((scale, 3))
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

# NOTE: DISTANCES ARE NOW APPENDED AS A 3RD VALUE, NOT 2ND DUE TO BIAS
# This method computes the hamming distance of each sequence to the optimum set
# and their PCR bias. It requires sequence length as 3rd argument
   def seqsHamming_and_Bias(self, optimseqs, pool_file, seqLen):
      seqs={}
      seqsfile = open(pool_file)
      seq=seqsfile.readline()
      while(seq):
         if seq not in seqs:
            seqs.setdefault(seq, []).append(1) #append seq count
            seqs.setdefault(seq, []).append(0) #append seq bias
            for nt in seq:
                if(nt == 'C') or (nt == 'T'):
                    seqs[seq][1]+=1 #increment no. of pyrimidines
            seqs[seq][1] = 0.1*(2*seqs[seq][1] - seqLen)/seqLen #compute bias
            for optimseq in optimseqs:
               hammdist = hamming_func(optimseq, seq) #compute distance
               seqs.setdefault(seq, []).append(hammdist) #append distance
         else:
            seqs[seq][0]+=1 #increment count
         seq=seqsfile.readline()
      seqsfile.close()
      print("distance calculations and bias scoring completed")
      return seqs
