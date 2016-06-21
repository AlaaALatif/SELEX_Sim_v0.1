import time
import random
import linecache
from itertools import izip, imap
import operator
from collections import OrderedDict
from scipy import stats
import numpy as np

class Mutation:

# randomly mutates each seq copy in mutatedPool using the corresponding
# random positions. random mutations can be circular (i.e. A to A)
# mutated seqs are added to the amplfied pool of wild-type seqs
    def mutGen(self, mutatedPool, amplfd_seqs):
        alpha = ['A', 'T', 'G', 'C'] #DNA alphabet
        for seq in mutatedPool:
            mutNum = mutatedPool[seq][0] #no. of mutation instances
            i=1
            while(mutatedPool[seq][0] > 0): #for each seq copy
                mut = mutatedPool[seq][i] # grab the i-th mutation instance for a given seq
                for j in mut: #for each mutation position in the instance
                    randNuc = alpha[np.random.randint(4)] #pick random nt
                    mutatedSeq = seq[:(j-1)] + randNuc + seq[j:] #mutate seq using random nt
                if mutatedSeq not in amplfd_seqs: # add mutated seq to the amplified pool
                    amplfd_seqs.setdefault(mutatedSeq, []).append(1)
                else:
                    amplfd_seqs[mutatedSeq][0] += 1
                mutatedPool[seq][0] -= 1 #decrement no. of seq copies to be mutated
                i += 1
        outPool = amplfd_seqs
        print("Mutation has been carried out")
        return outPool
