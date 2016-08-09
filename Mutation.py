import time
import random
import linecache
from itertools import izip, imap
import operator
from collections import OrderedDict
from scipy import stats
import numpy as np

class Mutation:
##### NEED TO INCLUDE GENERATION NUMBER DISTRIBUTION; ACCOUNT FOR PARTIAL AMPLIFICATIONS OF WILD TYPE AND MUTANTS
# randomly mutates each seq copy in mutatedPool using the corresponding
# random positions. random mutations can be circular (i.e. A to A)
# mutated seqs are added to the amplfied pool of wild-type seqs
    def mutGen(self, seqLength, mutatedPool, amplfdSeqs, aptamerSeqs, alphabetSet):
        alphabetSize = len(alphabetSet)
        for seqIdx in mutatedPool:
            mutatedSeq = apt.pseudoAptamerGenerator(seqIdx, alphabetSet, seqLength) #fetch sequence
            mutNum = mutatedPool[seqIdx][0] #no. of mutation instances
            i=1
            while(mutatedPool[seqIdx][0] > 0): #for each seq copy
                mut = mutatedPool[seqIdx][i] # grab the i-th mutation instance for a given seq
                for j in mut: #for each mutation position in the instance
                    randNuc = alphabetSet[random.randint(0, alphabetSize - 1)] #pick random nt
                    mutatedSeq = mutatedSeq[:(j-1)] + randNuc + mutatedSeq[j:] #mutate seq using random nt
                mutatedSeqIdx = apt.pseudoAptamerIndexGenerator(mutatedSeq, alphabetSet, seqLength)
                if mutatedSeqIdx not in amplfdSeqs: # add mutated seq to the amplified pool
                    amplfdSeqs.setdefault(mutatedSeqIdx, []).append(1)
                    amplfdSeqs.setdefault(mutatedSeqIdx, []).append(d.hamming_func(mutatedSeq, aptamerSeqs))

                    amplfdSeqs.setdefault(mutatedSeqIdx, []).append(d.bias_func(mutatedSeq, seqLength))

                else:
                    amplfdSeqs[mutatedSeqIdx][0] += 1
                mutatedPool[seq][0] -= 1 #decrement no. of seq copies to be mutated
                i += 1 #move to next randPos array
        print("Mutation has been carried out")
        return amplfdSeqs
