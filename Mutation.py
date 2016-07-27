import time
import random
import linecache
from itertools import izip, imap
import operator
from collections import OrderedDict
from scipy import stats
import numpy as np

class Mutation:
##### NEED TO IMPLEMENT METHOD TO APPEND NOVEL MUTATED SEQ INFO 
##### INTO AMPLIFIED POOL AND FETCH SEQ INFO OF MUTATED SEQ THAT
##### IS ALREADY PRESENT IN AMPLIFIED POOL
# randomly mutates each seq copy in mutatedPool using the corresponding
# random positions. random mutations can be circular (i.e. A to A)
# mutated seqs are added to the amplfied pool of wild-type seqs
    def mutGen(self, seqLength, mutatedPool, amplfdSeqs, alphabetSet):
        alphabetSize = len(alphabetSet)
        for seqInfo in mutatedPool:
            seq = apt.pseudoAptamerGenerator(seqInfo[0], alphabetSet, seqLength) #fetch sequence
            mutNum = mutatedPool[seqInfo][0] #no. of mutation instances
            i=1
            while(mutatedPool[seqInfo][0] > 0): #for each seq copy
                mut = mutatedPool[seqInfo][i] # grab the i-th mutation instance for a given seq
                for j in mut: #for each mutation position in the instance
                    randNuc = alphabetSet[random.randint(0, alphabetSize - 1)] #pick random nt
                    mutatedSeq = seq[:(j-1)] + randNuc + seq[j:] #mutate seq using random nt
                mutatedSeqIdx = apt.pseudoAptamerIndexGenerator(mutatedSeq, alphabetSet, seqLength)
                if mutatedSeqIdx not in amplfdSeqs[:,0]: # add mutated seq to the amplified pool
                    amplfdSeqs.setdefault(mutatedSeq, []).append(1)
                else:
                    amplfdSeqs[mutatedSeq][0] += 1
                mutatedPool[seq][0] -= 1 #decrement no. of seq copies to be mutated
                i += 1
        outPool = amplfd_seqs
        print("Mutation has been carried out")
        return outPool
