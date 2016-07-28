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
    def mutGen(self, seqLength, mutatedPool, amplfdSeqs, aptamerSeqs, alphabetSet):
        alphabetSize = len(alphabetSet)
        amplfdSeqs = amplfdSeqs.tolist()
        for seqInfo in mutatedPool:
            mutatedSeq = apt.pseudoAptamerGenerator(seqInfo[0], alphabetSet, seqLength) #fetch sequence
            mutNum = mutatedPool[seqInfo][0] #no. of mutation instances
            i=1
            while(mutatedPool[seqInfo][0] > 0): #for each seq copy
                mut = mutatedPool[seqInfo][i] # grab the i-th mutation instance for a given seq
                for j in mut: #for each mutation position in the instance
                    randNuc = alphabetSet[random.randint(0, alphabetSize - 1)] #pick random nt
                    mutatedSeq = mutatedSeq[:(j-1)] + randNuc + mutatedSeq[j:] #mutate seq using random nt
                mutatedSeqIdx = apt.pseudoAptamerIndexGenerator(mutatedSeq, alphabetSet, seqLength)
                if mutatedSeqIdx not in amplfdSeqs[:,0]: # add mutated seq to the amplified pool
                    mutatedSeqInfo = np.zeros(4)
                    mutatedSeqInfo[0] = mutatedSeqIdx
                    mutatedSeqInfo[1] = 1
                    mutatedSeqInfo[2] = d.hamming_func(mutatedSeq, aptamerSeqs)
                    mutatedSeqInfo[3] = d.bias_func(mutatedSeq, seqLength)
                    amplfdSeqs.append(mutatedSeqInfo)
                else:
                    for seqInfo in amplfdSeq:
                        if(seqInfo[0] == mutatedSeqIdx):
                            amplfdSeqs[amplfdSeqs.index(seqInfo)][1] += 1
                            break
                mutatedPool[seq][0] -= 1 #decrement no. of seq copies to be mutated
                i += 1
        outPool = np.asarray(amplfdSeqs)
        print("Mutation has been carried out")
        return outPool
