import time
import random
import numpy as np
import linecache
from itertools import izip, imap, islice, product
import operator
from collections import OrderedDict
from scipy import stats
import Aptamers, Distance

D = Distance.Distance()
Apt = Aptamers.Aptamers()

class Selection: 
    def stochasticHammingSelection(self, alphabetSet, seqLength, seqPool, aptPool, selectionThreshold, totalSeqNum, roundNum):
        if(roundNum == 0):
            slctdSeqs = {} #initialize seqInfo matrix
            selectedSeqs = 0
            print("parameters for selection have been initialized")
   # select top 20 percent in terms of hamm distance
            
            while(selectedSeqs <= selectionThreshold):
                randHammScore = random.randint(0, seqLength)
                randSeqIdx = random.randint(0, totalSeqNum - 1)
                randSeq = Apt.pseudoAptamerGenerator(randSeqIdx, alphabetSet, seqLength)  
                randSeqDist = D.hamming_func(randSeq, aptPool)
                if(randSeqDist <= randHammScore):
                    if(randSeqIdx in slctdSeqs):
                        slctdSeqs[randSeqIdx][0] += 1
                    else:

                        slctdSeqs.setdefault(randSeqIdx, []).append(1) #add to selected pool
                        slctdSeqs.setdefault(randSeqIdx, []).append(randSeqDist) #distance
                        slctdSeqs.setdefault(randSeqIdx, []).append(D.bias_func(randSeq, seqLength)) #bias

                    selectedSeqs += 1 #increment sampled no.
                else:
                    continue

            print("sequence selection has been carried out")
            return slctdSeqs
        else:
            slctdSeqs = {}
            sampledSeqs = 0
            print("parameters for selection have been initialized")
   # select top sequences in terms of hamm distance
            while(sampledSeqs <= selectionThreshold):
                randHammScore = random.randint(0, seqLength)
                randPoolIdx = random.randint(0, totalSeqNum-1)
                randSeqIdx = seqPool.keys()[randPoolIdx]
                if(seqPool[randSeqIdx][0] <= randHammScore):
       # THIS IS TAKING TOO LONG # UPDATE: NOT ANYMORE
                    if(randSeqIdx not in slctdSeqs):
                        slctdSeqs.setdefault(randSeqIdx, []).append(1) #add to selected pool
                        slctdSeqs[randSeqIdx][1] = seqPool[randSeqIdx][1]
                        slctdSeqs[randSeqIdx][2] = seqPool[randSeqIdx][2]
                    else:
                        slctdSeqs[randSeqIdx] += 1 #increment count
                            
                    sampledSeqs += 1 #increment sampled no.
                else:
                    continue
            print("sequence selection has been carried out")
            return slctdSeqs


   

    def definiteSelection(self, seqs, selection_rate, totalseqs):
        seqs_ordrd = OrderedDict(sorted(seqs.items(), key=lambda seqs: seqs[1][1], reverse=False)) #sort seqs by hammdist
        print("seqs have been sorted by distance value in descending order")
   #i=0
        slctd_seqs={}
        sampled_seqs=0
        print("parameters for selection have been initialized")
   # select top 20 percent in terms of hamm distance
        for seq in seqs_ordrd:
            if(sampled_seqs <= (selection_rate*totalseqs)):
       # THIS IS TAKING TOO LONG # UPDATE: NOT ANYMORE
                slctd_seqs[seq] = seqs_ordrd[seq] #add to selected pool
                sampled_seqs += slctd_seqs[seq][0] #add number of repeats
            else:
                break
       #i+=1 #increment index
   #seqs_ordrd.clear() #clear table of previous round
        print("sequence selection has been carried out")
        return slctd_seqs

    def stochasticSelection(self, seqLength, seqs, selection_rate, totalseqs):
        seqs_ordrd = OrderedDict(sorted(seqs.items(), key=lambda seqs: seqs[1][1], reverse=False)) #sort seqs by hammdist
        print("seqs have been sorted by distance value in descending order")
   #i=0
        slctd_seqs={}
        sampled_seqs=0
        print("parameters for selection have been initialized")
   # select top 20 percent in terms of hamm distance
        for seq in seqs_ordrd:
            if(sampled_seqs <= (selection_rate*totalseqs)):
                randHamm = random.randint(0, seqLength)
                if(seqs_ordrd[seq][1] <= randHamm):
       # THIS IS TAKING TOO LONG # UPDATE: NOT ANYMORE
                    slctd_seqs[seq] = seqs_ordrd[seq] #add to selected pool
                    sampled_seqs += slctd_seqs[seq][0] #add number of repeats
                else:
                    continue
            else:
                break
       #i+=1 #increment index
   #seqs_ordrd.clear() #clear table of previous round
        print("sequence selection has been carried out")
        return slctd_seqs

