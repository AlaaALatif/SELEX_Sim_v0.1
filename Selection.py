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
    def stochasticHammingSelection(self, alphabetSet, seqLength, seqPool, aptPool, selectionThreshold, totalSeqNum, uniqSeqNum, roundNum):
        # if this is the first SELEX round
        if(roundNum == 0):
            #initialize seqInfo matrix
            slctdSeqs = {} 
            selectedSeqs = 0
            print("parameters for selection have been initialized")
            #stochastic selection until threshold is met
            while(selectedSeqs <= selectionThreshold):
                #draw a random hamming distance (i.e. affinity score)
                randHammScore = random.randint(0, seqLength)
                #draw a random sequence index
                randSeqIdx = random.randint(0, int(totalSeqNum - 1))
                #generate the seq string from index
                randSeq = Apt.pseudoAptamerGenerator(randSeqIdx, alphabetSet, seqLength)
                #compute the hamming distance (affinity score) for the seq
                randSeqDist = D.hamming_func(randSeq, aptPool)
                #stochastic selection protocol
                if(randSeqDist <= randHammScore):
                    #if seq already been selected
                    if(randSeqIdx in slctdSeqs):
                        #increment its count
                        slctdSeqs[randSeqIdx][0] += 1
                    #otherwise, add to the selected pool
                    else:
                        slctdSeqs.setdefault(randSeqIdx, []).append(1) #add initial count
                        slctdSeqs.setdefault(randSeqIdx, []).append(randSeqDist) #distance
                        slctdSeqs.setdefault(randSeqIdx, []).append(D.bias_func(randSeq, seqLength)) #bias
                    selectedSeqs += 1 #increment no. of selected seqs
            print("sequence selection has been carried out")
            return slctdSeqs
        # If this is not the 1st SELEX round
        else:
            #initialize selected sequence pool
            slctdSeqs = {}
            selectedSeqs = 0
            print("parameters for selection have been initialized")
            #stochastic selection until threshold is met
            while(selectedSeqs <= selectionThreshold):
                randHammScore = random.randint(0, seqLength)
                randPoolIdx = random.randint(0, int(uniqSeqNum-1))
                randSeqIdx = seqPool.keys()[randPoolIdx]
                #stochatic selection protocol
                if(seqPool[randSeqIdx][1] <= randHammScore):
                    #if seq is selected for the first time
                    if(randSeqIdx not in slctdSeqs):
                        slctdSeqs.setdefault(randSeqIdx, []).append(1) #add to selected pool
                        slctdSeqs.setdefault(randSeqIdx, []).append(seqPool[randSeqIdx][1]) #distance
                        slctdSeqs.setdefault(randSeqIdx, []).append(seqPool[randSeqIdx][2]) #bias
                    else:
                        slctdSeqs[randSeqIdx][0] += 1 #increment count
                    selectedSeqs += 1 #increment sampled no.
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

