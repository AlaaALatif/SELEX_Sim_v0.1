import sys, time
from numpy import random
import numpy as np
import linecache
from itertools import izip, imap, islice, product
import operator
from collections import OrderedDict
from scipy import stats
import Aptamers, Distance, utils

D = Distance.Distance()
Apt = Aptamers.Aptamers()
#append path for ViennaRNA module
sys.path.append("/local/data/public/aaaa3/Simulations/ViennaRNA/lib/python2.7/site-packages/")
import RNA
from RNA import fold, bp_distance
## NEED TO CHANGE SAMPLING FOR SELECTION TO BE WEIGHTED BY COUNT OF EACH UNIQUE SEQ


class Selection: 
 
    def stochasticBasePairSelection_initial(self, alphabetSet, seqLength, aptPool, selectionThreshold, totalSeqNum):
        #initialize seqInfo matrix
        slctdSeqs = {} 
        selectedSeqs = 0
        aptStruct = fold(aptPool)[0]
        print("parameters for selection have been initialized")
        #stochastic selection until threshold is met
        while(selectedSeqs <= selectionThreshold):
            #draw a random hamming distance (i.e. affinity score)
            #randHammScore = random.randint(0, seqLength)
            randHammScore = random.randint(0, seqLength-3)
            #draw a random sequence index
            randSeqIdx = random.randint(0, int(totalSeqNum - 1))
            #generate the seq string from index
            randSeq = Apt.pseudoAptamerGenerator(randSeqIdx, alphabetSet, seqLength)
            #compute seq secondary structure
            randSeqStruct = fold(randSeq)[0]
            #compute the hamming distance (affinity score) for the seq
            randSeqDist = bp_distance(randSeqStruct, aptStruct)
            #stochastic selection protocol
            if(randSeqDist < randHammScore):
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

    def stochasticBasePairSelection(self, alphabetSet, seqLength, seqPool, selectionThreshold, uniqSeqNum, totalSeqNum):
        #initialize selected sequence pool
        slctdSeqs = {}
        selectedSeqs = 0
        print("seq length = "+str(seqLength))
        print("seq selection threshold = "+str(selectionThreshold))
        print("unique seq number = "+str(uniqSeqNum))
        print("sample distance = "+str(seqPool[seqPool.keys()[5]][1]))
        print("parameters for selection have been initialized")
        x = np.zeros((uniqSeqNum, 4))
        for i, seqIdx in enumerate(seqPool):
            x[i][0] = seqIdx
            x[i][1] = seqPool[seqIdx][0]
            x[i][2] = seqPool[seqIdx][1]
            x[i][3] = seqPool[seqIdx][2]
        print("Selection sample distribution being computed...")
        #compute sampling distribution for selection
        selectionDist = utils.rvd(x, totalSeqNum, "selectionDist")
        print("Selection sample distribution computed")
        for i, seqIdx in enumerate(seqPool):
            x[i][1] = 0
        while(selectedSeqs < selectionThreshold):
            rand_xIdx = selectionDist.rvs(size=1000000)
            print("Sample batch drawn...")
            for i, randIdx in enumerate(rand_xIdx):
                randHammScore = random.randint(0, seqLength-3)
                if( int(x[randIdx][2]) < randHammScore):
                    x[randIdx][1] += 1
                    selectedSeqs += 1
        del(rand_xIdx)
        x = x[x[:, 1] != 0]
        for seqInfo in x:
            #change it so that seqInfo are added as on np array, without setdefault
            slctdSeqs[int(seqInfo[0])] = seqInfo[1:]
        print("sequence selection has been carried out")
        return slctdSeqs


    def stochasticHammingSelection_initial(self, alphabetSet, seqLength, aptPool, selectionThreshold, totalSeqNum):
        #initialize seqInfo matrix
        slctdSeqs = {} 
        selectedSeqs = 0
        print("parameters for selection have been initialized")
        #stochastic selection until threshold is met
        while(selectedSeqs <= selectionThreshold):
            #draw a random hamming distance (i.e. affinity score)
            #randHammScore = random.randint(0, seqLength)
            randHammScore = random.randint(0, seqLength)
            #draw a random sequence index
            randSeqIdx = random.randint(0, int(totalSeqNum - 1))
            #generate the seq string from index
            randSeq = Apt.pseudoAptamerGenerator(randSeqIdx, alphabetSet, seqLength)
            #compute the hamming distance (affinity score) for the seq
            randSeqDist = D.hamming_func(randSeq, aptPool)
            #stochastic selection protocol
            if(randSeqDist < randHammScore):
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

    def stochasticHammingSelection(self, alphabetSet, seqLength, seqPool, selectionThreshold, uniqSeqNum, totalSeqNum):
        #initialize selected sequence pool
        slctdSeqs = {}
        selectedSeqs = 0
        print("seq length = "+str(seqLength))
        print("seq selection threshold = "+str(selectionThreshold))
        print("unique seq number = "+str(uniqSeqNum))
        print("eample distance = "+str(seqPool[seqPool.keys()[5]][1]))
        print("parameters for selection have been initialized")
        x = np.zeros((uniqSeqNum, 4))
        for i, seqIdx in enumerate(seqPool):
            x[i][0] = seqIdx
            x[i][1] = seqPool[seqIdx][0]
            x[i][2] = seqPool[seqIdx][1]
            x[i][3] = seqPool[seqIdx][2]
        print("Selection sample distribution being computed...")
        selectionDist = utils.rvd(x, totalSeqNum, "selectionDist")
        print("Selection sample distribution computed")
        for i, seqIdx in enumerate(seqPool):
            x[i][1] = 0
        while(selectedSeqs < selectionThreshold):
            rand_xIdx = selectionDist.rvs(size=1000000)
            print("Sample batch drawn...")
            for i, randIdx in enumerate(rand_xIdx):
                randHammScore = random.randint(0, seqLength)
                if( int(x[randIdx][2]) < randHammScore):
                    x[randIdx][1] += 1
                    selectedSeqs += 1
        del(rand_xIdx)
        x = x[x[:, 1] != 0]
        for seqInfo in x:
            #change it so that seqInfo are added as on np array, without setdefault
            slctdSeqs[int(seqInfo[0])] = seqInfo[1:]
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

