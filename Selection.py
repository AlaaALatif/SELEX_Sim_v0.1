import sys
from numpy import random
import numpy as np
import Aptamers
import Distance
import utils

D = Distance.Distance()
Apt = Aptamers.Aptamers()
#append path for ViennaRNA module
sys.path.append("/local/data/public/aaaa3/Simulations/ViennaRNA/lib/python2.7/site-packages/")
import RNA
## NEED TO CHANGE SAMPLING FOR SELECTION TO BE WEIGHTED BY COUNT OF EACH UNIQUE SEQ


class Selection:

    def stochasticBasePairSelection_initial(self, alphabetSet, seqLength, aptPool, selectionThreshold, totalSeqNum, samplingSize, outputFileNames, rnd):
        #sampling
        print("sampling from initial library...")
        randomSamples = random.randint(0, int(totalSeqNum-1), size=samplingSize)
        sampleFileName = outputFileNames+"_samples_R{:03d}".format(rnd+1)
        with open(sampleFileName, 'w') as s:
            for seqIdx in randomSamples:
                seq = Apt.pseudoAptamerGenerator(seqIdx, alphabetSet, seqLength)
                s.write(seq+'\n')
        print("Sampling completed")
        #initialize seqInfo matrix
        slctdSeqs = {} 
        selectedSeqs = 0
        aptStruct = RNA.fold(aptPool)[0]
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
            randSeqStruct = RNA.fold(randSeq)[0]
            #compute the hamming distance (affinity score) for the seq
            randSeqDist = RNA.bp_distance(randSeqStruct, aptStruct)
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

    def stochasticBasePairSelection(self, alphabetSet, seqLength, seqPool, selectionThreshold, uniqSeqNum, totalSeqNum, samplingSize, outputFileNames, rnd):
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
        print("Sampling has started...")
        randSamples = selectionDist.rvs(size=samplingSize)
        sampleFileName = outputFileNames+"_samples_R{:03d}".format(rnd+1)
        with open(sampleFileName, 'w') as s:
            for seqIdx in randSamples:
                seq = Apt.pseudoAptamerGenerator(x[seqIdx][0], alphabetSet, seqLength)
                s.write(str(seq)+'\t'+str(int(x[seqIdx][1]))+'\n')
        print("Sampling has completed")
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


    def stochasticHammingSelection_initial(self, alphabetSet, seqLength, aptPool, selectionThreshold, totalSeqNum, samplingSize, outputFileNames, rnd):
        #sampling
        print("sampling from initial library...")
        randomSamples = random.randint(0, int(totalSeqNum-1), size=samplingSize)
        sampleFileName = outputFileNames+"_samples_R{:03d}".format(rnd+1)
        with open(sampleFileName, 'w') as s:
            for seqIdx in randomSamples:
                seq = Apt.pseudoAptamerGenerator(seqIdx, alphabetSet, seqLength)
                s.write(seq+'\n')
        print("Sampling completed")
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

    def stochasticHammingSelection(self, alphabetSet, seqLength, seqPool, selectionThreshold, uniqSeqNum, totalSeqNum, samplingSize, outputFileNames, rnd):
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
        selectionDist = utils.rvd(x, totalSeqNum, "selectionDist")
        print("Selection sample distribution computed")
        print("Sampling has started...")
        randSamples = selectionDist.rvs(size=samplingSize)
        sampleFileName = outputFileNames+"_samples_R{:03d}".format(rnd+1)
        with open(sampleFileName, 'w') as s:
            for seqIdx in randSamples:
                seq = Apt.pseudoAptamerGenerator(x[seqIdx][0], alphabetSet, seqLength)
                s.write(str(seq)+'\t'+str(int(x[seqIdx][1]))+'\n')
        print("Sampling has completed")
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
