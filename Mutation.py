import time
import linecache
from itertools import izip, imap
import operator
from collections import OrderedDict
from scipy import stats
import numpy as np
from numpy import random, exp, log
from numpy.random import binomial as binom
from math import factorial as fact
from sklearn.preprocessing import normalize
from utils import binomCoeff, convert_to_distribution
from Aptamers import Aptamers
from Distance import Distance
class Mutation(object):
    #constructor
    def __init__(self, seqLength=0, 
                 mutatedPool=None, 
                 amplfdSeqs=None,
                 aptamerSeqs=None, 
                 alphabetSet=None, errorRate=0, 
                 pcrCycleNum=0, pcrYld=0, 
                 seqPop=None):
        #initialize parameters
        self.seqLength = seqLength
        self.mutatedPool = mutatedPool
        self.amplfdSeqs = amplfdSeqs
        self.aptamerSeqs = aptamerSeqs
        self.alphabetSet = alphabetSet
        self.errorRate = errorRate
        self.pcrCycleNum = pcrCycleNum
        self.pcrYld = pcrYld
        self.seqPop = seqPop
        #add error handling for invalid param values 
    #This method computes the probability of drawing a seq after each pcr cycle
    def get_cycleNumber_probabilities(self, seqPop):
        #normalize each element so that they all sum to one (i.e. probability measures)
        cycleNumProbs = normalize(seqPop.reshape(1, -1), norm='l1')[0]
        return cycleNumProbs
    #This method computes the distribution of drawing seqs after the different pcr cycles
    def get_cycleNumber_distribution(self, seqPop):
        N = self.pcrCycleNum
        cycleNumProbs = normalize(seqPop.reshape(1, -1), norm='l1')[0]
        cycleVec = np.arange(N)
        # compute discrete distribution 
        cycleNumDist = stats.discrete.rv_discrete(name='cycleNumDist', 
                                                  values=(cycleVec, cycleNumProbs))
        return cycleNumDist
#TEST AREA
#mut = Mutation()
#mutDist = mut.get_distribution(20, 0.000006, 0.85, 15)

    #This method computes the probabilities of each possible number of mutations (1-seqLength)
    #These probabilities can used to approximated the fraction of sequences that will undergo
    #certain numbers of mutation, assuming sequence count is sufficiently large
    def get_mutation_probabilities(self):
        L = self.seqLength
        N = self.pcrCycleNum
        e = self.errorRate
        y = self.pcrYld
        mutNumProbs = np.zeros(L+1)
        for m in xrange(L):
            for n in xrange(1, N+1, 1):
                mutNumProbs[m+1] += exp(-n*e*L)* \
                                    (n*e*L)**(m+1)* \
                                    fact(N)*y**(n)/ \
                                    (fact(m+1)*fact(n)*fact(N-n)*(1+y)**(n))
            mutNumProbs[0] += mutNumProbs[m+1]
        mutNumProbs[0] = 1 - mutNumProbs[0]
        return mutNumProbs


    #This method computes the discrete distribution of number of mutations (1-seqLength)
    #This distribution can be used to draw random numbers of mutations
    #The method is relatively slow but can be used when sequence count is small
    def get_mutation_distribution(self):
        L = self.seqLength
        N = self.pcrCycleNum
        e = self.errorRate
        y = self.pcrYld
        prob_m = np.zeros(L+1)
        # for each mutation number
        for m in xrange(L):
            # compute the probability for it to occur
            for n in xrange(1, N+1, 1):
                prob_m[m+1] += exp(-n*e*L)* \
                               (n*e*L)**(m+1)* \
                               fact(N)*y**(n)/ \
                               (fact(m+1)*fact(n)*fact(N-n)*(1+y)**(n))
            prob_m[0] += prob_m[m+1]
        prob_m[0] = 1 - prob_m[0]
        # initialize vector containing each possible number of mutations (1-seqLength)
        mut_m = np.arange(L+1)
        # compute mutation number distribution
        mutDist = stats.rv_discrete(name='mutDist', values=(mut_m, prob_m))
        return mutDist

# This method aims to carry out the mutations on the pool of sequences that are in 
# the given mutated pool. It also updates the counts of the wild-type sequence and their
# mutated variants to take into account pcr amplification during the process
    def generate_mutants(self,
                          mutatedPool, amplfdSeqs, 
                          aptamerSeqs, alphabetSet):
        pcrCycleNum = self.pcrCycleNum
        pcrYld = self.pcrYld
        seqLength = self.seqLength
        # compute size of alphabet (i.e. 4 for DNA/RNA, 20 for peptides)
        alphabetSize = len(alphabetSet)
        # initialize aptamers class
        apt = Aptamers()
        # initialize distance class
        d = Distance()
        # for each seq in the mutation pool
        for seqIdx in mutatedPool:
            # grab probabilities to draw it after each pcr cycle
            cycleNumProbs = amplfdSeqs[seqIdx][3]
            #print cycleNumProbs
            # compute a discrete distribution from probabilities
            cycleNumDist = convert_to_distribution(np.arange(pcrCycleNum), 
                                                    cycleNumProbs, 
                                                    'cycleNumDist')
            #print cycleNumDist.rvs(size=10)
            # for each mutation instance for the seq
            for mutNum, mutFreq in enumerate(mutatedPool[seqIdx][1:]):
                mutNum += 1
                mutFreq = int(mutatedPool[seqIdx][mutNum])
                if mutFreq == 0:
                    break
               # if the mutation is carried out on less than 10,000 copies, draw random numbers...:(
                elif mutFreq < 10000:
                    # draw random cycle numbers after which the sequences were drawn for mutation
                    cycleNums = cycleNumDist.rvs(size=mutFreq)
                    #generate the wild-type sequence string
                    wildTypeSeq = apt.pseudoAptamerGenerator(seqIdx, alphabetSet, seqLength)
                    #for each copy to be mutated
                    for mut in xrange(mutFreq):
                        wildTypeCount = 0
                        mutantCount = 0
                        #draw random positions on the seq to mutate
                        randPos = random.randint(1, seqLength+1, size=mutNum)
                        #draw a random nucleotide for each position
                        randNucs = random.randint(alphabetSize, size=mutNum)
                        mutatedSeq = wildTypeSeq
                        #for each position in seq, replace with random nucleotide
                        for posNum, pos in enumerate(randPos):
                            mutatedSeq = mutatedSeq[:(pos-1)] + alphabetSet[randNucs[posNum]] + mutatedSeq[pos:]
                        #generate index of mutant based on string
                        mutatedSeqIdx = apt.pseudoAptamerIndexGenerator(mutatedSeq, 
                                                                            alphabetSet, 
                                                                            seqLength)
                        #if mutant already in amplified pool
                        if mutatedSeqIdx in amplfdSeqs:
                            #mutantNum = (1+pcrYld)**(pcrCycleNum - cycleNums[mut])
                            for n in xrange(pcrCycleNum-cycleNums[mut]):
                                mutantCount += int(binom(1, (pcrYld+amplfdSeqs[mutatedSeqIdx][2])))
                                wildTypeCount += int(binom(1, (pcrYld+amplfdSeqs[seqIdx][2])))
                            amplfdSeqs[mutatedSeqIdx][0] += mutantCount
                            amplfdSeqs[seqIdx][0] -= wildTypeCount
                        else:
                            #mutantNum = (1+pcrYld)**(pcrCycleNum - cycleNums[mut])
                            amplfdSeqs.setdefault(mutatedSeqIdx, []).append(1)
                            amplfdSeqs.setdefault(mutatedSeqIdx, []).append(d.hamming_func(mutatedSeq, aptamerSeqs))
                            amplfdSeqs.setdefault(mutatedSeqIdx, []).append(d.bias_func(mutatedSeq, seqLength))
                            for n in xrange(pcrCycleNum-cycleNums[mut]):
                                mutantCount += int(binom(1, (pcrYld+amplfdSeqs[mutatedSeqIdx][2])))
                                wildTypeCount += int(binom(1, (pcrYld+amplfdSeqs[seqIdx][2])))
                            amplfdSeqs[mutatedSeqIdx][0] += mutantCount
                            amplfdSeqs[seqIdx][0] -= wildTypeCount
                # if mutation carried out on more than 10,000 copies, avoid drawing random nums
                elif mutFreq > 10000:
                    print("mutFreq is "+str(mutFreq))
                    initialMutCount = int(0.25*mutFreq/seqLength)
                    # for each possible position that mutation can occur
                    for seqPos in xrange(seqLength):
                        # grab the sequence encoding array 
                        seqArray = apt.get_seqArray(seqIdx, alphabetSet, seqLength)
                        # proportion of seq copies with back mutation (i.e. no effect)
                        amplfdSeqs[seqIdx][0] += initialMutCount
                        # if nucleotide in position is adenine
                        if seqArray[seqPos] == 0:
                            # mutate adenine to cytosine
                            mutatedSeqIdx = seqIdx+(4**(seqLength-seqPos))
                            # if the mutated seq is already in amplified pool
                            if mutatedSeqIdx in amplfdSeqs:
                                # increment seq count 
                                amplfdSeqs[mutatedSeqIdx][0] += initialMutCount
                                amplfdSeqs[seqIdx][0] -= initialMutCount
                                for cycleNum, cycleNumProb in enumerate(cycleNumProbs):
                                    # compute expected number of mutant copies after amplification
                                    amplfdSeqs[mutatedSeqIdx][0] += int(cycleNumProb*initialMutCount* \
                                                                    (1+pcrYld)** \
                                                                    (pcrCycleNum-cycleNum))
                                    # compute expected decrease in no. of wild type seq
                                    amplfdSeqs[seqIdx][0] -= int(cycleNumProb*initialMutCount* \
                                                             (1+pcrYld)**(pcrCycleNum-cycleNum))
                            else: # if not found in amplified pool
                                # add mutated seq to amplified pool
                                amplfdSeqs.setdefault(mutatedSeqIdx, []).append(initialMutCount)
                                for cycleNum, cycleNumProb in enumerate(cycleNumProbs):
                                    # compute expected number of mutant copies after amplification
                                    amplfdSeqs[mutatedSeqIdx][0] += int(cycleNumProb* \
                                                                    initialMutCount* \
                                                                    (1+pcrYld)** \
                                                                    (pcrCycleNum-cycleNum))
                                    # compute expected decrease in no. of wild type seq
                                    amplfdSeqs[seqIdx][0] -= int(cycleNumProb*initialMutCount* \
                                                             (1+pcrYld)**(pcrCycleNum-cycleNum))
                                # generate seq string using its index
                                mutatedSeq = apt.pseudoAptamerGenerator(mutatedSeqIdx, 
                                                                        alphabetSet, 
                                                                        seqLength)
                                # compute hamming distance of seq
                                amplfdSeqs.setdefault(mutatedSeqIdx, []).append(
                                        d.hamming_func(mutatedSeq, aptamerSeqs))
                                # compute bias score of seq
                                amplfdSeqs.setdefault(mutatedSeqIdx, []).append(
                                        d.bias_func(mutatedSeq, seqLength))
                            # mutate adenine to guanine
                            mutatedSeqIdx = seqIdx+(4**(seqLength-seqPos)*2)
                            # if the mutated seq is already in amplified pool
                            if mutatedSeqIdx in amplfdSeqs:
                                amplfdSeqs[mutatedSeqIdx][0] += initialMutCount
                                amplfdSeqs[seqIdx][0] -= initialMutCount
                                for cycleNum, cycleNumProb in enumerate(cycleNumProbs):
                                    # compute expected number of mutant copies after amplification
                                    amplfdSeqs[mutatedSeqIdx][0] += int(cycleNumProb*initialMutCount* \
                                                                    (1+pcrYld)** \
                                                                    (pcrCycleNum-cycleNum))
                                    # compute expected decrease in no. of wild type seq
                                    amplfdSeqs[seqIdx][0] -= int(cycleNumProb*initialMutCount* \
                                                             (1+pcrYld)**(pcrCycleNum-cycleNum))
                            else: # if not found in amplified pool
                                # add mutated seq to amplified pool
                                amplfdSeqs.setdefault(mutatedSeqIdx, []).append(initialMutCount)
                                amplfdSeqs[seqIdx][0] -= initialMutCount
                                for cycleNum, cycleNumProb in enumerate(cycleNumProbs):
                                    # compute expected number of mutant copies after amplification
                                    amplfdSeqs[mutatedSeqIdx][0] += int(cycleNumProb* \
                                                                    initialMutCount* \
                                                                    (1+pcrYld)** \
                                                                    (pcrCycleNum-cycleNum))
                                    # compute expected decrease in no. of wild type seq
                                    amplfdSeqs[seqIdx][0] -= int(cycleNumProb*initialMutCount* \
                                                             (1+pcrYld)**(pcrCycleNum-cycleNum))
                                # generate seq string using its index
                                mutatedSeq = apt.pseudoAptamerGenerator(mutatedSeqIdx, 
                                                                        alphabetSet, 
                                                                        seqLength)
                                # compute hamming distance of seq
                                amplfdSeqs.setdefault(mutatedSeqIdx, []).append(
                                        d.hamming_func(mutatedSeq, aptamerSeqs))
                                # compute bias score of seq
                                amplfdSeqs.setdefault(mutatedSeqIdx, []).append(
                                        d.bias_func(mutatedSeq, seqLength))                        
                            # mutate adenine to thymine
                            mutatedSeqIdx = seqIdx+(4**(seqLength-seqPos)*3)
                            # if the mutated seq is already in amplified pool
                            if mutatedSeqIdx in amplfdSeqs:
                                amplfdSeqs[mutatedSeqIdx][0] += initialMutCount
                                amplfdSeqs[seqIdx][0] -= initialMutCount
                                for cycleNum, cycleNumProb in enumerate(cycleNumProb[0]):
                                    # compute expected number of mutant copies after amplification
                                    amplfdSeqs[mutatedSeqIdx][0] += int(cycleNumProb*initialMutCount* \
                                                                    (1+pcrYld)** \
                                                                    (pcrCycleNum-cycleNum))
                                    # compute expected decrease in no. of wild type seq
                                    amplfdSeqs[seqIdx][0] -= int(cycleNumProb*initialMutCount* \
                                                             (1+pcrYld)**(pcrCycleNum-cycleNum))
                            else: # if not found in amplified pool
                                # add mutated seq to amplified pool
                                amplfdSeqs.setdefault(mutatedSeqIdx, []).append(initialMutCount)
                                amplfdSeqs[seqIdx][0] -= initialMutCount
                                for cycleNum, cycleNumProb in enumerate(cycleNumProbs):
                                    # compute expected number of mutant copies after amplification
                                    amplfdSeqs[mutatedSeqIdx][0] += int(cycleNumProb* \
                                                                    initialMutCount* \
                                                                    (1+pcrYld)** \
                                                                    (pcrCycleNum-cycleNum))
                                    # compute expected decrease in no. of wild type seq
                                    amplfdSeqs[seqIdx][0] -= int(cycleNumProb*initialMutCount* \
                                                             (1+pcrYld)**(pcrCycleNum-cycleNum))
                                # generate seq string using its index
                                mutatedSeq = apt.pseudoAptamerGenerator(mutatedSeqIdx, 
                                                                        alphabetSet, 
                                                                        seqLength)
                                # compute hamming distance of seq
                                amplfdSeqs.setdefault(mutatedSeqIdx, []).append(
                                        d.hamming_func(mutatedSeq, aptamerSeqs))
                                # compute bias score of seq
                                amplfdSeqs.setdefault(mutatedSeqIdx, []).append(
                                        d.bias_func(mutatedSeq, seqLength)) 
                       # if nucleotide in position is cytosine
                        elif seqArray[seqPos] == 1:
                            # mutate cytosine to adenine
                            mutatedSeqIdx = seqIdx-(4**(seqLength-seqPos))
                            initialMutCount = int(0.25*mutFreq/seqLength)
                            # if the mutated seq is already in amplified pool
                            if mutatedSeqIdx in amplfdSeqs:
                                # increment seq count 
                                amplfdSeqs[mutatedSeqIdx][0] += initialMutCount
                                amplfdSeqs[seqIdx][0] -= initialMutCount
                                for cycleNum, cycleNumProb in enumerate(cycleNumProbs):
                                    # compute expected number of mutant copies after amplification
                                    amplfdSeqs[mutatedSeqIdx][0] += int(cycleNumProb*initialMutCount* \
                                                                    (1+pcrYld)** \
                                                                    (pcrCycleNum-cycleNum))
                                    # compute expected decrease in no. of wild type seq
                                    amplfdSeqs[seqIdx][0] -= int(cycleNumProb*initialMutCount* \
                                                             (1+pcrYld)**(pcrCycleNum-cycleNum))
                            else: # if not found in amplified pool
                                # add mutated seq to amplified pool
                                amplfdSeqs.setdefault(mutatedSeqIdx, []).append(initialMutCount)
                                for cycleNum, cycleNumProb in enumerate(cycleNumProbs):
                                    # compute expected number of mutant copies after amplification
                                    amplfdSeqs[mutatedSeqIdx][0] += int(cycleNumProb* \
                                                                    initialMutCount* \
                                                                    (1+pcrYld)** \
                                                                    (pcrCycleNum-cycleNum))
                                    # compute expected decrease in no. of wild type seq
                                    amplfdSeqs[seqIdx][0] -= int(cycleNumProb*initialMutCount* \
                                                             (1+pcrYld)**(pcrCycleNum-cycleNum))
                                # generate seq string using its index
                                mutatedSeq = apt.pseudoAptamerGenerator(mutatedSeqIdx, 
                                                                        alphabetSet, 
                                                                        seqLength)
                                # compute hamming distance of seq
                                amplfdSeqs.setdefault(mutatedSeqIdx, []).append(
                                        d.hamming_func(mutatedSeq, aptamerSeqs))
                                # compute bias score of seq
                                amplfdSeqs.setdefault(mutatedSeqIdx, []).append(
                                        d.bias_func(mutatedSeq, seqLength))
                            # mutate cytosine to guanine
                            mutatedSeqIdx = seqIdx+(4**(seqLength-seqPos))
                            # if the mutated seq is already in amplified pool
                            if mutatedSeqIdx in amplfdSeqs:
                                amplfdSeqs[mutatedSeqIdx][0] += initialMutCount
                                amplfdSeqs[seqIdx][0] -= initialMutCount
                                for cycleNum, cycleNumProb in enumerate(cycleNumProbs):
                                    # compute expected number of mutant copies after amplification
                                    amplfdSeqs[mutatedSeqIdx][0] += int(cycleNumProb*initialMutCount* \
                                                                    (1+pcrYld)** \
                                                                    (pcrCycleNum-cycleNum))
                                    # compute expected decrease in no. of wild type seq
                                    amplfdSeqs[seqIdx][0] -= int(cycleNumProb*initialMutCount* \
                                                             (1+pcrYld)**(pcrCycleNum-cycleNum))
                            else: # if not found in amplified pool
                                # add mutated seq to amplified pool
                                amplfdSeqs.setdefault(mutatedSeqIdx, []).append(initialMutCount)
                                amplfdSeqs[seqIdx][0] -= initialMutCount
                                for cycleNum, cycleNumProb in enumerate(cycleNumProbs):
                                    # compute expected number of mutant copies after amplification
                                    amplfdSeqs[mutatedSeqIdx][0] += int(cycleNumProb* \
                                                                    initialMutCount* \
                                                                    (1+pcrYld)** \
                                                                    (pcrCycleNum-cycleNum))
                                    # compute expected decrease in no. of wild type seq
                                    amplfdSeqs[seqIdx][0] -= int(cycleNumProb*initialMutCount* \
                                                             (1+pcrYld)**(pcrCycleNum-cycleNum))
                                # generate seq string using its index
                                mutatedSeq = apt.pseudoAptamerGenerator(mutatedSeqIdx, 
                                                                        alphabetSet, 
                                                                        seqLength)
                                # compute hamming distance of seq
                                amplfdSeqs.setdefault(mutatedSeqIdx, []).append(
                                        d.hamming_func(mutatedSeq, aptamerSeqs))
                                # compute bias score of seq
                                amplfdSeqs.setdefault(mutatedSeqIdx, []).append(
                                        d.bias_func(mutatedSeq, seqLength))                        
                            # mutate cytosine to thymine
                            mutatedSeqIdx = seqIdx+(4**(seqLength-seqPos)*2)
                            # if the mutated seq is already in amplified pool
                            if mutatedSeqIdx in amplfdSeqs:
                                amplfdSeqs[mutatedSeqIdx][0] += initialMutCount
                                amplfdSeqs[seqIdx][0] -= initialMutCount
                                for cycleNum, cycleNumProb in enumerate(cycleNumProb[0]):
                                    # compute expected number of mutant copies after amplification
                                    amplfdSeqs[mutatedSeqIdx][0] += int(cycleNumProb*initialMutCount* \
                                                                    (1+pcrYld)** \
                                                                    (pcrCycleNum-cycleNum))
                                    # compute expected decrease in no. of wild type seq
                                    amplfdSeqs[seqIdx][0] -= int(cycleNumProb*initialMutCount* \
                                                             (1+pcrYld)**(pcrCycleNum-cycleNum))
                            else: # if not found in amplified pool
                                # add mutated seq to amplified pool
                                amplfdSeqs.setdefault(mutatedSeqIdx, []).append(initialMutCount)
                                amplfdSeqs[seqIdx][0] -= initialMutCount
                                for cycleNum, cycleNumProb in enumerate(cycleNumProbs):
                                    # compute expected number of mutant copies after amplification
                                    amplfdSeqs[mutatedSeqIdx][0] += int(cycleNumProb* \
                                                                    initialMutCount* \
                                                                    (1+pcrYld)** \
                                                                    (pcrCycleNum-cycleNum))
                                    # compute expected decrease in no. of wild type seq
                                    amplfdSeqs[seqIdx][0] -= int(cycleNumProb*initialMutCount* \
                                                             (1+pcrYld)**(pcrCycleNum-cycleNum))
                                # generate seq string using its index
                                mutatedSeq = apt.pseudoAptamerGenerator(mutatedSeqIdx, 
                                                                        alphabetSet, 
                                                                        seqLength)
                                # compute hamming distance of seq
                                amplfdSeqs.setdefault(mutatedSeqIdx, []).append(
                                        d.hamming_func(mutatedSeq, aptamerSeqs))
                                # compute bias score of seq
                                amplfdSeqs.setdefault(mutatedSeqIdx, []).append(
                                        d.bias_func(mutatedSeq, seqLength)) 
                        # if nucleotide in position is guanine
                        elif seqArray[seqPos] == 2:
                            # mutate guanine to adenine
                            mutatedSeqIdx = seqIdx-(4**(seqLength-seqPos)*2)
                            initialMutCount = int(0.25*mutFreq/seqLength)
                            # if the mutated seq is already in amplified pool
                            if mutatedSeqIdx in amplfdSeqs:
                                # increment seq count 
                                amplfdSeqs[mutatedSeqIdx][0] += initialMutCount
                                amplfdSeqs[seqIdx][0] -= initialMutCount
                                for cycleNum, cycleNumProb in enumerate(cycleNumProbs):
                                    # compute expected number of mutant copies after amplification
                                    amplfdSeqs[mutatedSeqIdx][0] += int(cycleNumProb*initialMutCount* \
                                                                    (1+pcrYld)** \
                                                                    (pcrCycleNum-cycleNum))
                                    # compute expected decrease in no. of wild type seq
                                    amplfdSeqs[seqIdx][0] -= int(cycleNumProb*initialMutCount* \
                                                             (1+pcrYld)**(pcrCycleNum-cycleNum))
                            else: # if not found in amplified pool
                                # add mutated seq to amplified pool
                                amplfdSeqs.setdefault(mutatedSeqIdx, []).append(initialMutCount)
                                for cycleNum, cycleNumProb in enumerate(cycleNumProbs):
                                    # compute expected number of mutant copies after amplification
                                    amplfdSeqs[mutatedSeqIdx][0] += int(cycleNumProb* \
                                                                    initialMutCount* \
                                                                    (1+pcrYld)** \
                                                                    (pcrCycleNum-cycleNum))
                                    # compute expected decrease in no. of wild type seq
                                    amplfdSeqs[seqIdx][0] -= int(cycleNumProb*initialMutCount* \
                                                             (1+pcrYld)**(pcrCycleNum-cycleNum))
                                # generate seq string using its index
                                mutatedSeq = apt.pseudoAptamerGenerator(mutatedSeqIdx, 
                                                                        alphabetSet, 
                                                                        seqLength)
                                # compute hamming distance of seq
                                amplfdSeqs.setdefault(mutatedSeqIdx, []).append(
                                        d.hamming_func(mutatedSeq, aptamerSeqs))
                                # compute bias score of seq
                                amplfdSeqs.setdefault(mutatedSeqIdx, []).append(
                                        d.bias_func(mutatedSeq, seqLength))
                            # mutate guanine to cytosine
                            mutatedSeqIdx = seqIdx-(4**(seqLength-seqPos))
                            # if the mutated seq is already in amplified pool
                            if mutatedSeqIdx in amplfdSeqs:
                                amplfdSeqs[mutatedSeqIdx][0] += initialMutCount
                                amplfdSeqs[seqIdx][0] -= initialMutCount
                                for cycleNum, cycleNumProb in enumerate(cycleNumProbs):
                                    # compute expected number of mutant copies after amplification
                                    amplfdSeqs[mutatedSeqIdx][0] += int(cycleNumProb*initialMutCount* \
                                                                    (1+pcrYld)** \
                                                                    (pcrCycleNum-cycleNum))
                                    # compute expected decrease in no. of wild type seq
                                    amplfdSeqs[seqIdx][0] -= int(cycleNumProb*initialMutCount* \
                                                             (1+pcrYld)**(pcrCycleNum-cycleNum))
                            else: # if not found in amplified pool
                                # add mutated seq to amplified pool
                                amplfdSeqs.setdefault(mutatedSeqIdx, []).append(initialMutCount)
                                amplfdSeqs[seqIdx][0] -= initialMutCount
                                for cycleNum, cycleNumProb in enumerate(cycleNumProbs):
                                    # compute expected number of mutant copies after amplification
                                    amplfdSeqs[mutatedSeqIdx][0] += int(cycleNumProb* \
                                                                    initialMutCount* \
                                                                    (1+pcrYld)** \
                                                                    (pcrCycleNum-cycleNum))
                                    # compute expected decrease in no. of wild type seq
                                    amplfdSeqs[seqIdx][0] -= int(cycleNumProb*initialMutCount* \
                                                             (1+pcrYld)**(pcrCycleNum-cycleNum))
                                # generate seq string using its index
                                mutatedSeq = apt.pseudoAptamerGenerator(mutatedSeqIdx, 
                                                                        alphabetSet, 
                                                                        seqLength)
                                # compute hamming distance of seq
                                amplfdSeqs.setdefault(mutatedSeqIdx, []).append(
                                        d.hamming_func(mutatedSeq, aptamerSeqs))
                                # compute bias score of seq
                                amplfdSeqs.setdefault(mutatedSeqIdx, []).append(
                                        d.bias_func(mutatedSeq, seqLength))                        
                            # mutate guanine to thymine
                            mutatedSeqIdx = seqIdx+(4**(seqLength-seqPos))
                            # if the mutated seq is already in amplified pool
                            if mutatedSeqIdx in amplfdSeqs:
                                amplfdSeqs[mutatedSeqIdx][0] += initialMutCount
                                amplfdSeqs[seqIdx][0] -= initialMutCount
                                for cycleNum, cycleNumProb in enumerate(cycleNumProb[0]):
                                    # compute expected number of mutant copies after amplification
                                    amplfdSeqs[mutatedSeqIdx][0] += int(cycleNumProb*initialMutCount* \
                                                                    (1+pcrYld)** \
                                                                    (pcrCycleNum-cycleNum))
                                    # compute expected decrease in no. of wild type seq
                                    amplfdSeqs[seqIdx][0] -= int(cycleNumProb*initialMutCount* \
                                                             (1+pcrYld)**(pcrCycleNum-cycleNum))
                            else: # if not found in amplified pool
                                # add mutated seq to amplified pool
                                amplfdSeqs.setdefault(mutatedSeqIdx, []).append(initialMutCount)
                                amplfdSeqs[seqIdx][0] -= initialMutCount
                                for cycleNum, cycleNumProb in enumerate(cycleNumProbs):
                                    # compute expected number of mutant copies after amplification
                                    amplfdSeqs[mutatedSeqIdx][0] += int(cycleNumProb* \
                                                                    initialMutCount* \
                                                                    (1+pcrYld)** \
                                                                    (pcrCycleNum-cycleNum))
                                    # compute expected decrease in no. of wild type seq
                                    amplfdSeqs[seqIdx][0] -= int(cycleNumProb*initialMutCount* \
                                                             (1+pcrYld)**(pcrCycleNum-cycleNum))
                                # generate seq string using its index
                                mutatedSeq = apt.pseudoAptamerGenerator(mutatedSeqIdx, 
                                                                        alphabetSet, 
                                                                        seqLength)
                                # compute hamming distance of seq
                                amplfdSeqs.setdefault(mutatedSeqIdx, []).append(
                                        d.hamming_func(mutatedSeq, aptamerSeqs))
                                # compute bias score of seq
                                amplfdSeqs.setdefault(mutatedSeqIdx, []).append(
                                        d.bias_func(mutatedSeq, seqLength)) 
                        # if nucleotide in position is thymine
                        elif seqArray[seqPos] == 3:
                            # mutate thymine to adenine
                            mutatedSeqIdx = seqIdx-(4**(seqLength-seqPos)*3)
                            initialMutCount = int(0.25*mutFreq/seqLength)
                            # if the mutated seq is already in amplified pool
                            if mutatedSeqIdx in amplfdSeqs:
                                # increment seq count 
                                amplfdSeqs[mutatedSeqIdx][0] += initialMutCount
                                amplfdSeqs[seqIdx][0] -= initialMutCount
                                for cycleNum, cycleNumProb in enumerate(cycleNumProbs):
                                    # compute expected number of mutant copies after amplification
                                    amplfdSeqs[mutatedSeqIdx][0] += int(cycleNumProb*initialMutCount* \
                                                                    (1+pcrYld)** \
                                                                    (pcrCycleNum-cycleNum))
                                    # compute expected decrease in no. of wild type seq
                                    amplfdSeqs[seqIdx][0] -= int(cycleNumProb*initialMutCount* \
                                                             (1+pcrYld)**(pcrCycleNum-cycleNum))
                            else: # if not found in amplified pool
                                # add mutated seq to amplified pool
                                amplfdSeqs.setdefault(mutatedSeqIdx, []).append(initialMutCount)
                                for cycleNum, cycleNumProb in enumerate(cycleNumProbs):
                                    # compute expected number of mutant copies after amplification
                                    amplfdSeqs[mutatedSeqIdx][0] += int(cycleNumProb* \
                                                                    initialMutCount* \
                                                                    (1+pcrYld)** \
                                                                    (pcrCycleNum-cycleNum))
                                    # compute expected decrease in no. of wild type seq
                                    amplfdSeqs[seqIdx][0] -= int(cycleNumProb*initialMutCount* \
                                                             (1+pcrYld)**(pcrCycleNum-cycleNum))
                                # generate seq string using its index
                                mutatedSeq = apt.pseudoAptamerGenerator(mutatedSeqIdx, 
                                                                        alphabetSet, 
                                                                        seqLength)
                                # compute hamming distance of seq
                                amplfdSeqs.setdefault(mutatedSeqIdx, []).append(
                                        d.hamming_func(mutatedSeq, aptamerSeqs))
                                # compute bias score of seq
                                amplfdSeqs.setdefault(mutatedSeqIdx, []).append(
                                        d.bias_func(mutatedSeq, seqLength))
                            # mutate thymine to cytosine
                            mutatedSeqIdx = seqIdx-(4**(seqLength-seqPos)*2)
                            # if the mutated seq is already in amplified pool
                            if mutatedSeqIdx in amplfdSeqs:
                                amplfdSeqs[mutatedSeqIdx][0] += initialMutCount
                                amplfdSeqs[seqIdx][0] -= initialMutCount
                                for cycleNum, cycleNumProb in enumerate(cycleNumProbs):
                                    # compute expected number of mutant copies after amplification
                                    amplfdSeqs[mutatedSeqIdx][0] += int(cycleNumProb*initialMutCount* \
                                                                    (1+pcrYld)** \
                                                                    (pcrCycleNum-cycleNum))
                                    # compute expected decrease in no. of wild type seq
                                    amplfdSeqs[seqIdx][0] -= int(cycleNumProb*initialMutCount* \
                                                             (1+pcrYld)**(pcrCycleNum-cycleNum))
                            else: # if not found in amplified pool
                                # add mutated seq to amplified pool
                                amplfdSeqs.setdefault(mutatedSeqIdx, []).append(initialMutCount)
                                amplfdSeqs[seqIdx][0] -= initialMutCount
                                for cycleNum, cycleNumProb in enumerate(cycleNumProbs):
                                    # compute expected number of mutant copies after amplification
                                    amplfdSeqs[mutatedSeqIdx][0] += int(cycleNumProb* \
                                                                    initialMutCount* \
                                                                    (1+pcrYld)** \
                                                                    (pcrCycleNum-cycleNum))
                                    # compute expected decrease in no. of wild type seq
                                    amplfdSeqs[seqIdx][0] -= int(cycleNumProb*initialMutCount* \
                                                             (1+pcrYld)**(pcrCycleNum-cycleNum))
                                # generate seq string using its index
                                mutatedSeq = apt.pseudoAptamerGenerator(mutatedSeqIdx, 
                                                                        alphabetSet, 
                                                                        seqLength)
                                # compute hamming distance of seq
                                amplfdSeqs.setdefault(mutatedSeqIdx, []).append(
                                        d.hamming_func(mutatedSeq, aptamerSeqs))
                                # compute bias score of seq
                                amplfdSeqs.setdefault(mutatedSeqIdx, []).append(
                                        d.bias_func(mutatedSeq, seqLength))                        
                            # mutate thymine to guanine
                            mutatedSeqIdx = seqIdx-(4**(seqLength-seqPos))
                            # if the mutated seq is already in amplified pool
                            if mutatedSeqIdx in amplfdSeqs:
                                amplfdSeqs[mutatedSeqIdx][0] += initialMutCount
                                amplfdSeqs[seqIdx][0] -= initialMutCount
                                for cycleNum, cycleNumProb in enumerate(cycleNumProb[0]):
                                    # compute expected number of mutant copies after amplification
                                    amplfdSeqs[mutatedSeqIdx][0] += int(cycleNumProb*initialMutCount* \
                                                                    (1+pcrYld)** \
                                                                    (pcrCycleNum-cycleNum))
                                    # compute expected decrease in no. of wild type seq
                                    amplfdSeqs[seqIdx][0] -= int(cycleNumProb*initialMutCount* \
                                                             (1+pcrYld)**(pcrCycleNum-cycleNum))
                            else: # if not found in amplified pool
                                # add mutated seq to amplified pool
                                amplfdSeqs.setdefault(mutatedSeqIdx, []).append(initialMutCount)
                                amplfdSeqs[seqIdx][0] -= initialMutCount
                                for cycleNum, cycleNumProb in enumerate(cycleNumProbs):
                                    # compute expected number of mutant copies after amplification
                                    amplfdSeqs[mutatedSeqIdx][0] += int(cycleNumProb* \
                                                                    initialMutCount* \
                                                                    (1+pcrYld)** \
                                                                    (pcrCycleNum-cycleNum))
                                    # compute expected decrease in no. of wild type seq
                                    amplfdSeqs[seqIdx][0] -= int(cycleNumProb*initialMutCount* \
                                                             (1+pcrYld)**(pcrCycleNum-cycleNum))
                                # generate seq string using its index
                                mutatedSeq = apt.pseudoAptamerGenerator(mutatedSeqIdx, 
                                                                        alphabetSet, 
                                                                        seqLength)
                                # compute hamming distance of seq
                                amplfdSeqs.setdefault(mutatedSeqIdx, []).append(
                                        d.hamming_func(mutatedSeq, aptamerSeqs))
                                # compute bias score of seq
                                amplfdSeqs.setdefault(mutatedSeqIdx, []).append(
                                        d.bias_func(mutatedSeq, seqLength))
                        else:
                            print("ERROR: seqPos integer does not correspond to any character(s) in the given alphabet set")
        print("Mutation has been carried out")
        return amplfdSeqs
##### NEED TO INCLUDE GENERATION NUMBER DISTRIBUTION; ACCOUNT FOR PARTIAL AMPLIFICATIONS OF WILD TYPE AND MUTANTS
# randomly mutates each seq copy in mutatedPool using the corresponding
# random positions. random mutations can be circular (i.e. A to A)
# mutated seqs are added to the amplfied pool of wild-type seqs
    def mutantGeneration(seqLength, mutatedPool, amplfdSeqs, aptamerSeqs, alphabetSet):
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

   
    
