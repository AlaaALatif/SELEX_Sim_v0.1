import time, sys
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

#append ViennaRNA package to python path
sys.path.append("/local/data/public/aaaa3/Simulations/ViennaRNA/lib/python2.7/site-packages/")
import RNA
from RNA import fold, bp_distance

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

    #This method computes the probabilities of each possible number of mutations (1-seqLength)
    #These probabilities can used to approximated the fraction of sequences that will undergo
    #certain numbers of mutation, assuming sequence count is sufficiently large
    def get_mutation_probabilities_original(self):
        L = self.seqLength
        e = self.errorRate
        mutNumProbs = np.zeros(L+1)
        lamb = L*e
        for m in xrange(L+1):
            mutNumProbs[m] = lamb**(m)*exp(-lamb)/(fact(m))
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


    #This method computes the probabilities of each possible number of mutations (1-seqLength)
    #These probabilities can used to approximated the fraction of sequences that will undergo
    #certain numbers of mutation, assuming sequence count is sufficiently large
    def get_mutation_distribution_original(self):
        L = self.seqLength
        e = self.errorRate
        mutNumProbs = np.zeros(L+1)
        lamb = L*e
        for m in xrange(L+1):
            mutNumProbs[m] = lamb**(m)*exp(-lamb)/(fact(m))
        mut_m = np.arange(L+1)
        mutDist = stats.rv_discrete(name='Poisson-based mutation distribution', 
                                    values=(mut_m, mutNumProbs))
        return mutDist

# This method aims to carry out the mutations on the pool of sequences that are in 
# the given mutated pool. It also updates the counts of the wild-type sequence and their
# mutated variants to take into account pcr amplification during the process
    def generate_mutants_2D(self,
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
        # compute 2D structure of aptamer(s)
        aptamerSeqsStruct = fold(aptamerSeqs)[0]
        # for each seq in the mutation pool
        for seqIdx in mutatedPool:
            # grab probabilities to draw it after each pcr cycle
            cycleNumProbs = amplfdSeqs[seqIdx][3:]
            #print cycleNumProbs
            # compute a discrete distribution from probabilities
            cycleNumDist = convert_to_distribution(np.arange(pcrCycleNum), 
                                                    cycleNumProbs, 
                                                    'cycleNumDist')
            #print cycleNumDist.rvs(size=10)
            # for each mutation instance for the seq
            for mutNum, mutFreq in enumerate(mutatedPool[seqIdx]): 
                mutFreq = int(mutatedPool[seqIdx][mutNum])
               # if the mutation is carried out on less than 10,000 copies, draw random numbers...:(
                if mutFreq < 10000:
                    # draw random cycle numbers after which the sequences were drawn for mutation
                    cycleNums = cycleNumDist.rvs(size=mutFreq)
                    #generate the wild-type sequence string
                    wildTypeSeq = apt.pseudoAptamerGenerator(seqIdx, alphabetSet, seqLength)
                    #for each copy to be mutated
                    for mut in xrange(mutFreq):
                        wildTypeCount = 0
                        mutantCount = 0
                        #draw random positions on the seq to mutate
                        randPos = random.randint(1, seqLength+1, size=mutNum+1)
                        #draw a random nucleotide for each position
                        randNucs = random.randint(alphabetSize, size=mutNum+1)
                        mutatedSeq = wildTypeSeq
                        #for each position in seq, replace with random nucleotide
                        for posNum, pos in enumerate(randPos):
                            mutatedSeq = mutatedSeq[:(pos-1)] + alphabetSet[randNucs[posNum]] + \
                                         mutatedSeq[pos:]
                        #generate index of mutant based on string
                        mutatedSeqIdx = apt.pseudoAptamerIndexGenerator(mutatedSeq, 
                                                                            alphabetSet, 
                                                                            seqLength)
                        #if mutant already in amplified pool
                        if mutatedSeqIdx in amplfdSeqs:
                            #mutantNum = (1+pcrYld)**(pcrCycleNum - cycleNums[mut])
                            #for each pcr cycle after mutation has occured
                            for n in xrange(pcrCycleNum-cycleNums[mut]):
                                #compute amplified mutant count
                                mutantCount += int(binom(1, (pcrYld+amplfdSeqs[mutatedSeqIdx][2])))
                                #compute loss of count from wild-type
                                wildTypeCount += int(binom(1, (pcrYld+amplfdSeqs[seqIdx][2])))
                            #increment mutant seq count in amplified pool
                            amplfdSeqs[mutatedSeqIdx][0] += mutantCount
                            #decrement wild-type seq count in amplfied pool
                            amplfdSeqs[seqIdx][0] -= wildTypeCount
                        #if mutant not found in amplified pool
                        else:
                            #mutantNum = (1+pcrYld)**(pcrCycleNum - cycleNums[mut])
                            #add seq and its info to the amplified pool
                            initialMutCount = 1
                            mutatedSeqStruct = fold(mutatedSeq)[0]
                            mutDist = bp_distance(mutatedSeqStruct, aptamerSeqsStruct)
                            mutBias = d.bias_func(mutatedSeq, seqLength)
                            amplfdSeqs[mutatedSeqIdx] = np.array([initialMutCount, 
                                                                  mutDist, mutBias])
                            #for each pcr cycle after mutation has occured
                            for n in xrange(pcrCycleNum-cycleNums[mut]):
                                mutantCount += int(binom(1, (pcrYld+amplfdSeqs[mutatedSeqIdx][2])))
                                wildTypeCount += int(binom(1, (pcrYld+amplfdSeqs[seqIdx][2])))
                            #increment mutant seq count in amplified pool
                            amplfdSeqs[mutatedSeqIdx][0] += mutantCount
                            #decrement wild-type seq count in amplified pool
                            amplfdSeqs[seqIdx][0] -= wildTypeCount
                # if mutation carried out on more than 10,000 copies, avoid drawing random nums
                elif mutFreq > 10000:
                    # calculate fraction of mutants for each possible mutation
                    initialMutCount = int(0.333*mutFreq/seqLength)
                    # for each possible position that mutation can occur
                    for seqPos in xrange(seqLength):
                        # grab the sequence encoding array 
                        seqArray = apt.get_seqArray(seqIdx, alphabetSet, seqLength)
                        # if nucleotide in the position is adenine
                        if seqArray[seqPos] == 0:
                            # mutate adenine to cytosine
                            mutatedSeqIdx = seqIdx+(4**(seqLength-seqPos-1))
                            # if the mutated seq is already in amplified pool
                            if mutatedSeqIdx in amplfdSeqs:
                                # increment seq count 
                                amplfdSeqs[mutatedSeqIdx][0] += initialMutCount
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
                            else: # if not found in amplified pool
                                amplifiedMutCount = 0
                                for cycleNum, cycleNumProb in enumerate(cycleNumProbs):
                                    # compute expected number of mutant copies after amplification
                                    amplifiedMutCount += int(cycleNumProb* \
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
                                #compute 2D structure of mutant
                                mutatedSeqStruct = fold(mutatedSeq)[0]
                                #compute bp distance from aptamer(s)
                                mutDist = bp_distance(mutatedSeqStruct, aptamerSeqsStruct)
                                # compute bias score of seq
                                mutBias = d.bias_func(mutatedSeq, seqLength)
                                # add to amplified pool
                                amplfdSeqs[mutatedSeqIdx] = np.array([amplifiedMutCount, 
                                                                      mutDist, mutBias])
                            # mutate adenine to guanine
                            mutatedSeqIdx = seqIdx+(4**(seqLength-seqPos-1)*2)
                            # if the mutated seq is already in amplified pool
                            if mutatedSeqIdx in amplfdSeqs:
                                # increment seq count 
                                amplfdSeqs[mutatedSeqIdx][0] += initialMutCount
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
                            else: # if not found in amplified pool
                                amplifiedMutCount = 0
                                for cycleNum, cycleNumProb in enumerate(cycleNumProbs):
                                    # compute expected number of mutant copies after amplification
                                    amplifiedMutCount += int(cycleNumProb* \
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
                                #compute 2D structure of mutant
                                mutatedSeqStruct = fold(mutatedSeq)[0]
                                #compute bp distance from aptamer(s)
                                mutDist = bp_distance(mutatedSeqStruct, aptamerSeqsStruct)
                                # compute bias score of seq
                                mutBias = d.bias_func(mutatedSeq, seqLength)
                                # add to amplified pool
                                amplfdSeqs[mutatedSeqIdx] = np.array([amplifiedMutCount, 
                                                                      mutDist, mutBias])
                            # mutate adenine to thymine
                            mutatedSeqIdx = seqIdx+(4**(seqLength-seqPos-1)*3)
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
                                amplifiedMutCount = 0
                                for cycleNum, cycleNumProb in enumerate(cycleNumProbs):
                                    # compute expected number of mutant copies after amplification
                                    amplifiedMutCount += int(cycleNumProb* \
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
                                #compute 2D structure of mutant
                                mutatedSeqStruct = fold(mutatedSeq)[0]
                                #compute bp distance from aptamer(s)
                                mutDist = bp_distance(mutatedSeqStruct, aptamerSeqsStruct)
                                # add to amplified pool
                                amplfdSeqs[mutatedSeqIdx] = np.array([amplifiedMutCount, 
                                                                      mutDist, mutBias])
                        # if nucleotide in the position is cytosine
                        elif seqArray[seqPos] == 1:
                            # mutate cytosine to adenine
                            mutatedSeqIdx = seqIdx-(4**(seqLength-seqPos-1))
                            # if the mutated seq is already in amplified pool
                            if mutatedSeqIdx in amplfdSeqs:
                                # increment seq count 
                                amplfdSeqs[mutatedSeqIdx][0] += initialMutCount
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
                            else: # if not found in amplified pool
                                amplifiedMutCount = 0
                                for cycleNum, cycleNumProb in enumerate(cycleNumProbs):
                                    # compute expected number of mutant copies after amplification
                                    amplifiedMutCount += int(cycleNumProb* \
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
                                #compute 2D structure of mutant
                                mutatedSeqStruct = fold(mutatedSeq)[0]
                                #compute bp distance from aptamer(s)
                                mutDist = bp_distance(mutatedSeqStruct, aptamerSeqsStruct)
                                # compute bias score of seq
                                mutBias = d.bias_func(mutatedSeq, seqLength)
                                # add to amplified pool
                                amplfdSeqs[mutatedSeqIdx] = np.array([amplifiedMutCount, 
                                                                      mutDist, mutBias])
                            # mutate cytosine to guanine
                            mutatedSeqIdx = seqIdx+(4**(seqLength-seqPos-1))
                            # if the mutated seq is already in amplified pool
                            if mutatedSeqIdx in amplfdSeqs:
                                # increment seq count 
                                amplfdSeqs[mutatedSeqIdx][0] += initialMutCount
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
                            else: # if not found in amplified pool
                                amplifiedMutCount = 0
                                for cycleNum, cycleNumProb in enumerate(cycleNumProbs):
                                    # compute expected number of mutant copies after amplification
                                    amplifiedMutCount += int(cycleNumProb* \
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
                                #compute 2D structure of mutant
                                mutatedSeqStruct = fold(mutatedSeq)[0]
                                #compute bp distance from aptamer(s)
                                mutDist = bp_distance(mutatedSeqStruct, aptamerSeqsStruct)
                                # compute bias score of seq
                                mutBias = d.bias_func(mutatedSeq, seqLength)
                                # add to amplified pool
                                amplfdSeqs[mutatedSeqIdx] = np.array([amplifiedMutCount, 
                                                                      mutDist, mutBias])
                            # mutate cytosine to thymine
                            mutatedSeqIdx = seqIdx+(4**(seqLength-seqPos-1)*2)
                            # if the mutated seq is already in amplified pool
                            if mutatedSeqIdx in amplfdSeqs:
                                # increment seq count 
                                amplfdSeqs[mutatedSeqIdx][0] += initialMutCount
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
                            else: # if not found in amplified pool
                                amplifiedMutCount = 0
                                for cycleNum, cycleNumProb in enumerate(cycleNumProbs):
                                    # compute expected number of mutant copies after amplification
                                    amplifiedMutCount += int(cycleNumProb* \
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
                                #compute 2D structure of mutant
                                mutatedSeqStruct = fold(mutatedSeq)[0]
                                #compute bp distance from aptamer(s)
                                mutDist = bp_distance(mutatedSeqStruct, aptamerSeqsStruct)
                                # compute bias score of seq
                                mutBias = d.bias_func(mutatedSeq, seqLength)
                                # add to amplified pool
                                amplfdSeqs[mutatedSeqIdx] = np.array([amplifiedMutCount, 
                                                                      mutDist, mutBias])
                        # if nucleotide in the position is guanine
                        elif seqArray[seqPos] == 2:
                            # mutate guanine to adenine
                            mutatedSeqIdx = seqIdx-(4**(seqLength-seqPos-1)*2)
                            # if the mutated seq is already in amplified pool
                            if mutatedSeqIdx in amplfdSeqs:
                                # increment mutant seq count 
                                amplfdSeqs[mutatedSeqIdx][0] += initialMutCount
                                # decrement wild-type seq count
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
                            else: # if not found in amplified pool
                                amplifiedMutCount = 0
                                for cycleNum, cycleNumProb in enumerate(cycleNumProbs):
                                    # compute expected number of mutant copies after amplification
                                    amplifiedMutCount += int(cycleNumProb* \
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
                                #compute 2D structure of mutant
                                mutatedSeqStruct = fold(mutatedSeq)[0]
                                #compute bp distance from aptamer(s)
                                mutDist = bp_distance(mutatedSeqStruct, aptamerSeqsStruct)
                                # compute bias score of seq
                                mutBias = d.bias_func(mutatedSeq, seqLength)
                                # add to amplified pool
                                amplfdSeqs[mutatedSeqIdx] = np.array([amplifiedMutCount, 
                                                                      mutDist, mutBias])
                            # mutate guanine to cytosine
                            mutatedSeqIdx = seqIdx-(4**(seqLength-seqPos-1))
                            # if the mutated seq is already in amplified pool
                            if mutatedSeqIdx in amplfdSeqs:
                                # increment seq count 
                                amplfdSeqs[mutatedSeqIdx][0] += initialMutCount
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
                            else: # if not found in amplified pool
                                amplifiedMutCount = 0
                                for cycleNum, cycleNumProb in enumerate(cycleNumProbs):
                                    # compute expected number of mutant copies after amplification
                                    amplifiedMutCount += int(cycleNumProb* \
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
                                #compute 2D structure of mutant
                                mutatedSeqStruct = fold(mutatedSeq)[0]
                                #compute bp distance from aptamer(s)
                                mutDist = bp_distance(mutatedSeqStruct, aptamerSeqsStruct)
                                # compute bias score of seq
                                mutBias = d.bias_func(mutatedSeq, seqLength)
                                # add to amplified pool
                                amplfdSeqs[mutatedSeqIdx] = np.array([amplifiedMutCount, 
                                                                      mutDist, mutBias])
                            # mutate guanine to thymine
                            mutatedSeqIdx = seqIdx+(4**(seqLength-seqPos-1))
                            # if the mutated seq is already in amplified pool
                            if mutatedSeqIdx in amplfdSeqs:
                                # increment seq count 
                                amplfdSeqs[mutatedSeqIdx][0] += initialMutCount
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
                            else: # if not found in amplified pool
                                amplifiedMutCount = 0
                                for cycleNum, cycleNumProb in enumerate(cycleNumProbs):
                                    # compute expected number of mutant copies after amplification
                                    amplifiedMutCount += int(cycleNumProb* \
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
                                #compute 2D structure of mutant
                                mutatedSeqStruct = fold(mutatedSeq)[0]
                                #compute bp distance from aptamer(s)
                                mutDist = bp_distance(mutatedSeqStruct, aptamerSeqsStruct)
                                # compute bias score of seq
                                mutBias = d.bias_func(mutatedSeq, seqLength)
                                # add to amplified pool
                                amplfdSeqs[mutatedSeqIdx] = np.array([amplifiedMutCount, 
                                                                      mutDist, mutBias])
                        # if nucleotide in the position is thymine
                        elif seqArray[seqPos] == 3:
                            # mutate thymine to adenine
                            mutatedSeqIdx = seqIdx-(4**(seqLength-seqPos-1)*3)
                            # if the mutated seq is already in amplified pool
                            if mutatedSeqIdx in amplfdSeqs:
                                # increment seq count 
                                amplfdSeqs[mutatedSeqIdx][0] += initialMutCount
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
                            else: # if not found in amplified pool
                                amplifiedMutCount = 0
                                for cycleNum, cycleNumProb in enumerate(cycleNumProbs):
                                    # compute expected number of mutant copies after amplification
                                    amplifiedMutCount += int(cycleNumProb* \
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
                                #compute 2D structure of mutant
                                mutatedSeqStruct = fold(mutatedSeq)[0]
                                #compute bp distance from aptamer(s)
                                mutDist = bp_distance(mutatedSeqStruct, aptamerSeqsStruct)
                                # compute bias score of seq
                                mutBias = d.bias_func(mutatedSeq, seqLength)
                                # add to amplified pool
                                amplfdSeqs[mutatedSeqIdx] = np.array([amplifiedMutCount, 
                                                                      mutDist, mutBias])
                            # mutate thymine to cytosine
                            mutatedSeqIdx = seqIdx-(4**(seqLength-seqPos-1)*2)
                            # if the mutated seq is already in amplified pool
                            if mutatedSeqIdx in amplfdSeqs:
                                # increment seq count 
                                amplfdSeqs[mutatedSeqIdx][0] += initialMutCount
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
                            else: # if not found in amplified pool
                                amplifiedMutCount = 0
                                for cycleNum, cycleNumProb in enumerate(cycleNumProbs):
                                    # compute expected number of mutant copies after amplification
                                    amplifiedMutCount += int(cycleNumProb* \
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
                                #compute 2D structure of mutant
                                mutatedSeqStruct = fold(mutatedSeq)[0]
                                #compute bp distance from aptamer(s)
                                mutDist = bp_distance(mutatedSeqStruct, aptamerSeqsStruct)
                                # compute bias score of seq
                                mutBias = d.bias_func(mutatedSeq, seqLength)
                                # add to amplified pool
                                amplfdSeqs[mutatedSeqIdx] = np.array([amplifiedMutCount, 
                                                                      mutDist, mutBias])
                            # mutate thymine to guanine
                            mutatedSeqIdx = seqIdx-(4**(seqLength-seqPos-1))
                            # if the mutated seq is already in amplified pool
                            if mutatedSeqIdx in amplfdSeqs:
                                # increment seq count 
                                amplfdSeqs[mutatedSeqIdx][0] += initialMutCount
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
                            else: # if not found in amplified pool
                                amplifiedMutCount = 0
                                for cycleNum, cycleNumProb in enumerate(cycleNumProbs):
                                    # compute expected number of mutant copies after amplification
                                    amplifiedMutCount += int(cycleNumProb* \
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
                                #compute 2D structure of mutant
                                mutatedSeqStruct = fold(mutatedSeq)[0]
                                #compute bp distance from aptamer(s)
                                mutDist = bp_distance(mutatedSeqStruct, aptamerSeqsStruct)
                                # compute bias score of seq
                                mutBias = d.bias_func(mutatedSeq, seqLength)
                                # add to amplified pool
                                amplfdSeqs[mutatedSeqIdx] = np.array([amplifiedMutCount, 
                                                                      mutDist, mutBias])
                        else:
                            print("ERROR: seqPos integer does not correspond to any character(s) in the given alphabet set")
        print("Mutation has been carried out")
        return amplfdSeqs



# This method aims to carry out the mutations on the pool of sequences that are in 
# the given mutated pool. It also updates the counts of the wild-type sequence and their
# mutated variants to take into account pcr amplification during the process
    def generate_mutants_1D(self,
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
            cycleNumProbs = amplfdSeqs[seqIdx][3:]
            #print cycleNumProbs
            # compute a discrete distribution from probabilities
            cycleNumDist = convert_to_distribution(np.arange(pcrCycleNum), 
                                                    cycleNumProbs, 
                                                    'cycleNumDist')
            #print cycleNumDist.rvs(size=10)
            # for each mutation instance for the seq
            for mutNum, mutFreq in enumerate(mutatedPool[seqIdx]): 
                mutFreq = int(mutatedPool[seqIdx][mutNum])
               # if the mutation is carried out on less than 10,000 copies, draw random numbers...:(
                if mutFreq < 10000:
                    # draw random cycle numbers after which the sequences were drawn for mutation
                    cycleNums = cycleNumDist.rvs(size=mutFreq)
                    #generate the wild-type sequence string
                    wildTypeSeq = apt.pseudoAptamerGenerator(seqIdx, alphabetSet, seqLength)
                    #for each copy to be mutated
                    for mut in xrange(mutFreq):
                        wildTypeCount = 0
                        mutantCount = 0
                        #draw random positions on the seq to mutate
                        randPos = random.randint(1, seqLength+1, size=mutNum+1)
                        #draw a random nucleotide for each position
                        randNucs = random.randint(alphabetSize, size=mutNum+1)
                        mutatedSeq = wildTypeSeq
                        #for each position in seq, replace with random nucleotide
                        for posNum, pos in enumerate(randPos):
                            mutatedSeq = mutatedSeq[:(pos-1)] + alphabetSet[randNucs[posNum]] + \
                                         mutatedSeq[pos:]
                        #generate index of mutant based on string
                        mutatedSeqIdx = apt.pseudoAptamerIndexGenerator(mutatedSeq, 
                                                                            alphabetSet, 
                                                                            seqLength)
                        #if mutant already in amplified pool
                        if mutatedSeqIdx in amplfdSeqs:
                            #mutantNum = (1+pcrYld)**(pcrCycleNum - cycleNums[mut])
                            #for each pcr cycle after mutation has occured
                            for n in xrange(pcrCycleNum-cycleNums[mut]):
                                #compute amplified mutant count
                                mutantCount += int(binom(1, (pcrYld+amplfdSeqs[mutatedSeqIdx][2])))
                                #compute loss of count from wild-type
                                wildTypeCount += int(binom(1, (pcrYld+amplfdSeqs[seqIdx][2])))
                            #increment mutant seq count in amplified pool
                            amplfdSeqs[mutatedSeqIdx][0] += mutantCount
                            #decrement wild-type seq count in amplfied pool
                            amplfdSeqs[seqIdx][0] -= wildTypeCount
                        #if mutant not found in amplified pool
                        else:
                            #mutantNum = (1+pcrYld)**(pcrCycleNum - cycleNums[mut])
                            #add seq and its info to the amplified pool
                            initialMutCount = 1
                            mutDist = d.hamming_func(mutatedSeq, aptamerSeqs)
                            mutBias = d.bias_func(mutatedSeq, seqLength)
                            amplfdSeqs[mutatedSeqIdx] = np.array([initialMutCount, 
                                                                  mutDist, mutBias])
                            #for each pcr cycle after mutation has occured
                            for n in xrange(pcrCycleNum-cycleNums[mut]):
                                mutantCount += int(binom(1, (pcrYld+amplfdSeqs[mutatedSeqIdx][2])))
                                wildTypeCount += int(binom(1, (pcrYld+amplfdSeqs[seqIdx][2])))
                            #increment mutant seq count in amplified pool
                            amplfdSeqs[mutatedSeqIdx][0] += mutantCount
                            #decrement wild-type seq count in amplified pool
                            amplfdSeqs[seqIdx][0] -= wildTypeCount
                # if mutation carried out on more than 10,000 copies, avoid drawing random nums
                elif mutFreq > 10000:
                    # calculate fraction of mutants for each possible mutation
                    initialMutCount = int(0.333*mutFreq/seqLength)
                    # for each possible position that mutation can occur
                    for seqPos in xrange(seqLength):
                        # grab the sequence encoding array 
                        seqArray = apt.get_seqArray(seqIdx, alphabetSet, seqLength)
                        # if nucleotide in the position is adenine
                        if seqArray[seqPos] == 0:
                            # mutate adenine to cytosine
                            mutatedSeqIdx = seqIdx+(4**(seqLength-seqPos-1))
                            # if the mutated seq is already in amplified pool
                            if mutatedSeqIdx in amplfdSeqs:
                                # increment seq count 
                                amplfdSeqs[mutatedSeqIdx][0] += initialMutCount
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
                            else: # if not found in amplified pool
                                amplifiedMutCount = 0
                                for cycleNum, cycleNumProb in enumerate(cycleNumProbs):
                                    # compute expected number of mutant copies after amplification
                                    amplifiedMutCount += int(cycleNumProb* \
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
                                mutDist = d.hamming_func(mutatedSeq, aptamerSeqs)
                                # compute bias score of seq
                                mutBias = d.bias_func(mutatedSeq, seqLength)
                                # add to amplified pool
                                amplfdSeqs[mutatedSeqIdx] = np.array([amplifiedMutCount, 
                                                                      mutDist, mutBias])
                            # mutate adenine to guanine
                            mutatedSeqIdx = seqIdx+(4**(seqLength-seqPos-1)*2)
                            # if the mutated seq is already in amplified pool
                            if mutatedSeqIdx in amplfdSeqs:
                                # increment seq count 
                                amplfdSeqs[mutatedSeqIdx][0] += initialMutCount
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
                            else: # if not found in amplified pool
                                amplifiedMutCount = 0
                                for cycleNum, cycleNumProb in enumerate(cycleNumProbs):
                                    # compute expected number of mutant copies after amplification
                                    amplifiedMutCount += int(cycleNumProb* \
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
                                mutDist = d.hamming_func(mutatedSeq, aptamerSeqs)
                                # compute bias score of seq
                                mutBias = d.bias_func(mutatedSeq, seqLength)
                                # add to amplified pool
                                amplfdSeqs[mutatedSeqIdx] = np.array([amplifiedMutCount, 
                                                                      mutDist, mutBias])
                            # mutate adenine to thymine
                            mutatedSeqIdx = seqIdx+(4**(seqLength-seqPos-1)*3)
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
                                amplifiedMutCount = 0
                                for cycleNum, cycleNumProb in enumerate(cycleNumProbs):
                                    # compute expected number of mutant copies after amplification
                                    amplifiedMutCount += int(cycleNumProb* \
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
                                mutDist = d.hamming_func(mutatedSeq, aptamerSeqs)
                                # compute bias score of seq
                                mutBias = d.bias_func(mutatedSeq, seqLength)
                                # add to amplified pool
                                amplfdSeqs[mutatedSeqIdx] = np.array([amplifiedMutCount, 
                                                                      mutDist, mutBias])
                        # if nucleotide in the position is cytosine
                        elif seqArray[seqPos] == 1:
                            # mutate cytosine to adenine
                            mutatedSeqIdx = seqIdx-(4**(seqLength-seqPos-1))
                            # if the mutated seq is already in amplified pool
                            if mutatedSeqIdx in amplfdSeqs:
                                # increment seq count 
                                amplfdSeqs[mutatedSeqIdx][0] += initialMutCount
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
                            else: # if not found in amplified pool
                                amplifiedMutCount = 0
                                for cycleNum, cycleNumProb in enumerate(cycleNumProbs):
                                    # compute expected number of mutant copies after amplification
                                    amplifiedMutCount += int(cycleNumProb* \
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
                                mutDist = d.hamming_func(mutatedSeq, aptamerSeqs)
                                # compute bias score of seq
                                mutBias = d.bias_func(mutatedSeq, seqLength)
                                # add to amplified pool
                                amplfdSeqs[mutatedSeqIdx] = np.array([amplifiedMutCount, 
                                                                      mutDist, mutBias])
                            # mutate cytosine to guanine
                            mutatedSeqIdx = seqIdx+(4**(seqLength-seqPos-1))
                            # if the mutated seq is already in amplified pool
                            if mutatedSeqIdx in amplfdSeqs:
                                # increment seq count 
                                amplfdSeqs[mutatedSeqIdx][0] += initialMutCount
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
                            else: # if not found in amplified pool
                                amplifiedMutCount = 0
                                for cycleNum, cycleNumProb in enumerate(cycleNumProbs):
                                    # compute expected number of mutant copies after amplification
                                    amplifiedMutCount += int(cycleNumProb* \
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
                                mutDist = d.hamming_func(mutatedSeq, aptamerSeqs)
                                # compute bias score of seq
                                mutBias = d.bias_func(mutatedSeq, seqLength)
                                # add to amplified pool
                                amplfdSeqs[mutatedSeqIdx] = np.array([amplifiedMutCount, 
                                                                      mutDist, mutBias])
                            # mutate cytosine to thymine
                            mutatedSeqIdx = seqIdx+(4**(seqLength-seqPos-1)*2)
                            # if the mutated seq is already in amplified pool
                            if mutatedSeqIdx in amplfdSeqs:
                                # increment seq count 
                                amplfdSeqs[mutatedSeqIdx][0] += initialMutCount
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
                            else: # if not found in amplified pool
                                amplifiedMutCount = 0
                                for cycleNum, cycleNumProb in enumerate(cycleNumProbs):
                                    # compute expected number of mutant copies after amplification
                                    amplifiedMutCount += int(cycleNumProb* \
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
                                mutDist = d.hamming_func(mutatedSeq, aptamerSeqs)
                                # compute bias score of seq
                                mutBias = d.bias_func(mutatedSeq, seqLength)
                                # add to amplified pool
                                amplfdSeqs[mutatedSeqIdx] = np.array([amplifiedMutCount, 
                                                                      mutDist, mutBias])
                        # if nucleotide in the position is guanine
                        elif seqArray[seqPos] == 2:
                            # mutate guanine to adenine
                            mutatedSeqIdx = seqIdx-(4**(seqLength-seqPos-1)*2)
                            # if the mutated seq is already in amplified pool
                            if mutatedSeqIdx in amplfdSeqs:
                                # increment mutant seq count 
                                amplfdSeqs[mutatedSeqIdx][0] += initialMutCount
                                # decrement wild-type seq count
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
                            else: # if not found in amplified pool
                                amplifiedMutCount = 0
                                for cycleNum, cycleNumProb in enumerate(cycleNumProbs):
                                    # compute expected number of mutant copies after amplification
                                    amplifiedMutCount += int(cycleNumProb* \
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
                                mutDist = d.hamming_func(mutatedSeq, aptamerSeqs)
                                # compute bias score of seq
                                mutBias = d.bias_func(mutatedSeq, seqLength)
                                # add to amplified pool
                                amplfdSeqs[mutatedSeqIdx] = np.array([amplifiedMutCount, 
                                                                      mutDist, mutBias])
                            # mutate guanine to cytosine
                            mutatedSeqIdx = seqIdx-(4**(seqLength-seqPos-1))
                            # if the mutated seq is already in amplified pool
                            if mutatedSeqIdx in amplfdSeqs:
                                # increment seq count 
                                amplfdSeqs[mutatedSeqIdx][0] += initialMutCount
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
                            else: # if not found in amplified pool
                                amplifiedMutCount = 0
                                for cycleNum, cycleNumProb in enumerate(cycleNumProbs):
                                    # compute expected number of mutant copies after amplification
                                    amplifiedMutCount += int(cycleNumProb* \
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
                                mutDist = d.hamming_func(mutatedSeq, aptamerSeqs)
                                # compute bias score of seq
                                mutBias = d.bias_func(mutatedSeq, seqLength)
                                # add to amplified pool
                                amplfdSeqs[mutatedSeqIdx] = np.array([amplifiedMutCount, 
                                                                      mutDist, mutBias])
                            # mutate guanine to thymine
                            mutatedSeqIdx = seqIdx+(4**(seqLength-seqPos-1))
                            # if the mutated seq is already in amplified pool
                            if mutatedSeqIdx in amplfdSeqs:
                                # increment seq count 
                                amplfdSeqs[mutatedSeqIdx][0] += initialMutCount
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
                            else: # if not found in amplified pool
                                amplifiedMutCount = 0
                                for cycleNum, cycleNumProb in enumerate(cycleNumProbs):
                                    # compute expected number of mutant copies after amplification
                                    amplifiedMutCount += int(cycleNumProb* \
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
                                mutDist = d.hamming_func(mutatedSeq, aptamerSeqs)
                                # compute bias score of seq
                                mutBias = d.bias_func(mutatedSeq, seqLength)
                                # add to amplified pool
                                amplfdSeqs[mutatedSeqIdx] = np.array([amplifiedMutCount, 
                                                                      mutDist, mutBias])
                        # if nucleotide in the position is thymine
                        elif seqArray[seqPos] == 3:
                            # mutate thymine to adenine
                            mutatedSeqIdx = seqIdx-(4**(seqLength-seqPos-1)*3)
                            # if the mutated seq is already in amplified pool
                            if mutatedSeqIdx in amplfdSeqs:
                                # increment seq count 
                                amplfdSeqs[mutatedSeqIdx][0] += initialMutCount
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
                            else: # if not found in amplified pool
                                amplifiedMutCount = 0
                                for cycleNum, cycleNumProb in enumerate(cycleNumProbs):
                                    # compute expected number of mutant copies after amplification
                                    amplifiedMutCount += int(cycleNumProb* \
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
                                mutDist = d.hamming_func(mutatedSeq, aptamerSeqs)
                                # compute bias score of seq
                                mutBias = d.bias_func(mutatedSeq, seqLength)
                                # add to amplified pool
                                amplfdSeqs[mutatedSeqIdx] = np.array([amplifiedMutCount, 
                                                                      mutDist, mutBias])
                            # mutate thymine to cytosine
                            mutatedSeqIdx = seqIdx-(4**(seqLength-seqPos-1)*2)
                            # if the mutated seq is already in amplified pool
                            if mutatedSeqIdx in amplfdSeqs:
                                # increment seq count 
                                amplfdSeqs[mutatedSeqIdx][0] += initialMutCount
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
                            else: # if not found in amplified pool
                                amplifiedMutCount = 0
                                for cycleNum, cycleNumProb in enumerate(cycleNumProbs):
                                    # compute expected number of mutant copies after amplification
                                    amplifiedMutCount += int(cycleNumProb* \
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
                                mutDist = d.hamming_func(mutatedSeq, aptamerSeqs)
                                # compute bias score of seq
                                mutBias = d.bias_func(mutatedSeq, seqLength)
                                # add to amplified pool
                                amplfdSeqs[mutatedSeqIdx] = np.array([amplifiedMutCount, 
                                                                      mutDist, mutBias])
                            # mutate thymine to guanine
                            mutatedSeqIdx = seqIdx-(4**(seqLength-seqPos-1))
                            # if the mutated seq is already in amplified pool
                            if mutatedSeqIdx in amplfdSeqs:
                                # increment seq count 
                                amplfdSeqs[mutatedSeqIdx][0] += initialMutCount
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
                            else: # if not found in amplified pool
                                amplifiedMutCount = 0
                                for cycleNum, cycleNumProb in enumerate(cycleNumProbs):
                                    # compute expected number of mutant copies after amplification
                                    amplifiedMutCount += int(cycleNumProb* \
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
                                mutDist = d.hamming_func(mutatedSeq, aptamerSeqs)
                                # compute bias score of seq
                                mutBias = d.bias_func(mutatedSeq, seqLength)
                                # add to amplified pool
                                amplfdSeqs[mutatedSeqIdx] = np.array([amplifiedMutCount, 
                                                                      mutDist, mutBias])
                        else:
                            print("ERROR: seqPos integer does not correspond to any character(s) in the given alphabet set")
        print("Mutation has been carried out")
        return amplfdSeqs
