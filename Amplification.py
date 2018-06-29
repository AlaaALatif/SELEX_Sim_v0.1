import numpy as np
from numpy.random import binomial as binom
from numpy.random import poisson

from Mutation import Mutation


# Initiate class
class Amplification:
    def randomPCR_with_ErrorsAndBias(self, slctdSeqs,
                                     seqLength, pcrCycleNum,
                                     pcrYld, errorRate,
                                     aptamerSeqs, alphabetSet, distance):
        return self.randomPCR_with_ErrorsAndBias_FASTv2(slctdSeqs,
                                                        seqLength, pcrCycleNum,
                                                        pcrYld, errorRate,
                                                        aptamerSeqs, alphabetSet, distance)

    def randomPCR_with_ErrorsAndBias_FASTv2(self, slctdSeqs,
                                            seqLength, pcrCycleNum,
                                            pcrYld, errorRate,
                                            aptamerSeqs, alphabetSet, distance):
        # initialize Mutation object from class
        mut = Mutation(seqLength=seqLength, errorRate=errorRate,
                       pcrCycleNum=pcrCycleNum, pcrYld=pcrYld)
        # count number of seqs in selected pool
        totalseqs = 0
        uniqSeqs = 0
        # compute total seq num, unique seq num, and transfer info to x
        for i, seqIdx in enumerate(slctdSeqs):
            uniqSeqs += 1
            totalseqs += int(slctdSeqs[seqIdx][0])
        for seqIdx in slctdSeqs:
            slctdSeqs[seqIdx].resize(pcrCycleNum+3)
        print("number of unique seqs in selected pool prior to amplification: "+str(uniqSeqs))
        print("number of seqs in selected pool prior to amplification: "+str(totalseqs))
        # calculate probabilities of different possible mutation numbers
        mutNumProbs = mut.get_mutation_probabilities_original()
        # # compute a discrete distribution of mutation numbers
        # mutDist = mut.get_mutation_distribution_original()
        print("Discrete Mutation Distribution has been computed")
    # PCR Amplification
        totalseqs = 0
        # initialize dictionary to keep info on seqs to be mutated
        mutatedPool = {}
        # keep track of sequence count after each pcr cycle (except last one)
        seqPop = np.zeros(pcrCycleNum)
        # compute cycle number probabilities for this seq
        cycleNumProbs = np.zeros(pcrCycleNum)
        print("Amplification has started...")
        # for each sequence in the selected pool
        for i, seqIdx in enumerate(slctdSeqs):
            mutatedPool[seqIdx] = np.zeros(seqLength)
            sn = slctdSeqs[seqIdx][0]
            # random PCR with bias using brute force
            for n in range(pcrCycleNum):
                # sequence count after n cycles
                seqPop[n] = sn
                # amplify count using initial count, polymerase yield, and bias score
                sn += int(binom(sn, min(0.9999, pcrYld+slctdSeqs[seqIdx][2])))
            slctdSeqs[seqIdx][0] = sn
            # compute cycle number probabilities
            for s, seqNum in enumerate(seqPop):
                cycleNumProbs[s] = seqNum/np.sum(seqPop)
            # transfer info to x
            for j, cycleNumProb in enumerate(cycleNumProbs):
                slctdSeqs[seqIdx][j+3] = cycleNumProb
            # update total num of seqs
            totalseqs += slctdSeqs[seqIdx][0]
            # tranfer seq index to matrix y
            # if accumulated seq count is greater than 10,000
            if np.sum(seqPop) > 10000:
                # for each possible number of mutations in any seq copy (1-seqLength)
                for mutNum in range(seqLength):
                    # approximate the proportion of copies that will be mutated using
                    # corresponding probability p(M=mutNum)
                    mutatedPool[seqIdx][mutNum] = mutNumProbs[mutNum+1]*np.sum(seqPop)
            # if seq count is less than 10,000
            else:
                # draw random mutNum from the mutation distribution for each seq copy
                muts = poisson(errorRate*seqLength, int(np.sum(seqPop)))  # SLOW STEP
                # remove all drawn numbers equal to zero
                muts = muts[muts != 0]
                # for each non-zero mutation number
                for mutNum in muts:
                    # increment copy number to be mutated
                    mutatedPool[seqIdx][mutNum] += 1
        print("Amplification carried out")
        print("Sequence selection for mutation has started...")
        # remove all mutation numbers with zero copies to be mutated
        for kd in [k for k, v in mutatedPool.items() if v[0] == 0]:
            del mutatedPool[kd]
        for k, mi in mutatedPool.items():
            mutatedPool[k] = mi[mi != 0]
        print("Mutation selection has been carried out")
        print("Mutant generation has started...")
        print("Mutating {} sequences...".format(len(mutatedPool)))
        # generate mutants and add to the amplfied sequence pool
        mut.generate_mutants(mutatedPool=mutatedPool,
                             amplfdSeqs=slctdSeqs,
                             aptamerSeqs=aptamerSeqs,
                             alphabetSet=alphabetSet,
                             distname=distance)
        return slctdSeqs
