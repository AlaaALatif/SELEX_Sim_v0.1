import numpy as np
import Aptamers
import Distance
import utils

# append path for ViennaRNA module
import RNA


D = Distance.Distance()
Apt = Aptamers.Aptamers()

# NEED TO CHANGE SAMPLING FOR SELECTION TO BE WEIGHTED BY COUNT OF EACH UNIQUE SEQ

# number of random samples to draw at a time
Nrsamples = 10**4


class Selection:
    def __init__(self, distname):
        self.distances = ("hamming", "basepair", "loop", "random")
        self.distname = distname
        if self.distname not in self.distances:
            print("Invalid argument for distance measure")
            raise
        self.distance = None
        if self.distname == "hamming":
            self.distance = D.hamming_func
        elif self.distname == "basepair":
            self.distance = D.bp_func
        elif self.distname == "loop":
            self.distance = D.loop_func
        elif self.distname == "random":
            self.distance = D.nodist_func

    # This function takes an empty selected pool, aptamer sequence structure and loop,
    # number of target binding sites, the alphabet set of the molecule, length,
    # total sequence number and stringency factor and returns full selected pool
    # their Lavenshtein distance
    # Input: dict(), str(), str(), str(), int(), str(), int(), int(), int()
    # Output: dict()
    def selectionProcess_loop_initial(self, slctdSeqs, aptSeq,
                                      aptStruct, aptLoop,
                                      selectionThreshold,
                                      alphabetSet, seqLength,
                                      totalSeqNum, stringency):
        selectedSeqs = 0
        print("Drawing sample batch")
        while(selectedSeqs < selectionThreshold):
            print("{}% completed".format(100.0*selectedSeqs/selectionThreshold))
            randIdxs = utils.randint(0, int(totalSeqNum-1), size=Nrsamples)
            randHamms = utils.randint(0, seqLength-stringency, size=Nrsamples)
            for i, randIdx in enumerate(randIdxs):
                randSeq = Apt.pseudoAptamerGenerator(randIdx, alphabetSet, seqLength)
                randSeqDist = D.loop_func(aptSeq, aptStruct, aptLoop, seqLength, randSeq)
                if(selectedSeqs == selectionThreshold):
                    return slctdSeqs
                elif(randSeqDist < randHamms[i]):
                    if(randIdx in slctdSeqs):
                        slctdSeqs[randIdx][0] += 1
                    else:
                        randSeqBias = D.bias_func(randSeq, seqLength)
                        slctdSeqs[randIdx] = np.array([1, randSeqDist, randSeqBias])
                    selectedSeqs += 1

    def selectionProcess_1D_initial(self, slctdSeqs, aptPool,
                                    selectionThreshold,
                                    alphabetSet, seqLength,
                                    totalSeqNum, stringency,
                                    distf=D.hamming_func):
        selectedSeqs = 0
        print("Drawing sample batch of {} sequences".format(selectionThreshold))
        while(selectedSeqs < selectionThreshold):
            randIdxs = utils.randint(0, int(totalSeqNum-1), size=Nrsamples)
            randHamms = utils.randint(0, seqLength-stringency, size=Nrsamples)
            for i, randIdx in enumerate(randIdxs):
                randSeq = Apt.pseudoAptamerGenerator(randIdx, alphabetSet, seqLength)
                # distance to optimal aptamer (stored in aptPool)
                randSeqDist = distf(aptPool, randSeq)
                if(selectedSeqs == selectionThreshold):
                    return slctdSeqs
                elif(randSeqDist < randHamms[i]):
                    if(randIdx in slctdSeqs):
                        slctdSeqs[randIdx][0] += 1
                    else:
                        randSeqBias = D.bias_func(randSeq, seqLength)
                        slctdSeqs[randIdx] = np.array([1, randSeqDist, randSeqBias])
                    selectedSeqs += 1
            print("{:6.2f}% completed".format(100.0*selectedSeqs/selectionThreshold), flush=True)
        return

    def stochasticSelection_initial(self, alphabetSet, seqLength,
                                    aptPool, selectionThreshold,
                                    totalSeqNum, samplingSize,
                                    outputFileNames, rnd, stringency):
        slctdSeqs = {}
        ref = aptPool
        if self.distname == "basepair":
            ref = RNA.fold(aptPool)[0]
            print("Optimum aptamer structure: {}".format(ref))
        print("Selection has started...", flush=True)
        # stochastic selection until threshold is met
        if self.distname == "loop":
            aptStruct = RNA.fold(aptPool)[0]
            print("Optimum aptamer structure: {}".format(aptStruct))
            aptLoop = utils.apt_loopFinder(aptPool, aptStruct, seqLength)
            self.selectionProcess_loop_initial(slctdSeqs, aptPool,
                                               aptStruct, aptLoop,
                                               selectionThreshold,
                                               alphabetSet, seqLength,
                                               totalSeqNum, stringency)
        else:
            self.selectionProcess_1D_initial(slctdSeqs, ref,
                                             selectionThreshold,
                                             alphabetSet, seqLength,
                                             totalSeqNum, stringency,
                                             distf=self.distance)
        print("sequence selection has been carried out")
        # sampling
        selectionDist = utils.rv_int(slctdSeqs, "selectionDist")
        print("sampling from initial round...")
        self.samplingProcess(alphabetSet, seqLength, slctdSeqs,
                             selectionDist, samplingSize,
                             outputFileNames, rnd)
        print("Sampling completed")
        return slctdSeqs

    def stochasticSelection(self, alphabetSet, seqLength,
                            seqPool, selectionThreshold,
                            uniqSeqNum, totalSeqNum, samplingSize,
                            outputFileNames, rnd, stringency):
        # initialize selected sequence pool
        print("seq length = "+str(seqLength))
        print("seq selection threshold = "+str(selectionThreshold))
        # compute sampling distribution for selection
        # using count of each unique seq
        selectionDist = utils.rv_int(seqPool, "selectionDist")
        print("Sampling has started...")
        self.samplingProcess(alphabetSet, seqLength, seqPool,
                             selectionDist, samplingSize,
                             outputFileNames, rnd)
        print("Sampling has completed")
        # reset all seq counts prior to selection
        for k in seqPool:
            seqPool[k][0] = 0
        # draw a bunch of random seqs
        self.selectionProcess(seqPool, selectionThreshold, selectionDist,
                              seqLength, stringency)
        # remove all seqs that haven't been selected
        for ki in [k for k, v in seqPool.items() if v[0] == 0]:
            del seqPool[ki]
        print("sequence selection has been carried out")
        return seqPool

    def samplingProcess(self, alphabetSet, seqLength,
                        seqPool, selectionDist, samplingSize,
                        outputFileNames, rnd):
        samps = dict()
        # draw random samples from distribution
        for seqIdx in selectionDist.rvs(size=samplingSize):
            if seqIdx in samps:
                samps[seqIdx] += 1
            else:
                samps[seqIdx] = 1
        sampleFileName = outputFileNames+"_samples_R{:03d}".format(rnd+1)
        # write to samples file
        with open(sampleFileName, 'w') as s:
            for seqIdx, N in samps.items():
                seq = Apt.pseudoAptamerGenerator(seqIdx, alphabetSet, seqLength)
                s.write(str(seq)+'\t'+str(int(seqPool[seqIdx][1]))+'\t'+str(N)+'\n')
        print("Sampling has completed")
        return

    # This function takes an empty selected pool, aptamer sequence structure and loop,
    # number of target binding sites, the alphabet set of the molecule, length,
    # total sequence number and stringency factor and returns full selected pool
    # Input: np.array(), int(), stats.obj(), int()
    # Output: np.array()
    def selectionProcess(self, seqPool, selectionThreshold,
                         selectionDist, seqLength,
                         stringency):
        selectedSeqs = 0
        # until all sites are occupied
        print("Drawing sample batch")
        while(selectedSeqs < selectionThreshold):
            # draw random sequences
            # warning: looping here causes large memory consumption
            for randIdx in selectionDist.rvs(size=Nrsamples):
                # carry out stochastic selection
                # draw random affinities
                if(int(seqPool[randIdx][1]) < utils.randint(0, seqLength-stringency)):
                    seqPool[randIdx][0] += 1
                    selectedSeqs += 1
                    if selectedSeqs % Nrsamples == 0:
                        print("{}% completed".format(100.0*selectedSeqs/selectionThreshold))
        return
