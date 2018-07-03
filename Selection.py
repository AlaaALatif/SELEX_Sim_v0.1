import numpy as np
import Distance
import utils

# append path for ViennaRNA module
import RNA


D = Distance.Distance()

# NEED TO CHANGE SAMPLING FOR SELECTION TO BE WEIGHTED BY COUNT OF EACH UNIQUE SEQ

# number of random samples to draw at a time
Nrsamples = 10**4


class Selection:
    def __init__(self, distname, selectionThreshold, initialSize, samplingSize, stringency):
        self.distances = ("hamming", "basepair", "loop", "random")
        self.distname = distname
        self.selectionThreshold = selectionThreshold
        self.initialSize = initialSize
        self.samplingSize = samplingSize
        self.stringency = stringency
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
                                      apt,
                                      totalSeqNum):
        selectedSeqs = 0
        print("Drawing sample batch of {} sequences".format(self.initialSize))
        while(selectedSeqs < self.initialSize):
            randIdxs = utils.randint(0, int(totalSeqNum-1), size=Nrsamples)
            randHamms = utils.randint(0, apt.seqLength-self.stringency, size=Nrsamples)
            for i, randIdx in enumerate(randIdxs):
                randSeq = apt.pseudoAptamerGenerator(randIdx)
                randSeqDist = D.loop_func(aptSeq, aptStruct, aptLoop, apt.seqLength, randSeq)
                if(selectedSeqs == self.initialSize):
                    return slctdSeqs
                elif(randSeqDist < randHamms[i]):
                    if(randIdx in slctdSeqs):
                        slctdSeqs[randIdx][0] += 1
                    else:
                        randSeqBias = D.bias_func(randSeq, apt.seqLength)
                        slctdSeqs[randIdx] = np.array([1, randSeqDist, randSeqBias])
                    selectedSeqs += 1
            print("{:6.2f}% completed".format(100.0*selectedSeqs/self.initialSize), flush=True)
        return

    def selectionProcess_1D_initial(self, slctdSeqs, aptPool,
                                    apt,
                                    totalSeqNum,
                                    distf=D.hamming_func):
        selectedSeqs = 0
        print("Drawing sample batch of {} sequences".format(self.initialSize))
        while(selectedSeqs < self.initialSize):
            randIdxs = utils.randint(0, int(totalSeqNum-1), size=Nrsamples)
            randHamms = utils.randint(0, apt.seqLength-self.stringency, size=Nrsamples)
            for i, randIdx in enumerate(randIdxs):
                randSeq = apt.pseudoAptamerGenerator(randIdx)
                # distance to optimal aptamer (stored in aptPool)
                randSeqDist = distf(aptPool, randSeq)
                if(selectedSeqs == self.initialSize):
                    return slctdSeqs
                elif(randSeqDist < randHamms[i]):
                    if(randIdx in slctdSeqs):
                        slctdSeqs[randIdx][0] += 1
                    else:
                        randSeqBias = D.bias_func(randSeq, apt.seqLength)
                        slctdSeqs[randIdx] = np.array([1, randSeqDist, randSeqBias])
                    selectedSeqs += 1
            print("{:6.2f}% completed".format(100.0*selectedSeqs/self.initialSize), flush=True)
        return

    def stochasticSelection_initial(self, apt, aptPool,
                                    totalSeqNum,
                                    outputFileNames, rnd):
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
            aptLoop = utils.apt_loopFinder(aptPool, aptStruct, apt.seqLength)
            self.selectionProcess_loop_initial(slctdSeqs, aptPool,
                                               aptStruct, aptLoop,
                                               apt, totalSeqNum)
        else:
            self.selectionProcess_1D_initial(slctdSeqs, ref,
                                             apt, totalSeqNum,
                                             distf=self.distance)
        print("sequence selection has been carried out")
        # sampling
        selectionDist = utils.rv_int(slctdSeqs, "selectionDist")
        print("sampling from initial round...")
        self.samplingProcess(apt, slctdSeqs, selectionDist, self.samplingSize,
                             outputFileNames, rnd)
        print("Sampling completed")
        return slctdSeqs

    def stochasticSelection(self, apt, seqPool,
                            outputFileNames, rnd):
        # initialize selected sequence pool
        print("seq selection threshold = "+str(self.selectionThreshold))
        # compute sampling distribution for selection
        # using count of each unique seq
        selectionDist = utils.rv_int(seqPool, "selectionDist")
        print("Sampling has started...")
        self.samplingProcess(apt, seqPool, selectionDist, self.samplingSize,
                             outputFileNames, rnd)
        print("Sampling has completed")
        # reset all seq counts prior to selection
        for k in seqPool:
            seqPool[k][0] = 0
        # draw a bunch of random seqs
        self.selectionProcess(seqPool, selectionDist, apt.seqLength)
        # remove all seqs that haven't been selected
        for ki in [k for k, v in seqPool.items() if v[0] == 0]:
            del seqPool[ki]
        print("sequence selection has been carried out")
        return seqPool

    def samplingProcess(self, apt,
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
                seq = apt.pseudoAptamerGenerator(seqIdx)
                s.write(str(seq)+'\t'+str(int(seqPool[seqIdx][1]))+'\t'+str(N)+'\n')
        return

    # This function takes an empty selected pool, aptamer sequence structure and loop,
    # number of target binding sites, the alphabet set of the molecule, length,
    # total sequence number and stringency factor and returns full selected pool
    # Input: np.array(), int(), stats.obj(), int()
    # Output: np.array()
    def selectionProcess(self, seqPool, selectionDist, seqLength):
        selectedSeqs = 0
        # until all sites are occupied
        print("Drawing sample batch")
        while(selectedSeqs < self.selectionThreshold):
            # draw random sequences
            # warning: looping here causes large memory consumption
            for randIdx in selectionDist.rvs(size=Nrsamples):
                # carry out stochastic selection
                # draw random affinities
                if(int(seqPool[randIdx][1]) < utils.randint(0, seqLength-self.stringency)):
                    seqPool[randIdx][0] += 1
                    selectedSeqs += 1
                    if selectedSeqs % Nrsamples == 0:
                        print("{}% completed".format(100.0*selectedSeqs/self.selectionThreshold))
        return
