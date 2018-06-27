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
Nrsamples = 10**3


class Selection:
    def __init__(self):
        self.distances = ("hamming", "basepair", "loop", "random")
        self.select_init = dict(zip(self.distances, (self.stochasticHammingSelection_initial,
                                                     self.stochasticBasePairSelection_initial,
                                                     self.stochasticLoopSelection_initial,
                                                     self.randomSelection_initial)))
        self.select = dict(zip(self.distances, (self.stochasticHammingSelection,
                                                self.stochasticBasePairSelection,
                                                self.stochasticLoopSelection,
                                                self.randomSelection)))

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

    # This function takes an empty selected pool, aptamer sequence structure and loop,
    # number of target binding sites, the alphabet set of the molecule, length,
    # total sequence number and stringency factor and returns full selected pool
    # their Lavenshtein distance
    # Input: dict(), str(), str(), str(), int(), str(), int(), int(), int()
    # Output: dict()
    def stochasticLoopSelection_initial(self, alphabetSet, seqLength,
                                        aptPool, selectionThreshold,
                                        totalSeqNum, samplingSize,
                                        outputFileNames, rnd, stringency):
        # sampling
        print("sampling from initial library...")
        randomSamples = utils.randint(0, int(totalSeqNum-1), size=samplingSize)
        sampleFileName = outputFileNames+"_samples_R{:03d}".format(rnd+1)
        with open(sampleFileName, 'w') as s:
            for seqIdx in randomSamples:
                seq = Apt.pseudoAptamerGenerator(seqIdx, alphabetSet, seqLength)
                s.write(seq+'\n')
        print("Sampling completed")
        # initialize seqInfo matrix
        slctdSeqs = {}
        aptStruct = RNA.fold(aptPool)[0]
        print("Optimum aptamer structure: {}".format(aptStruct))
        aptLoop = utils.apt_loopFinder(aptPool, aptStruct, seqLength)
        print("Selection has started")
        # stochastic selection until threshold is met
        slctdSeqs = self.selectionProcess_loop_initial(slctdSeqs, aptPool,
                                                       aptStruct, aptLoop,
                                                       selectionThreshold,
                                                       alphabetSet, seqLength,
                                                       totalSeqNum, stringency)
        print("sequence selection has been carried out")
        return slctdSeqs

    def stochasticLoopSelection(self, alphabetSet, seqLength,
                                seqPool, selectionThreshold,
                                uniqSeqNum, totalSeqNum, samplingSize,
                                outputFileNames, rnd, stringency):
        # initialize selected sequence pool
        slctdSeqs = {}
        print("seq length = "+str(seqLength))
        print("seq selection threshold = "+str(selectionThreshold))
        print("parameters for selection have been initialized")
        x = np.zeros((uniqSeqNum, 4), dtype="object")
        for i, seqIdx in enumerate(seqPool):
            x[i][0] = seqIdx
            x[i][1] = seqPool[seqIdx][0]
            x[i][2] = seqPool[seqIdx][1]
            x[i][3] = seqPool[seqIdx][2]
        print("Selection sample distribution being computed...")
        # compute sampling distribution for selection
        selectionDist = utils.rvd(x, totalSeqNum, "selectionDist")
        print("Selection sample distribution computed")
        print("Sampling has started...")
        sampleFileName = outputFileNames+"_samples_R{:03d}".format(rnd+1)
        with open(sampleFileName, 'w') as s:
            for seqIdx in utils.rvs_iter(selectionDist, samplingSize):
                seq = Apt.pseudoAptamerGenerator(x[seqIdx][0], alphabetSet, seqLength)
                s.write(str(seq)+'\t'+str(int(x[seqIdx][2]))+'\n')
        print("Sampling has completed")
        for i, seqIdx in enumerate(x):
            x[i][1] = 0
        x = self.selectionProcess(x, selectionThreshold,
                                  selectionDist, seqLength,
                                  stringency)
        x = x[x[:, 1] != 0]
        for seqInfo in x:
            # change it so that seqInfo are added as on np array, without setdefault
            slctdSeqs[int(seqInfo[0])] = seqInfo[1:]
        print("sequence selection has been carried out")
        return slctdSeqs

    def selectionProcess_2D_initial(self, slctdSeqs, aptStruct,
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
                randSeqDist = D.bp_func(aptStruct, randSeq)
                if(selectedSeqs == selectionThreshold):
                    return slctdSeqs
                elif(randSeqDist < randHamms[i]):
                    if(randIdx in slctdSeqs):
                        slctdSeqs[randIdx][0] += 1
                    else:
                        randSeqBias = D.bias_func(randSeq, seqLength)
                        slctdSeqs[randIdx] = np.array([1, randSeqDist, randSeqBias])
                    selectedSeqs += 1

    def stochasticBasePairSelection_initial(self, alphabetSet, seqLength,
                                            aptPool, selectionThreshold,
                                            totalSeqNum, samplingSize,
                                            outputFileNames, rnd, stringency):
        # sampling
        print("sampling from initial library...")
        randomSamples = utils.randint(0, int(totalSeqNum-1), size=samplingSize)
        sampleFileName = outputFileNames+"_samples_R{:03d}".format(rnd+1)
        with open(sampleFileName, 'w') as s:
            for seqIdx in randomSamples:
                seq = Apt.pseudoAptamerGenerator(seqIdx, alphabetSet, seqLength)
                s.write(seq+'\n')
        print("Sampling completed")
        # initialize seqInfo matrix
        slctdSeqs = {}
        aptStruct = RNA.fold(aptPool)[0]
        print("Optimum aptamer structure: {}".format(aptStruct))
        print("Selection has started")
        # stochastic selection until threshold is met
        slctdSeqs = self.selectionProcess_2D_initial(slctdSeqs, aptStruct,
                                                     selectionThreshold,
                                                     alphabetSet, seqLength,
                                                     totalSeqNum, stringency)
        print("sequence selection has been carried out")
        return slctdSeqs

    def stochasticBasePairSelection(self, alphabetSet, seqLength,
                                    seqPool, selectionThreshold,
                                    uniqSeqNum, totalSeqNum, samplingSize,
                                    outputFileNames, rnd, stringency):
        # initialize selected sequence pool
        slctdSeqs = {}
        print("seq length = "+str(seqLength))
        print("seq selection threshold = "+str(selectionThreshold))
        print("parameters for selection have been initialized")
        x = np.zeros((uniqSeqNum, 4), dtype="object")
        for i, seqIdx in enumerate(seqPool):
            x[i][0] = seqIdx
            x[i][1] = seqPool[seqIdx][0]
            x[i][2] = seqPool[seqIdx][1]
            x[i][3] = seqPool[seqIdx][2]
        print("Selection sample distribution being computed...")
        # compute sampling distribution for selection
        selectionDist = utils.rvd(x, totalSeqNum, "selectionDist")
        print("Selection sample distribution computed")
        print("Sampling has started...")
        sampleFileName = outputFileNames+"_samples_R{:03d}".format(rnd+1)
        with open(sampleFileName, 'w') as s:
            for seqIdx in utils.rvs_iter(selectionDist, samplingSize):
                seq = Apt.pseudoAptamerGenerator(x[seqIdx][0], alphabetSet, seqLength)
                s.write(str(seq)+'\t'+str(int(x[seqIdx][2]))+'\n')
        print("Sampling has completed")
        for i, seqIdx in enumerate(x):
            x[i][1] = 0
        x = self.selectionProcess(x, selectionThreshold,
                                  selectionDist, seqLength,
                                  stringency)
        x = x[x[:, 1] != 0]
        for seqInfo in x:
            # change it so that seqInfo are added as on np array, without setdefault
            slctdSeqs[int(seqInfo[0])] = seqInfo[1:]
        print("sequence selection has been carried out")
        return slctdSeqs

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
                randSeqDist = distf(randSeq, aptPool)
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
        return slctdSeqs

    def stochasticHammingSelection_initial(self, alphabetSet, seqLength,
                                           aptPool, selectionThreshold,
                                           totalSeqNum, samplingSize,
                                           outputFileNames, rnd, stringency):
        # initialize seqInfo matrix
        slctdSeqs = {}
        print("sampling from initial library...")
        print("Selection has started...", flush=True)
        # sampling
        slctdSeqs = self.selectionProcess_1D_initial(slctdSeqs,
                                                     aptPool, selectionThreshold,
                                                     alphabetSet, seqLength,
                                                     totalSeqNum, stringency)
        print("sequence selection has been carried out")
        sampleFileName = outputFileNames+"_samples_R{:03d}".format(rnd+1)
        with open(sampleFileName, 'w') as s:
            for sidx, x in slctdSeqs.items():
                seq = Apt.pseudoAptamerGenerator(sidx, alphabetSet, seqLength)
                s.write(str(seq)+'\t'+str(int(x[1]))+'\t'+str(int(x[0]))+'\n')
        print("Sampling completed")
        return slctdSeqs

    def stochasticHammingSelection(self, alphabetSet, seqLength,
                                   seqPool, selectionThreshold,
                                   uniqSeqNum, totalSeqNum, samplingSize,
                                   outputFileNames, rnd, stringency):
        # initialize selected sequence pool
        slctdSeqs = {}
        print("seq length = "+str(seqLength))
        print("seq selection threshold = "+str(selectionThreshold))
        print("parameters for selection have been initialized")
        x = np.zeros((uniqSeqNum, 4), dtype="object")
        # transfer selected pool to matrix x
        for i, seqIdx in enumerate(seqPool):
            x[i][0] = seqIdx
            # seq count
            x[i][1] = seqPool[seqIdx][0]
            # seq distance
            x[i][2] = seqPool[seqIdx][1]
            # seq bias
            x[i][3] = seqPool[seqIdx][2]
        print("Selection sample distribution being computed...")
        # distribution computed using count of each unique seq
        selectionDist = utils.rvd(x, totalSeqNum, "selectionDist")
        print("Selection sample distribution computed")
        print("Sampling has started...")
        samps = dict()
        # draw random samples from distribution
        for seqIdx in utils.rvs_iter(selectionDist, samplingSize):
            if seqIdx in samps:
                samps[seqIdx] += 1
            else:
                samps[seqIdx] = 1
        sampleFileName = outputFileNames+"_samples_R{:03d}".format(rnd+1)
        # write to samples file
        with open(sampleFileName, 'w') as s:
            for seqIdx, N in samps.items():
                seq = Apt.pseudoAptamerGenerator(x[seqIdx][0], alphabetSet, seqLength)
                s.write(str(seq)+'\t'+str(int(x[seqIdx][2]))+'\t'+str(N)+'\n')
        print("Sampling has completed")
        # reset all seq counts prior to selection
        for i, seqIdx in enumerate(x):
            x[i][1] = 0
        # draw a bunch of random seqs
        x = self.selectionProcess(x, selectionThreshold, selectionDist,
                                  seqLength, stringency)
        # remove all seqs that haven't been selected
        x = x[x[:, 1] != 0]
        # transfer info back to selected pool
        for seqInfo in x:
            # change it so that seqInfo are added as on np array, without setdefault
            slctdSeqs[int(seqInfo[0])] = seqInfo[1:]
        print("sequence selection has been carried out")
        return slctdSeqs

    # This function takes an empty selected pool, aptamer sequence structure and loop,
    # number of target binding sites, the alphabet set of the molecule, length,
    # total sequence number and stringency factor and returns full selected pool
    # Input: np.array(), int(), stats.obj(), int()
    # Output: np.array()
    def selectionProcess(self, x, selectionThreshold,
                         selectionDist, seqLength,
                         stringency):
        selectedSeqs = 0
        # until all sites are occupied
        print("Drawing sample batch")
        while(selectedSeqs < selectionThreshold):
            print("{}% completed".format(100.0*selectedSeqs/selectionThreshold))
            # draw random sequences
            randIdxs = selectionDist.rvs(size=Nrsamples)
            # draw random affinities
            randHamms = utils.randint(0, seqLength-stringency, size=Nrsamples)
            # carry out stochastic selection
            for i, randIdx in enumerate(randIdxs):
                if(selectedSeqs == selectionThreshold):
                    return x
                elif(int(x[randIdx][2]) < randHamms[i]):
                    x[randIdx][1] += 1
                    selectedSeqs += 1
        return x

    def randomSelection_initial(self, alphabetSet, seqLength,
                                aptPool, selectionThreshold,
                                totalSeqNum, samplingSize,
                                outputFileNames, rnd, stringency):
        # sampling
        print("sampling from initial library...")
        randomSamples = utils.randint(0, int(totalSeqNum-1), size=samplingSize)
        sampleFileName = outputFileNames+"_samples_R{:03d}".format(rnd+1)
        with open(sampleFileName, 'w') as s:
            for seqIdx in randomSamples:
                seq = Apt.pseudoAptamerGenerator(seqIdx, alphabetSet, seqLength)
                s.write(seq+'\n')
        print("Sampling completed")
        # initialize seqInfo matrix
        slctdSeqs = {}
        print("Selection has started...", flush=True)
        slctdSeqs = self.selectionProcess_1D_initial(slctdSeqs,
                                                     aptPool, selectionThreshold,
                                                     alphabetSet, seqLength,
                                                     totalSeqNum, stringency,
                                                     distf=D.nodist_func)
        print("sequence selection has been carried out")
        return slctdSeqs

    def randomSelection(self, alphabetSet, seqLength,
                        seqPool, selectionThreshold,
                        uniqSeqNum, totalSeqNum, samplingSize,
                        outputFileNames, rnd, stringency):
        # initialize selected sequence pool
        slctdSeqs = {}
        print("seq length = "+str(seqLength))
        print("seq selection threshold = "+str(selectionThreshold))
        print("parameters for selection have been initialized")
        x = np.zeros((uniqSeqNum, 4), dtype="object")
        # transfer selected pool to matrix x
        for i, seqIdx in enumerate(seqPool):
            x[i][0] = seqIdx
            # seq count
            x[i][1] = seqPool[seqIdx][0]
            # seq distance
            x[i][2] = seqPool[seqIdx][1]
            # seq bias
            x[i][3] = seqPool[seqIdx][2]
        print("Selection sample distribution being computed...")
        # distribution computed using count of each unique seq
        selectionDist = utils.rvd(x, totalSeqNum, "selectionDist")
        print("Selection sample distribution computed")
        print("Sampling has started...")
        # write to samples file
        sampleFileName = outputFileNames+"_samples_R{:03d}".format(rnd+1)
        with open(sampleFileName, 'w') as s:
            # draw random samples from distribution
            for seqIdx in utils.rvs_iter(selectionDist, samplingSize):
                seq = Apt.pseudoAptamerGenerator(x[seqIdx][0], alphabetSet, seqLength)
                s.write(str(seq)+'\t'+str(int(x[seqIdx][2]))+'\n')
        print("Sampling has completed")
        # reset all seq counts prior to selection
        for i, seqIdx in enumerate(x):
            x[i][1] = 0
        # draw a bunch of random seqs
        x = self.selectionProcess(x, selectionThreshold, selectionDist,
                                  seqLength, stringency)
        # remove all seqs that haven't been selected
        x = x[x[:, 1] != 0]
        # transfer info back to selected pool
        for seqInfo in x:
            # change it so that seqInfo are added as on np array, without setdefault
            slctdSeqs[int(seqInfo[0])] = seqInfo[1:]
        print("sequence selection has been carried out")
        return slctdSeqs
