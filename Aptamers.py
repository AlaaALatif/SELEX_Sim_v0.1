import math
import numpy as np
import random
from itertools import islice, product


class Aptamers:
    # Add __init__ constructor here
    def __init__(self, alphabetSet, seqLength):
        self.alphabetSet = alphabetSet
        self.seqLength = seqLength
        self.La = len(self.alphabetSet)
        self.td = str.maketrans(dict(zip(self.alphabetSet, "0123")))

    # Generate any sequence given it's index, length and the alphabet set
    # Sequences are indexed in order of alphabet set provided
    # Ex1:   if alphabetSet = 'ATCG', length = 4, index = 0 --> 'AAAA'
    # Ex2:   index = 15 --> 'GGGG'
    def pseudoAptamerGenerator_(self, seqIdx):
        seq = str()  # initialize seq
        seqArray = np.zeros(self.seqLength)
        assert seqIdx >= 0
        assert seqIdx <= self.La**(self.seqLength) - 1
        while seqIdx > 0:
            charIdx = int(math.floor(math.log(seqIdx, self.La)))
            if charIdx > 0:
                seqArray[charIdx] += 1  # mutate lexicographically
                seqIdx -= self.La**charIdx  # next seqIdx
            else:
                seqArray[charIdx] = seqIdx
                break
        for charCode in seqArray:
            for char in self.alphabetSet:
                if(charCode == self.alphabetSet.index(char)):
                    seq += self.alphabetSet[int(charCode)]
        seq = seq[::-1]  # reverse string
        assert len(seq) == self.seqLength
        return seq

    def pseudoAptamerGenerator(self, sn):
        sn = int(sn)
        sl = ""
        for i in range(self.seqLength):
            sl += self.alphabetSet[sn % self.La]
            sn = sn // self.La
        return sl[::-1]

    # method to get seqArray given seq index
    def get_seqArray(self, seqIdx):
        seqArray = np.zeros(self.seqLength)
        assert seqIdx > 0
        assert seqIdx <= self.La**(self.seqLength) - 1
        while seqIdx > 0:
            charIdx = int(math.floor(math.log(seqIdx, self.La)))
            if charIdx > 0:
                seqArray[charIdx] += 1  # mutate lexicographically
                seqIdx -= self.La**charIdx  # next seqIdx
            else:
                seqArray[charIdx] = seqIdx
                break
        seqArray = seqArray[::-1]  # reverse string
        return seqArray

    # Generate any sequence index given the sequence, length and the alphabet set
    # Sequences are indexed in order of alphabet set provided
    # Ex1:   if alphabetSet = 'ATCG', length = 4, seq = 'AAAA' --> 0
    # Ex2:   seq = 'GGGG' --> 15
    def pseudoAptamerIndexGenerator_(self, seq):
        assert len(seq) == self.seqLength
        seq = seq[::-1]  # reverse seq
        seqIdx = 0
        for ntPos, nt in enumerate(seq):
            seqIdx += self.alphabetSet.index(nt)*(self.La)**(ntPos)
        return seqIdx

    def pseudoAptamerIndexGenerator(self, seq):
        return int(seq.translate(self.td), self.La)

    def pseudoAptamerIterator(self):
        initLibrary = product(self.alphabetSet, repeat=self.seqLength)
        return initLibrary


# choose a random subset of sequences to be aptamers from the initial pool
# input params are the number of aptamers to choose and the initial pool
# sequence file
# returns total number of seqs in initial pool and set of aptamers
    def randomAptamerChooser(self, aptamerNum, initLib):
        initialSeqNum = 4**(self.seqLength)
        optimumAptamers = np.chararray(aptamerNum, itemsize=self.seqLength)
        for aptNum in range(aptamerNum):
            aptIdx = random.randint(0, initialSeqNum)  # random seq index
            aptSeqList = list(islice(initLib, aptIdx, aptIdx+1))[0]
            aptSeq = 0
            for residue in aptSeqList:
                aptSeq += residue
            optimumAptamers[aptNum] = aptSeq
        return optimumAptamers, initialSeqNum

    # Generate all possible sequences
    def aptamerGenerator(self, start, finish, outFile):
        initialLibrary = product(self.alphabetSet, repeat=self.seqLength)
        # poolFraction = list(islice(initialLibrary, start, finish))
        with open(outFile, 'w') as o:
            for i in range(start, finish, 1000000):
                seqList = list(islice(initialLibrary, (start+i), (start+i+1000000)))
                for seq in seqList:
                    sequence = str()
                    for charIdx in range(self.seqLength):
                        sequence += seq[charIdx]
                    o.write(sequence+'\n')
        o.close()
        print("first generation completed")
        return initialLibrary

###NEED TO MODIFY TO ALLOW MULTIPLE OPTIMUM APTAMERS

    def optimumAptamerGenerator(self, aptamerNum):
        seq = str()  # initialize seq
        seqArray = np.zeros(self.seqLength)
        initialSeqNum = self.La**(self.seqLength)
        seqIdx = random.randint(0, initialSeqNum - 1)
        while seqIdx > 0:
            charIdx = int(math.floor(math.log(seqIdx, self.La)))
            if(charIdx > 0):
                seqArray[charIdx] += 1  # mutate lexicographically
                seqIdx -= self.La**charIdx  # next seqIdx
            else:
                seqArray[charIdx] = seqIdx
                break
        for charCode in seqArray:
            for char in self.alphabetSet:
                if(charCode == self.alphabetSet.index(char)):
                    seq += self.alphabetSet[int(charCode)]
        seq = seq[::-1]  # reverse string
        return seq, initialSeqNum
