import random
from math import factorial
import numpy as np
import numpy.random as nr


def seqNumberCounter(seqPool):
    totalSeqNum = int(0)
    uniqSeqNum = int(0)
    for seqIdx in seqPool:
        totalSeqNum += seqPool[seqIdx][0]
        uniqSeqNum += 1
    return int(totalSeqNum), int(uniqSeqNum)


# This computes the binomial coefficient (not used)
def binomCoeff(n, k):
    binom = factorial(n)/(factorial(k)*factorial(n-k))
    return binom


# This finds the loop region in a given sequence
def apt_loopFinder(apt_seq, apt_struct, seqLength):
    base = None
    baseIdx = 0
    while(base != ')' and baseIdx < seqLength):
        base = apt_struct[baseIdx]
        baseIdx += 1
    if(baseIdx == seqLength):
        while(base != '('and baseIdx > 1):
            baseIdx -= 1
            base = apt_struct[baseIdx-1]
        if(baseIdx == 1):
            apt_loop = apt_seq
            return apt_loop
        else:
            apt_loop = apt_seq[baseIdx:]
            return apt_loop
    else:
        loop_end = baseIdx-1
        while(base != '(' and baseIdx > 1):
            baseIdx -= 1
            base = apt_struct[baseIdx-1]
        if(baseIdx == 1):
            apt_loop = apt_seq[:loop_end]
            return apt_loop
        else:
            apt_loop = apt_seq[baseIdx:loop_end]
            return apt_loop

# Add method for computing the binomial coefficient

# Add method for computing the L1 norm

# Add method to convert probability vectors to discrete distributions


class rv_int():
    def __init__(self, seqPool, distName):
        self.si = list(seqPool)
        self.probas = np.array([seqPool[i][0] for i in self.si], dtype=np.float64)
        npb = self.probas[self.probas < -0.1]
        if len(npb) > 0:
            print("ERROR", npb)
        self.probas[self.probas < 0] = 0
        self.probas /= self.probas.sum()

    def rvs(self, size=1):
        return nr.choice(self.si, p=self.probas, size=size)


def batch_size(size, Nbatch):
    i = 0
    while size-i > Nbatch:
        yield(Nbatch)
        i += Nbatch
    else:
        yield(size-i)


# long random numbers, via random
def randint(a, b, size=1):
    return np.array([random.randint(a, b) for i in range(size)])
