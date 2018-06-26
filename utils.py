import random
from math import factorial
from scipy import stats
import numpy as np


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


# This converts an array of probabilities into a
# discrete probability distribution
def convert_to_distribution(x, y, distName):
    xDist = stats.rv_discrete(name=distName, values=(x, y))
    return xDist


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


def rvd(X, X_sum, distName):
    seqIdxs = np.zeros(X.shape[0])
    probs = np.zeros(X.shape[0])
    for i, seq in enumerate(X):
        seqIdxs[i] = i
        probs[i] = seq[1]/X_sum
    return stats.rv_discrete(name=distName, values=(seqIdxs, probs))


def batch_size(size, Nbatch):
    i = 0
    while size-i > Nbatch:
        yield(Nbatch)
        i += Nbatch
    else:
        yield(size-i)


def rvs_iter(dist, size, Nbatch=1000):
    for S in batch_size(size, Nbatch):
        for v in dist.rvs(size=Nbatch):
            yield v


# long random numbers, via random
def randint(a, b, size=1):
    return np.array([random.randint(a, b) for i in range(size)])
