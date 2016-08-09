import time
import random
import linecache
from itertools import izip, imap
import operator
from collections import OrderedDict
from scipy import stats
import numpy as np




def seqNum(self, seqPool):
    for seqIdx in seqPool:
        totalSeqNum += seqPool[seqIdx][0]
        uniqSeqNum += 1

    return totalSeqNum, uniqSeqNum

def convert_to_distribution(x, y, distName):
    xDist = stats.rv_discrete(name=distName, values=(x, y))
    return xDist

# Add method for computing the binomial coefficient

# Add method for computing the L1 norm 

# Add method to convert probability vectors to discrete distributions
