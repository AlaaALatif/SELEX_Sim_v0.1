import time
import random
import linecache
from itertools import izip, imap
import operator
from collections import OrderedDict
from scipy import stats
import numpy as np



class utils:

    def seqNum(self, seqPool):
        for seqIdx in seqPool:
            totalSeqNum += seqPool[seqIdx][0]
            uniqSeqNum += 1

        return totalSeqNum, uniqSeqNum
