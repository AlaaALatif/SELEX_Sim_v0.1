import time
import random
import linecache
from itertools import izip, imap, islice
import operator
from collections import OrderedDict
from scipy import stats

class Distance:

    def hamming_func(self, str1, str2):
        assert len(str1) == len(str2)
        ne = operator.ne
        return sum(imap(ne, str1, str2))
    
    def bias_func(self, seq, seqLen):
        pyrNum = 0
        for nt in seq[:-1]:
            if(nt == 'C') or (nt == 'T'):
                pyrNum += 1 #increment no. of pyrimidines
        biasScore = 0.1*(2*pyrNum - seqLen)/seqLen #compute bias
        return biasScore

