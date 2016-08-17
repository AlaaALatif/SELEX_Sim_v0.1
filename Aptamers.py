import time
import math
import numpy as np
import random
import linecache
from itertools import izip, imap, islice, product
import operator
from collections import OrderedDict
from scipy import stats 

class Aptamers:
    #Add __init__ constructor here

#Generate any sequence given it's index, length and the alphabet set
#Sequences are indexed in order of alphabet set provided
#Ex1:   if alphabetSet = 'ATCG', length = 4, index = 0 --> 'AAAA'
#Ex2:   index = 15 --> 'GGGG'
    def pseudoAptamerGenerator(self, seqIdx, alphabetSet, seqLen):
        seq = str() #initialize seq 
        seqArray = np.zeros(seqLen) 
        alphabetSetSize = len(alphabetSet)
        assert seqIdx >= 0
        assert seqIdx <= alphabetSetSize**(seqLen) - 1
        while(seqIdx>0):
            charIdx = int(math.floor(math.log(seqIdx, alphabetSetSize)))
            if(charIdx > 0):
                seqArray[charIdx] += 1 #mutate lexicographically
                seqIdx -= alphabetSetSize**charIdx #next seqIdx
            else:
                seqArray[charIdx] = seqIdx
                break
        for charCode in seqArray:
            for char in alphabetSet:
                if(charCode == alphabetSet.index(char)):
                    seq += alphabetSet[int(charCode)]
        seq = seq[::-1] #reverse string
        assert len(seq) == seqLen
        return seq

#TEST AREA - DELETE
#apt = Aptamers()
#testSeq = apt.pseudoAptamerGenerator(1099511627775, 'ACGT', 20)

# method to get seqArray given seq index
    def get_seqArray(self, seqIdx, alphabetSet, seqLen):
        seq = str() #initialize seq 
        seqArray = np.zeros(seqLen) 
        alphabetSetSize = len(alphabetSet)
        assert seqIdx > 0
        assert seqIdx <= alphabetSetSize**(seqLen) - 1
        while(seqIdx>0):
            charIdx = int(math.floor(math.log(seqIdx, alphabetSetSize)))
            if(charIdx > 0):
                seqArray[charIdx] += 1 #mutate lexicographically
                seqIdx -= alphabetSetSize**charIdx #next seqIdx
            else:
                seqArray[charIdx] = seqIdx
                break
        seqArray = seqArray[::-1] #reverse string
        return seqArray

#Generate any sequence index given the sequence, length and the alphabet set
#Sequences are indexed in order of alphabet set provided
#Ex1:   if alphabetSet = 'ATCG', length = 4, seq = 'AAAA' --> 0
#Ex2:   seq = 'GGGG' --> 15
    def pseudoAptamerIndexGenerator(self, seq, alphabetSet, seqLen):
        assert len(seq) == seqLen
        alphabetSize = len(alphabetSet)
        seq = seq[::-1] #reverse seq
        seqIdx = 0
        for ntPos, nt in enumerate(seq):
            seqIdx += alphabetSet.index(nt)*(alphabetSize)**(ntPos) 
        return seqIdx
#TEST AREA - DELETE    
#apt = Aptamers()
#apt.pseudoAptamerIndexGenerator('TTTTTTTTTTTTTTTTTTTT', 'ACGT', 20)

    def pseudoAptamerIterator(self, alphabetSet, seqLen):
        initLibrary = product(alphabetSet, repeat=seqLen)
        return initLibrary

###NEED TO MODIFY TO ALLOW MULTIPLE OPTIMUM APTAMERS

    def optimumAptamerGenerator(self, aptamerNum, alphabetSet, seqLen):
        seq = str() #initialize seq 
        seqArray = np.zeros(seqLen) 
        alphabetSetSize = len(alphabetSet)
        initialSeqNum = alphabetSetSize**(seqLen)
        seqIdx = random.randint(0, initialSeqNum - 1)
        while(seqIdx>0):
            charIdx = int(math.floor(math.log(seqIdx, alphabetSetSize)))
            if(charIdx > 0):
                seqArray[charIdx] += 1 #mutate lexicographically
                seqIdx -= alphabetSetSize**charIdx #next seqIdx
            else:
                seqArray[charIdx] = seqIdx
                break
        for charCode in seqArray:
            for char in alphabetSet:
                if(charCode == alphabetSet.index(char)):
                    seq += alphabetSet[int(charCode)]
        seq = seq[::-1] #reverse string
        return seq, initialSeqNum

