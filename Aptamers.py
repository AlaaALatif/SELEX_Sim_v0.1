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


# choose a random subset of sequences to be aptamers from the initial pool
# input params are the number of aptamers to choose and the initial pool 
# sequence file
# returns total number of seqs in initial pool and set of aptamers 
    def randomAptamerChooser(self, aptamerNum, seqLen, initLib):                                 
        initialSeqNum = 4**(seqLen)
        optimumAptamers = np.chararray(aptamerNum, itemsize=seqLen)
        for aptNum in xrange(aptamerNum):
            aptIdx = random.randint(0, initialSeqNum) #random seq index
            aptSeqList = list(islice(initLib, aptIdx, aptIdx+1))[0]
            for residue in aptSeqList:
                aptSeq += residue
            optimumAptamers[aptNum] = aptSeq
        return optimumAptamers, initialSeqNum
# Generate all possible sequences 
    def aptamerGenerator(self, alphabetSet, seqLength, start, finish, outFile):
        initialLibrary = product(alphabetSet, repeat=seqLength)
        #poolFraction = list(islice(initialLibrary, start, finish))
        with open(outFile, 'w') as o:
            for i in xrange(start, finish, 1000000):
                seqList = list(islice(initialLibrary, (start+i), (start+i+1000000)))
                for seq in seqList:
                    sequence = str()
                    for charIdx in range(seqLength):
                        sequence += seq[charIdx]
                    o.write(sequence+'\n')
        o.close()
        print("first generation completed")
        return initialLibrary

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



# choose a random subset of sequences to be aptamers from the initial pool
# input params are the number of aptamers to choose and the initial pool 
# sequence file
# returns total number of seqs in initial pool and set of aptamers 
    def optimumSeqs(self, num_optm_seqs, init_pool_file):                                 
        with open(init_pool_file) as f: #open initial library
            for totalseqs, l in enumerate(f): #calculate total num of seqs
                pass
            totalseqs+=1
        f.close()
        print(totalseqs)
        optimseqs_idx={}
        optimseqs={}
        i=1
        while(i <= num_optm_seqs):
            optimseq_idx=random.randint(0,totalseqs) #random seq index
            if optimseq_idx not in optimseqs_idx:
                optimseqs_idx.setdefault(optimseq_idx, [])                          
                optimseq = str(open(init_pool_file).readlines()[optimseq_idx])
                optimseqs.setdefault(optimseq, [])
            i+=1
        print(optimseqs)
        return optimseqs, totalseqs
