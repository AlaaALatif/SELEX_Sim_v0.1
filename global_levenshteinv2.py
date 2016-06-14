import time
import random
import linecache
from itertools import izip, imap
import operator
from collections import OrderedDict


def optimumseqs(num_optm_seqs, init_pool_file):
    with open(init_pool_file) as f: #open initial library 
        for totalseqs, l in enumerate(f): #calculate total num of seqs
            pass
        totalseqs+=1
    print("optimumseqs 1st loop passed")
    print(totalseqs)
    optimseqs_idx={}
    optimseqs={}
    i=1
    while(i <= num_optm_seqs):
       optimseq_idx=random.randint(0,totalseqs) #random seq index
       if optimseq_idx not in optimseqs_idx:
           optimseqs_idx.setdefault(optimseq_idx, [])
           optimseq = str(linecache.getline(init_pool_file, optimseq_idx))
           optimseqs.setdefault(optimseq, [])
       i+=1
    print("optimumseqs 2nd loop passed")
    print(optimseqs)
    return optimseqs, totalseqs

def hamming_func(str1, str2):
   assert len(str1) == len(str2)
   ne = operator.ne
   return sum(imap(ne, str1, str2))


def seqs_hamming(optimseqs, pool_file):
    seqs={}
    seqsfile = open(pool_file)
    seq=seqsfile.readline()
    while(seq):
        if seq not in seqs:
            seqs.setdefault(seq, []).append(1)
            for optimseq in optimseqs:
               hammdist = hamming_func(optimseq, seq)
               seqs.setdefault(seq, []).append(hammdist)
        else:
            seqs[seq][0]+=1
        seq=seqsfile.readline()
    seqsfile.close()
    print("distance calculations passed")
    return seqs

def selection(seqs, selection_rate, totalseqs):
   seqs_ordrd = OrderedDict(sorted(seqs.items(), key=lambda seqs: seqs[1][1], reverse=False)) #sort seqs by hammdist
   print("seqs have been sorted by distance value in descending order")
   #i=0
   slctd_seqs={}
   sampled_seqs=0
   print("parameters for selection have been initialized")
   # select top 20 percent in terms of hamm distance
   for seq in seqs_ordrd:
       if(sampled_seqs <= (selection_rate*totalseqs)):
       # THIS IS TAKING TOO LONG - NOT ANYMORE
           slctd_seqs[seq] = seqs_ordrd[seq] #add to selected pool
           sampled_seqs += slctd_seqs[seq][0] #add number of repeats
       else:
           break
       #i+=1 #increment index
   #seqs_ordrd.clear() #clear table of previous round
   print("sequence selection has been carried out")
   return slctd_seqs


def amplification(slctd_seqs, pcr_cycles):
   for seq in slctd_seqs:
      slctd_seqs[seq][0]*=(2**pcr_cycles) #assume no bias and 100% efficiency
   amplfd_seqs = slctd_seqs
   print("sequence amplification has been carried out")
   return amplfd_seqs


##TO DO
# amplify selected sequences randomly with bias

# repeat for required number of rounds

# select sequences stochastically

# add random mutations during PCR
        
          

