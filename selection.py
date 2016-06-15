import time
import random
import linecache
from itertools import izip, imap
import operator
from collections import OrderedDict
from scipy import stats

class selection

def definite_selection(seqs, selection_rate, totalseqs):
   seqs_ordrd = OrderedDict(sorted(seqs.items(), key=lambda seqs: seqs[1][1], reverse=False)) #sort seqs by hammdist
   print("seqs have been sorted by distance value in descending order")
   #i=0
   slctd_seqs={}
   sampled_seqs=0
   print("parameters for selection have been initialized")
   # select top 20 percent in terms of hamm distance
   for seq in seqs_ordrd:
       if(sampled_seqs <= (selection_rate*totalseqs)):
       # THIS IS TAKING TOO LONG # UPDATE: NOT ANYMORE
           slctd_seqs[seq] = seqs_ordrd[seq] #add to selected pool
           sampled_seqs += slctd_seqs[seq][0] #add number of repeats
       else:
           break
       #i+=1 #increment index
   #seqs_ordrd.clear() #clear table of previous round
   print("sequence selection has been carried out")
   return slctd_seqs

def stochastic_selection(seqs, selection_rate, totalseqs):
   seqs_ordrd = OrderedDict(sorted(seqs.items(), key=lambda seqs: seqs[1][1], reverse=False)) #sort seqs by hammdist
   print("seqs have been sorted by distance value in descending order")
   #i=0
   slctd_seqs={}
   sampled_seqs=0
   print("parameters for selection have been initialized")
   # select top 20 percent in terms of hamm distance
   for seq in seqs_ordrd:
       if(sampled_seqs <= (selection_rate*totalseqs)):
           randHamm = random.randint(0, seqLength)
           if(seqs_ordrd[seq][1] <= randHamm):
       # THIS IS TAKING TOO LONG # UPDATE: NOT ANYMORE
               slctd_seqs[seq] = seqs_ordrd[seq] #add to selected pool
               sampled_seqs += slctd_seqs[seq][0] #add number of repeats
           else:
               continue
       else:
           break
       #i+=1 #increment index
   #seqs_ordrd.clear() #clear table of previous round
   print("sequence selection has been carried out")
   return slctd_seqs

