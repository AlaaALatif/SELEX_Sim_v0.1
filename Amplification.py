import time
import random
import linecache
from itertools import izip, imap
import operator
from collections import OrderedDict
from scipy import stats
import numpy as np

class Amplification:

#pcr modelled as ideal process that doubles each sequence per cycle
#assume no bias and 100% efficiency
   def ampIdeal(self, slctd_seqs, pcr_cycles):
      for seq in slctd_seqs:
         slctd_seqs[seq][0]*=(2**pcr_cycles) 
      amplfd_seqs = slctd_seqs
      print("sequence amplification has been carried out")
      return amplfd_seqs

##pcr modelled as bernoulli test per cycle
   def ampEfficiencyBrute(self, slctd_seqs, pcr_cycles, pcr_yld):
      n=0
      while(n<=pcr_cycles):
         for seq in slctd_seqs:
            slctd_seqs[seq][0] += np.random.binomial(slctd_seqs[seq][0], pcr_yld)
         n+=1
      amplfd_seqs = slctd_seqs
      print("sequence amplification has been carried out")
      return amplfd_seqs



##pcr modelled as Bernoulli process after n iterations
# assume no bias or mutations
# compute expectation after n cycles for each selected sequence
# estimate amplified number as expectation of pgf 
   def ampEfficiencyDefinite(self, slctd_seqs, pcr_cycles, pcr_yld):
      for seq in slctd_seqs:
         slctd_seqs[seq][0] = slctd_seqs[seq][0]*(1+pcr_yld)**(pcr_cycles)
      amplfd_seqs = slctd_seqs
      print("sequence amplification has been carried out")
      return amplfd_seqs



##pcr modelled as Bernoulli process after n iterations
# assume no bias or mutations
# compute expectation and variance after n cycles for each selected sequence
# estimate amplified number as draw from Gaussian distribution 
# Law of Large Numbers
# Binomial distribution approximated as Gaussian due to large copy numbers
   def ampEfficiencyStochastic(self, slctd_seqs, pcr_cycles, pcr_yld):
      for seq in slctd_seqs:
         mean_seq_num = (slctd_seqs[seq][0]*(1+pcr_yld))**(pcr_cycles) #expectation
         var_seq_num = (slctd_seqs[seq][0]*(1+pcr_yld)**(pcr_cycles - 1))*pcr_yld*(1-pcr_yld)*((slctd_seqs[seq][0]*(1+pcr_yld))**(pcr_cycles)-1) #variance
         slctd_seqs[seq][0] = int(stats.norm.rvs(mean_seq_num, var_seq_num)) #draw number from Gaussian dist  
      amplfd_seqs = slctd_seqs
      print("sequence amplification has been carried out")
      return amplfd_seqs



