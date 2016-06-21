import time
import random
import linecache
from itertools import izip, imap
import operator
from collections import OrderedDict
from scipy import stats
import numpy as np
import math

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
         var_seq_num = ((slctd_seqs[seq][0]*(1+pcr_yld)**(pcr_cycles - 1))*slctd_seqs[seq][0]*pcr_yld*(1-pcr_yld)*((slctd_seqs[seq][0]*(1+pcr_yld))**(pcr_cycles)-1))/(slctd_seqs[seq][0] - 1) #variance
         slctd_seqs[seq][0] = int(stats.norm.rvs(mean_seq_num, var_seq_num)) #draw number from Gaussian dist  
      amplfd_seqs = slctd_seqs
      print("sequence amplification has been carried out")
      return amplfd_seqs



# Assume all mutations happened at end of PCR (i.e. mutants are not amplified till next round)
## ADD seqLength as arg to function
## EDITED variance calculation for amplification
   def ampEffMutStochastic(self, slctd_seqs, seqLength, pcr_cycles, pcr_yld, errorRate):
      mutatedPool = {} #initialize pool of seqs to be mutated
      gamma = np.arange(seqLength) #initialize vector of gamma values 
      mut_m = np.arange(seqLength) #initialize vector of mut no. (0-seqLen)
      prob_m = np.arange(seqLength) #initialize corresponding vector of probs
# Mutation Distribution
      for i in mut_m: # for each possible no. of muts, compute corresponding probability
        for j in range(pcr_cycles): # for each pcr cycle
          gamma[i] += np.exp(-j*errorRate*seqLength)*math.factorial(pcr_cycles)/(math.factorial(j)*math.factorial(pcr_cycles - j))*(j)**(i)*pcr_yld**(j) #calculate sum (SEE PAPER)
          prob_m[i] = (errorRate*seqLength)**(i)*gamma[i]/(math.factorial(i)*(1+pcr_yld)**(pcr_cycles)) #compute prob of m[i] muts
      mutDist = stats.rv_discrete(name='mutDist', values=(mut_m, prob_m)) #compute discrete mutation distribution
      print("Discrete Mutation Distribution has been computed")
# PCR Amplification
      for seq in slctd_seqs:
        mean_seq_num = slctd_seqs[seq][0]*(1+pcr_yld)**(pcr_cycles) #expectation of amplification
        var_seq_num = (1-pcr_yld)*(1+pcr_yld)**(pcr_cycles-1)*((1+pcr_yld)**(pcr_cycles)-1) #variance of amplification
        slctd_seqs[seq][0] = int(stats.norm.rvs(mean_seq_num, var_seq_num)) #draw int from Gaussian dist
        m = mutDist.rvs(size=slctd_seqs[seq][0]) #draw no. of mutations for each seq copy
        muts = m[m != 0] #keep only copies to be mutated (i.e. m >= 1)
        for mutNum, i in enumerate(muts): # for each mutation instance
          randPos = np.random.randint(seqLength, size=mutNum) #pick random nt positions for mutation  ##add seqLen as argumant to function
          if seq not in mutatedPool: #if seq never mutated before
            slctd_seqs[seq][0] -= 1 #decrement seq count of wild type
            mutatedPool.setdefault(seq, []).append(1) #add to mutated pool
            mutatedPool.setdefault(seq, []).append(randPos) #append mutation positions
          else: #if seq previously mutated
            slctd_seqs[seq][0] -= 1 #decrement no. of wild type seq
            mutatedPool[seq][0]+=1 #increment no. of seq to be mutated
            mutatedPool.setdefault(seq, []).append(randPos) #append mutation positions
            
      amplfd_seqs = slctd_seqs
      print("sequence amplification has been carried out")
      return amplfd_seqs, mutatedPool


# Assume all mutations happened at end of PCR (i.e. mutants are not amplified till next round)
## ADD seqLength as arg to function
## EDITED variance calculation for amplification
   def ampEffMutDefinite(self, slctd_seqs, seqLength, pcr_cycles, pcr_yld, errorRate):
      mutatedPool = {} #initialize pool of seqs to be mutated
      gamma = np.arange(seqLength) #initialize vector of gamma values 
      mut_m = np.arange(seqLength) #initialize vector of mut no. (0-seqLen)
      prob_m = np.arange(seqLength) #initialize corresponding vector of probs
# Mutation Distribution
      for i in mut_m: # for each possible no. of muts, compute corresponding probability
        for j in range(pcr_cycles): # for each pcr cycle
          gamma[i] += np.exp(-j*errorRate*seqLength)*math.factorial(pcr_cycles)/(math.factorial(j)*math.factorial(pcr_cycles - j))*(j)**(i)*pcr_yld**(j) #calculate sum (SEE PAPER)
          prob_m[i] = (errorRate*seqLength)**(i)*gamma[i]/(math.factorial(i)*(1+pcr_yld)**(pcr_cycles)) #compute prob of m[i] muts
      mutDist = stats.rv_discrete(name='mutDist', values=(mut_m, prob_m)) #compute discrete mutation distribution
      print("Discrete Mutation Distribution has been computed")
# PCR Amplification
      for seq in slctd_seqs:
        slctd_seqs[seq][0] = int(slctd_seqs[seq][0]*(1+pcr_yld)**(pcr_cycles)) #expectation of amplification
        m = mutDist.rvs(size=slctd_seqs[seq][0]) #draw no. of mutations for each seq copy
        muts = m[m != 0] #keep only copies to be mutated (i.e. m >= 1)
        for mutNum, i in enumerate(muts): # for each mutation instance
          randPos = np.random.randint(seqLength, size=mutNum) #pick random nt positions for mutation  ##add seqLen as argumant to function
          if seq not in mutatedPool: #if seq never mutated before
            slctd_seqs[seq][0] -= 1 #decrement seq count of wild type
            mutatedPool.setdefault(seq, []).append(1) #add to mutated pool
            mutatedPool.setdefault(seq, []).append(randPos) #append mutation positions
          else: #if seq previously mutated
            slctd_seqs[seq][0] -= 1 #decrement no. of wild type seq
            mutatedPool[seq][0]+=1 #increment no. of seq to be mutated
            mutatedPool.setdefault(seq, []).append(randPos) #append mutation positions
            
      amplfd_seqs = slctd_seqs
      print("sequence amplification has been carried out")
      return amplfd_seqs, mutatedPool




