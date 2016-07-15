import time
import random
import linecache
from itertools import izip, imap
import operator
from collections import OrderedDict
from scipy import stats
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
import math
from sklearn.mixture import GMM
from sklearn.mixture import VBGMM
from sklearn.mixture import DPGMM
from matplotlib import rcParams

#Allow LaTeX in plots
rcParams['text.usetex'] = True
rcParams['text.latex.preamble'] = r'\usepackage{amsmath}'


##PASTE THIS INTO TRAINING SET GENERATOR V6
## This script is used to generate a data set of frequency dists after PCR under different 
#   initial conditions. It then approximates each dist using a Gaussian Mixture model (GMM).
#   A histogram plot of the dist and GMM fitting is generated.
#   It follows the design of trainingSetGenerator.py but normalises
#   the gmm parameters based on expectation (see BruteGMMnormed). It also sorts the 
#   parameters of each gmm according to mean magnitudes



#Initiate class
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



# This method is intended to generate sequence population data after PCR
# This uses a brute-force method that models each individual cycle as Bernoulli test
# A histogram is generated from the data to represent a probability density function
# The pdf is then fitted with a Gaussian Mixture Model 
# Training carried out using expectation-maximization algorithm
############################################################################
# PARAMS:                                                                  #
# initialCount = the initial frequence of each sequence prior to PCR       #
# pcrCycles = the number of cycles to be performed prior to pdf generation #
# pcrYield = the average yield of PCR                                      #
# dataPoints = the initial number of unique sequences prior to PCR         #
# gaussNum = the max number of gaussians used to fit the pdf               #
# OUTPUT:                                                                  #
# bestModel = the best performing GMM model that approximated the pdf      #
# A plot of the hist of the actual population (green)  overlayed onto the  #
# GMM model (blue) and showing the individual gaussian components (red)    #
############################################################################
   def BruteGMM(self, initialCount, pcrCycles, pcrYield, dataPoints, gaussNum):
      amplfdSeqs = np.zeros((dataPoints, 1)) #create vector for seqs
      amplfdSeqsHist = np.zeros(dataPoints)
      for i in range(dataPoints):
        amplfdSeqs[i] = initialCount #assign init count to each seq
        amplfdSeqsHist[i] = initialCount
        for n in range(pcrCycles): #model each cycle as Bernoulli test
            amplfdSeqs[i] += np.random.binomial(amplfdSeqs[i], pcrYield) 

            amplfdSeqsHist[i] += np.random.binomial(amplfdSeqsHist[i], pcrYield) 
      
    #transform training data into appropriate dimensions   
      amplfdSeqs = amplfdSeqs.reshape(-1, 1)
      N = np.arange(1, gaussNum)

    #generate models with varying number of gaussian components (i.e. from 1 to N gaussians)  
      gmmModels = [None for i in range(len(N))]

    #train each GMM model on the generated data
      for i in range(len(N)):
          gmmModels[i] = GMM(n_components=N[i], n_init=10, n_iter=100).fit(amplfdSeqs)
      
    #calculate AIC and BIC for each trained model
      gmmAIC = [m.aic(amplfdSeqs) for m in gmmModels]
      gmmBIC = [m.bic(amplfdSeqs) for m in gmmModels]


    #pick best trained model based on AIC
      bestModel = gmmModels[np.argmin(gmmAIC)]


      
    #calculate exact solutions to the expectation and variance of the dist
      expctdSeqNum = initialCount*(1+pcrYield)**(pcrCycles)
      varSeqNum = pcrYield*(1-pcrYield)*(initialCount*(1+pcrYield)**(pcrCycles - 1))*((initialCount*(1+pcrYield)**(pcrCycles))-1)/(initialCount*(1+pcrYield)-1)

      dataRange =  expctdSeqNum + np.sqrt(varSeqNum)
    #declare sample space 
      space = np.linspace(1, (expctdSeqNum + (np.sqrt(varSeqNum)*10)), dataPoints/100)
      space = space.reshape(-1, 1)
    #calculate log likelihood of sample space using trained model
      logProbs, resps = bestModel.score_samples(space)
    #calculate prob density func of the model
      pdf = np.exp(logProbs)
    #calculate prob density func of each component gaussian
      individualPDFs = resps * pdf[:, np.newaxis]

     # this was added for histogram visualization purposes
      binsNum = dataPoints*5/1000
      weights = np.ones_like(amplfdSeqsHist)/dataPoints
      
      fig = plt.figure()
      ax = fig.add_subplot(111)


      ax.hist(amplfdSeqsHist, bins=50, normed=True, facecolor='green', alpha=0.75)
      ax.plot(space, pdf, '-k', color='b')
      ax.plot(space, individualPDFs, '--k', color='r')
      ax.set_xlabel('Sequence Count')
      ax.set_ylabel('p(x)')
      ax.set_title('GMM Best-fit to Population Distribution')
      # create annotation
      annot = " \mu & \sigma^{2} & \omega \\"+"\\"
      for i, mu in enumerate(bestModel.means_):
          annot += str(np.round(bestModel.means_[i][0], 2))+" & "+str(np.round(bestModel.covars_[i][0], 2))+" & "+str(np.round(bestModel.weights_[i], 2))+" \\"+"\\ "
      
      #add plot annotations
      ax.text(0.95, 0.95, r"Initial count = "+str(initialCount)+'\n'+"No. of cycles = "+str(pcrCycles)+'\n'+"Yield = "+str(pcrYield), verticalalignment='top', horizontalalignment='right', transform=ax.transAxes, color='black', fontsize=10)
      ax.text(0.935, 0.65, r"$ \begin{pmatrix} %s  \end{pmatrix}$" % annot, verticalalignment='top', horizontalalignment='right', transform=ax.transAxes, color='black', fontsize=10)


      # save plot
      plt.grid()
      plt.savefig('pcrDistEst_n'+str(pcrCycles)+'_i'+str(initialCount)+'_y'+str(pcrYield)+'.pdf', format='pdf')
      #plt.close()
      #plt.show()
      return bestModel

# TEST AREA - TO BE DELETED
#amp = Amplification()
#amp.BruteGMM(3, 5, 0.85, 10000, 20)


   def BruteGMMnormed(self, initialCount, pcrCycles, pcrYield, dataPoints, gaussNum):
      amplfdSeqs = np.zeros((dataPoints, 1)) #create vector for seqs
      amplfdSeqsHist = np.zeros(dataPoints)

    
    #calculate exact solutions to the expectation and variance of the dist
      expctdSeqNum = initialCount*(1+pcrYield)**(pcrCycles)
      varSeqNum = pcrYield*(1-pcrYield)*(initialCount*(1+pcrYield)**(pcrCycles - 1))*((initialCount*(1+pcrYield)**(pcrCycles))-1)/(initialCount*(1+pcrYield)-1)


      for i in range(dataPoints):
        amplfdSeqs[i] = initialCount #assign init count to each seq
        amplfdSeqsHist[i] = initialCount
        for n in range(pcrCycles): #model each cycle as Bernoulli test
            amplfdSeqs[i] += np.random.binomial(amplfdSeqs[i], pcrYield) 

            amplfdSeqsHist[i] += np.random.binomial(amplfdSeqsHist[i], pcrYield) 
      
      for i in range(dataPoints):
          amplfdSeqs[i] = amplfdSeqs[i]/(2*expctdSeqNum)

    #transform training data into appropriate dimensions   
      amplfdSeqs = amplfdSeqs.reshape(-1, 1)
      #N = np.arange(1, gaussNum)

    #generate models with varying number of gaussian components (i.e. from 1 to N gaussians)  
      #gmmModels = [None for i in range(len(N))]

    #train each GMM model on the generated data
      #for i in range(len(N)):
      gmmModel = GMM(n_components=gaussNum, n_init=10, n_iter=100).fit(amplfdSeqs)

      modelParams = np.zeros((gaussNum, 3))

      for i, mu in enumerate(gmmModel.means_):
          modelParams[i][0] = gmmModel.means_[i][0]
          modelParams[i][1] = gmmModel.covars_[i][0]
          modelParams[i][2] = gmmModel.weights_[i]

      modelParams = modelParams[modelParams[:, 0].argsort()]
    #calculate AIC and BIC for each trained model
      gmmAIC = gmmModel.aic(amplfdSeqs)
      gmmBIC = gmmModel.bic(amplfdSeqs)


    #pick best trained model based on AIC
      #bestModel = gmmModels[np.argmin(gmmAIC)]


      
      dataRange =  expctdSeqNum + np.sqrt(varSeqNum)
    #declare sample space 
      space = np.linspace(1, (expctdSeqNum + (np.sqrt(varSeqNum)*10)), dataPoints/100)
      normedSpace = np.linspace((1/(2*expctdSeqNum)), ((expctdSeqNum + (np.sqrt(varSeqNum)*10))/(2*expctdSeqNum)), dataPoints/100)
      space = space.reshape(-1, 1)
      normedSpace = normedSpace.reshape(-1, 1)
    #calculate log likelihood of sample space using trained model
      logProbs, resps = gmmModel.score_samples(normedSpace)
    #calculate prob density func of the model
      pdf = np.exp(logProbs)/(2*expctdSeqNum)
    #calculate prob density func of each component gaussian
      individualPDFs = resps * pdf[:, np.newaxis]

     # this was added for histogram visualization purposes
      #binsNum = dataPoints*5/1000
      #weights = np.ones_like(amplfdSeqsHist)/dataPoints
      
      #fig = plt.figure()
      #ax = fig.add_subplot(111)


      #ax.hist(amplfdSeqsHist, bins=50, normed=True, facecolor='green', alpha=0.75)
      #ax.plot(space, pdf, '-k', color='b')
      #ax.plot(space, individualPDFs, '--k', color='r')
      #ax.set_xlabel('Sequence Count')
      #ax.set_ylabel('p(x)')
      #ax.set_title('GMM Best-fit to Population Distribution')
      # create annotation
      #annot = " \mu & \sigma^{2} & \omega \\"+"\\"
      #for i, mu in enumerate(gmmModel.means_):
          #annot += str(np.round(gmmModel.means_[i][0], 2))+" & "+str(np.round(gmmModel.covars_[i][0], 2))+" & "+str(np.round(gmmModel.weights_[i], 2))+" \\"+"\\ "
      
      #add plot annotations
      #ax.text(0.95, 0.95, r"Initial count = "+str(initialCount)+'\n'+"No. of cycles = "+str(pcrCycles)+'\n'+"Yield = "+str(pcrYield), verticalalignment='top', horizontalalignment='right', transform=ax.transAxes, color='black', fontsize=10)
      #ax.text(0.935, 0.65, r"$ \begin{pmatrix} %s  \end{pmatrix}$" % annot, verticalalignment='top', horizontalalignment='right', transform=ax.transAxes, color='black', fontsize=10)


      # save plot
      #plt.grid()
      #plt.savefig('pcrDistEst_n'+str(pcrCycles)+'_i'+str(initialCount)+'_y'+str(pcrYield)+'.pdf', format='pdf')
      #plt.close()
      #plt.show()
      return modelParams

# TEST AREA - TO BE DELETED
#amp = Amplification()
#gmm, pdf, params = amp.BruteGMMnormed(3, 5, 0.85, 10000, 8)



   def GMMTest(self, initialCount, pcrCycles, pcrYield, dataPoints, gaussNum):
      amplfdSeqs = np.zeros((dataPoints, 1)) #create vector for seqs
      amplfdSeqsHist = np.zeros(dataPoints)
      for i in range(dataPoints):
        amplfdSeqs[i] = initialCount #assign init count to each seq
        amplfdSeqsHist[i] = initialCount
        for n in range(pcrCycles): #model each cycle as Bernoulli test
            amplfdSeqs[i] += np.random.binomial(amplfdSeqs[i], pcrYield) 

            amplfdSeqsHist[i] += np.random.binomial(amplfdSeqsHist[i], pcrYield) 
      
    #transform training data into appropriate dimensions   
      amplfdSeqs = amplfdSeqs.reshape(-1, 1)
      #N = np.arange(1, gaussNum)

    #generate models with varying number of gaussian components (i.e. from 1 to N gaussians)  
      #gmmModels = [None for i in range(len(N))]

    #train each GMM model on the generated data
      #for i in range(len(N)):
      gmmModel = GMM(n_components=gaussNum, n_init=10, n_iter=100).fit(amplfdSeqs)
      
    #calculate AIC and BIC for each trained model
      gmmAIC = gmmModel.aic(amplfdSeqs)
      gmmBIC = gmmModel.bic(amplfdSeqs)


    #pick best trained model based on AIC
      #bestModel = gmmModels[np.argmin(gmmAIC)]


      
    #calculate exact solutions to the expectation and variance of the dist
      expctdSeqNum = initialCount*(1+pcrYield)**(pcrCycles)
      varSeqNum = pcrYield*(1-pcrYield)*(initialCount*(1+pcrYield)**(pcrCycles - 1))*((initialCount*(1+pcrYield)**(pcrCycles))-1)/(initialCount*(1+pcrYield)-1)

      dataRange =  expctdSeqNum + np.sqrt(varSeqNum)
    #declare sample space 
      space = np.linspace(1, (expctdSeqNum + (np.sqrt(varSeqNum)*10)), dataPoints/100)
      space = space.reshape(-1, 1)
    #calculate log likelihood of sample space using trained model
      logProbs, resps = gmmModel.score_samples(space)
    #calculate prob density func of the model
      pdf = np.exp(logProbs)
    #calculate prob density func of each component gaussian
      individualPDFs = resps * pdf[:, np.newaxis]

     # this was added for histogram visualization purposes
      binsNum = dataPoints*5/1000
      weights = np.ones_like(amplfdSeqsHist)/dataPoints
      
      fig = plt.figure()
      ax = fig.add_subplot(111)


      ax.hist(amplfdSeqsHist, bins=50, normed=True, facecolor='green', alpha=0.75)
      ax.plot(space, pdf, '-k', color='b')
      ax.plot(space, individualPDFs, '--k', color='r')
      ax.set_xlabel('Sequence Count')
      ax.set_ylabel('p(x)')
      ax.set_title('GMM Best-fit to Population Distribution')
      # create annotation
      annot = " \mu & \sigma^{2} & \omega \\"+"\\"
      for i, mu in enumerate(gmmModel.means_):
          annot += str(np.round(gmmModel.means_[i][0], 2))+" & "+str(np.round(gmmModel.covars_[i][0], 2))+" & "+str(np.round(gmmModel.weights_[i], 2))+" \\"+"\\ "
      
      #add plot annotations
      ax.text(0.95, 0.95, r"Initial count = "+str(initialCount)+'\n'+"No. of cycles = "+str(pcrCycles)+'\n'+"Yield = "+str(pcrYield), verticalalignment='top', horizontalalignment='right', transform=ax.transAxes, color='black', fontsize=10)
      ax.text(0.935, 0.65, r"$ \begin{pmatrix} %s  \end{pmatrix}$" % annot, verticalalignment='top', horizontalalignment='right', transform=ax.transAxes, color='black', fontsize=10)


      # save plot
      plt.grid()
      plt.savefig('pcrDistEst_n'+str(pcrCycles)+'_i'+str(initialCount)+'_y'+str(pcrYield)+'.pdf', format='pdf')
      #plt.close()
      #plt.show()
      return gmmModel

# TEST AREA - TO BE DELETED
#amp = Amplification()
#amp.GMMTest(3, 5, 0.85, 10000, 8)






   def bruteHist(self, initialCount, pcrCycles, pcrYield, dataPoints, binsNum):
      amplfdSeqs = np.zeros((dataPoints, 1)) #create vector for seqs
      amplfdSeqsHist = np.zeros(dataPoints)
      for i in range(dataPoints):
        amplfdSeqs[i] = initialCount #assign init count to each seq
        amplfdSeqsHist[i] = initialCount
        for n in range(pcrCycles): #model each cycle as Bernoulli test
            amplfdSeqs[i] += np.random.binomial(amplfdSeqs[i], pcrYield) 

            amplfdSeqsHist[i] += np.random.binomial(amplfdSeqsHist[i], pcrYield) 
      
      n, bins, patches = plt.hist(amplfdSeqsHist, binsNum, normed=True, facecolor='green', alpha=0.75)
      #plt.plot(space, pdf, '-k', color='b')
      #plt.plot(space, individualPDFs, '--k')
      plt.xlabel('Sequence Count')
      plt.ylabel('p(x)')
      plt.title('GMM Best-fit')
      plt.savefig('pcrDistEst_n'+str(pcrCycles)+'_i'+str(initialCount)+'_y'+str(pcrYield)+'.pdf', format='pdf')
      plt.show()
# TEST AREA - TO BE DELETED
#amp = Amplification()
#amp.BruteHist(2, 20, 0.65, 10000, 20)




   def ampEffBruteTest(self, initialCount, pcrCycles, pcrYield, dataPoints):
      amplfdSeqs = np.zeros(dataPoints)
      for i in range(dataPoints):
        amplfdSeqs[i] = initialCount
        for n in range(pcrCycles):
            amplfdSeqs[i] += np.random.binomial(amplfdSeqs[i], pcrYield)
      n, bins, patches = plt.hist(amplfdSeqs, 10000, normed=1, facecolor='green', alpha=0.75)
      plt.xlabel('Sequence Count')
      plt.ylabel('Probability')
      plt.title('A histogram plot of '+str(dataPoints)+' sequences after '+str(pcrCycles)+' PCR ciycles with '+str(pcrYield)+'% yield.\n'+'Initial count = '+str(initialCount))
      plt.savefig('pcrDistEst_n'+str(pcrCycles)+'_i'+str(initialCount)+'_y'+str(pcrYield)+'.png')
      plt.close()
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
      mut_m = np.arange(seqLength) #initialize vector of mut no.s (0-seqLen)
      prob_m = np.arange(seqLength) #initialize corresponding vector of probs
      # Generate distribution for generation number        
      prob_n = np.zeros(pcr_cycles+1)
      genNum = np.arange(pcr_cycles+1)
      for n in range(pcr_cycles+1): #for each generation
        prob_n[n] = math.factorial(pcr_cycles)*(pcrYield)**(n)/(math.factorial(n)*math.factorial(pcr_cycles - n)*(1+pcrYield)**(pcr_cycles)) #compute prob of picking seq from gen n
      genDist = stats.rv_discrete(name='genDist', values=(genNum, prob_n)) #construct gen Distribution
      print("Discrete Generation Number Distribution has been computed")
# Mutation Distribution
      for i in mut_m: # for each possible no. of muts, compute corresponding probability
        for j in range(pcr_cycles+1): # for each pcr generation
          gamma[i] += np.exp(-j*errorRate*seqLength)*math.factorial(pcr_cycles)/(math.factorial(j)*math.factorial(pcr_cycles - j))*(j)**(i)*pcr_yld**(j) #calculate sum (SEE PAPER)
        prob_m[i] = (errorRate*seqLength)**(i)*gamma[i]/(math.factorial(i)*(1+pcr_yld)**(pcr_cycles)) #compute prob of m[i] muts
      mutDist = stats.rv_discrete(name='mutDist', values=(mut_m, prob_m)) #construct discrete mutation distribution
      print("Discrete Mutation Distribution has been computed")
# PCR Amplification
      for seq in slctd_seqs: #for each seq in selected pool
        slctd_seqs[seq][0] = int(slctd_seqs[seq][0]*(1+pcr_yld)**(pcr_cycles)) #expectation of amplification
        m = mutDist.rvs(size=slctd_seqs[seq][0]) #draw no. of mutations for each seq copy
        muts = m[m != 0] #keep only copies to be mutated (i.e. m >= 1)
        g = genDist.rvs(muts.shape) #draw gen no. for each mut'd seq copy
        for i, mutNum in enumerate(muts): # for each mutation instance
          randPos = np.random.randint(seqLength, size=mutNum) #pick random nt positions for mutation  ##add seqLen as argumant to function
          if seq not in mutatedPool: #if seq never mutated before
            slctd_seqs[seq][0] -= (1+pcr_yld)**(pcr_cycles-g[i]) #decrement seq count of wild type
            mutatedPool.setdefault(seq, []).append(np.zeros(muts.shape)) # add number of copies to be mutated
            mutatedPool[seq][0][i] = ((1+pcr_yld)**(pcr_cycles-g[i])) #add expected amplfd count 
            mutatedPool.setdefault(seq, []).append(randPos) #append mutation positions
          else: #if seq previously mutated
            slctd_seqs[seq][0] -= (1+pcr_yld)**(pcr_cycles-g[i]) #decrement no. of wild type seq
            mutatedPool[seq][0][i] = (1+pcr_yld)**(pcr_cycles-g[i]) #increment no. of seq to be mutated
            mutatedPool.setdefault(seq, []).append(randPos) #append mutation positions
            
      amplfd_seqs = slctd_seqs
      print("sequence amplification has been carried out")
      return amplfd_seqs, mutatedPool



   def BruteVBGMM(self, initialCount, pcrCycles, pcrYield, dataPoints, gaussNum, a, r):
      amplfdSeqs = np.zeros((dataPoints, 1)) #create vector for seqs
      amplfdSeqsHist = np.zeros(dataPoints)
      for i in range(dataPoints):
        amplfdSeqs[i] = initialCount #assign init count to each seq
        amplfdSeqsHist[i] = initialCount
        for n in range(pcrCycles): #model each cycle as Bernoulli test
            amplfdSeqs[i] += np.random.binomial(amplfdSeqs[i], pcrYield) 

            amplfdSeqsHist[i] += np.random.binomial(amplfdSeqsHist[i], pcrYield) 
      
      amplfdSeqs = amplfdSeqs.reshape(-1, 1)
      N = np.arange(1, gaussNum)
      gmmModels = [None for i in range(len(N))]

      for i in range(len(N)):
          gmmModels[i] = VBGMM(N[i], n_iter=1000, alpha=a, random_state=r).fit(amplfdSeqs)
      

      gmmAIC = [m.aic(amplfdSeqs) for m in gmmModels]
      gmmBIC = [m.bic(amplfdSeqs) for m in gmmModels]


      bestModel = gmmModels[np.argmin(gmmAIC)]
      

      expctdSeqNum = initialCount*(1+pcrYield)**(pcrCycles)
      varSeqNum = pcrYield*(1-pcrYield)*(initialCount*(1+pcrYield)**(pcrCycles - 1))*((initialCount*(1+pcrYield)**(pcrCycles))-1)/(initialCount*(1+pcrYield)-1)


      space = np.linspace(1, (expctdSeqNum + (np.sqrt(varSeqNum)*7)), dataPoints/100)
      space = space.reshape(-1, 1)
      logProbs, resps = bestModel.score_samples(space)

      pdf = np.exp(logProbs)
      individualPDFs = resps * pdf[:, np.newaxis]

      n, bins, patches = plt.hist(amplfdSeqsHist, (dataPoints/100), normed=True, facecolor='green', alpha=0.75)
      plt.plot(space, pdf, '-k', color='b')
      plt.plot(space, individualPDFs, '--k')
      plt.xlabel('Sequence Count')
      plt.ylabel('p(x)')
      plt.title('GMM Best-fit')
      plt.show()
      return bestModel

# TEST AREA - TO BE DELETED
#amp = Amplification()
#amp.BruteVBGMM(2, 20, 0.65, 10000, 20, 0.2, 1)




# THIS DOES NOT YET WORK PROPERLY
   def BruteDPGMM(self, initialCount, pcrCycles, pcrYield, dataPoints, gaussNum):
      amplfdSeqs = np.zeros((dataPoints, 1)) #create vector for seqs
      for i in range(dataPoints):
        amplfdSeqs[i] = initialCount #assign init count to each seq
        for n in range(pcrCycles): #model each cycle as Bernoulli test
            amplfdSeqs[i] += np.random.binomial(amplfdSeqs[i], pcrYield) 
      
      amplfdSeqs = amplfdSeqs.reshape(-1, 1)

      
      dpgmmModels = DPGMM(n_components=gaussNum, covariance_type='full').fit(amplfdSeqs)
      

      dpgmmAIC = dpgmmModels.aic(amplfdSeqs)
      dpgmmBIC = dpgmmModels.bic(amplfdSeqs)


      bestModel = dpgmmModels
      

      expctdSeqNum = initialCount*(1+pcrYield)**(pcrCycles)
      varSeqNum = pcrYield*(1-pcrYield)*(initialCount*(1+pcrYield)**(pcrCycles - 1))*((initialCount*(1+pcrYield)**(pcrCycles))-1)/(initialCount*(1+pcrYield)-1)


      space = np.linspace(1, (expctdSeqNum + (np.sqrt(varSeqNum)*7)), 1000)
      space = space.reshape(-1, 1)

      logProbs, weights = bestModel.score_samples(space)

      pdf = np.exp(logProbs)
      individualPDFs = weights * pdf[:, np.newaxis]

      plt.plot(space, pdf, '-k', color='b')
      plt.plot(space, individualPDFs, '--k')
      plt.xlabel('Sequence Count')
      plt.ylabel('p(x)')
      plt.title('GMM Best-fit')
      plt.show()
      return bestModel

# TEST AREA - TO BE DELETED
#amp = Amplification()
#amp.BruteDPGMM(2, 20, 0.65, 10000, 20)


