import time
import gc
from scipy import stats
from scipy.stats import norm
import numpy as np
from numpy.random import binomial as binom
from numpy.random import poisson
import matplotlib
matplotlib.use('Agg')
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
import math
from sklearn.mixture import GMM
from sklearn.mixture import VBGMM
from sklearn.mixture import DPGMM
from matplotlib import rcParams
from Mutation import Mutation

#Allow LaTeX in plots
rcParams['text.usetex'] = True
rcParams['text.latex.preamble'] = r'\usepackage{amsmath}'
#Turn off interactive mode for plots
plt.ioff()


#Initiate class
class Amplification:
    # add __init__ contructor here

    def randomPCR_with_ErrorsAndBias(self, slctdSeqs, 
                                     seqLength, pcrCycleNum, 
                                     pcrYld, errorRate, 
                                     aptamerSeqs, alphabetSet):
        # initialize Mutation object from class
        mut = Mutation(seqLength=seqLength, errorRate=errorRate, 
                        pcrCycleNum=pcrCycleNum, pcrYld=pcrYld)
        # count number of seqs in selected pool
        totalseqs = 0
        uniqSeqs = 0
        for i, seqIdx in enumerate(slctdSeqs):
            uniqSeqs += 1
            totalseqs += slctdSeqs[seqIdx][0]
        print("number of unique seqs in selected pool prior to amplification: "+str(uniqSeqs))
        print("number of seqs in selected pool prior to amplification: "+str(totalseqs))
        # calculate probabilities of different possible mutation numbers
        mutNumProbs = mut.get_mutation_probabilities_original()
        print mutNumProbs
        mutDist = mut.get_mutation_distribution_original()
        print mutDist.rvs(size=1000000)
        print("Discrete Mutation Distribution has been computed")
    # PCR Amplification
        totalseqs = 0
        # initialize dictionary to keep info on seqs to be mutated
        mutatedPool = {}
        # keep track of sequence count after each pcr cycle (except last one)
        seqPop = np.zeros(pcrCycleNum)
        print("Amplification has started...")
        # for each sequence in the selected pool
        for i, seqIdx in enumerate(slctdSeqs):
            # random PCR with bias using brute force        
            for n in xrange(pcrCycleNum):
                # sequence count after n cycles
                seqPop[n] = slctdSeqs[seqIdx][0]
                # amplify count using initial count, polymerase yield, and bias score
                slctdSeqs[seqIdx][0] += int(binom(slctdSeqs[seqIdx][0], 
                                             pcrYld + slctdSeqs[seqIdx][2]))
            # compute cycle number probabilities for this seq
            cycleNumProbs = mut.get_cycleNumber_probabilities(seqPop=seqPop)
            # append the information to the dictionary
            slctdSeqs.setdefault(seqIdx, []).append(cycleNumProbs) #index [3]
            # update total num of seqs
            totalseqs += slctdSeqs[seqIdx][0]
            if slctdSeqs[seqIdx][0] > 10000:
            # Note{The following should give accurate results only if the initial
            # seq count (slctdSeqs[seqIdx][0] is sufficiently large). If not, 
            # you have to implement an alternative method that relies on
            # random numbers :( }
            # for each possible number of mutations in one instance (1-seqLength)
                #print("seq num is greater than 10000")
                mutatedPool.setdefault(seqIdx, []).append(0)
                for mutNum in xrange(seqLength):
                    mutFreq = int(mutNumProbs[mutNum+1]*slctdSeqs[seqIdx][0])
                    mutatedPool[seqIdx][0] += mutFreq
                    mutatedPool.setdefault(seqIdx, []).append(mutFreq)
            else:
                muts = mutDist.rvs(size=slctdSeqs[seqIdx][0])
                muts = muts[muts != 0]
                if seqIdx in mutatedPool:
                    for mutNum in muts:
                        mutatedPool[seqIdx][mutNum] += 1
                else:
                    mutatedPool.setdefault(seqIdx, []).append(0)
                    for mutNum in xrange(seqLength):
                        mutatedPool.setdefault(seqIdx, []).append(0)
                    for mutNum in muts:
                        mutatedPool[seqIdx][mutNum] += 1
        print("Amplification has been carried out")
        print("Mutation has started...")
        # generate mutants and add to the amplfied sequence pool 
        amplfdSeqs = mut.generate_mutants(mutatedPool=mutatedPool, 
                                            amplfdSeqs=slctdSeqs, 
                                            aptamerSeqs=aptamerSeqs, 
                                            alphabetSet=alphabetSet)
        print("sequence amplification has been carried out")
        return amplfdSeqs

    def randomPCR_with_ErrorsAndBias_FAST(self, slctdSeqs, 
                                     seqLength, pcrCycleNum, 
                                     pcrYld, errorRate, 
                                     aptamerSeqs, alphabetSet):
        # initialize Mutation object from class
        mut = Mutation(seqLength=seqLength, errorRate=errorRate, 
                        pcrCycleNum=pcrCycleNum, pcrYld=pcrYld)
        # count number of seqs in selected pool
        totalseqs = 0
        uniqSeqs = 0
        for i, seqIdx in enumerate(slctdSeqs):
            uniqSeqs += 1
            totalseqs += slctdSeqs[seqIdx][0]
        print("number of unique seqs in selected pool prior to amplification: "+str(uniqSeqs))
        print("number of seqs in selected pool prior to amplification: "+str(totalseqs))
        # calculate probabilities of different possible mutation numbers
        mutNumProbs = mut.get_mutation_probabilities_original()
        # compute a discrete distribution of mutation numbers
        mutDist = mut.get_mutation_distribution_original()
        print("Discrete Mutation Distribution has been computed")
    # PCR Amplification
        totalseqs = 0
        # initialize dictionary to keep info on seqs to be mutated
        mutatedPool = {}
        # keep track of sequence count after each pcr cycle (except last one)
        seqPop = np.zeros(pcrCycleNum)
        print("Amplification has started...")
        #initialize matrix to hold info for amplified pool
        x = np.zeros((uniqSeqs, 2))
        # for each sequence in the selected pool
        for i, seqIdx in enumerate(slctdSeqs):
            # random PCR with bias using brute force        
            for n in xrange(pcrCycleNum):
                # sequence count after n cycles
                seqPop[n] = slctdSeqs[seqIdx][0]
                # amplify count using initial count, polymerase yield, and bias score
                slctdSeqs[seqIdx][0] += int(binom(slctdSeqs[seqIdx][0], 
                                             pcrYld + slctdSeqs[seqIdx][2]))
            # compute cycle number probabilities for this seq
            cycleNumProbs = mut.get_cycleNumber_probabilities(seqPop=seqPop)
            # append the information to the dictionary
            slctdSeqs.setdefault(seqIdx, []).append(cycleNumProbs) #index [3]
            x[i][0] = seqIdx
            x[i][1] = slctdSeqs[seqIdx][0]
            # update total num of seqs
            totalseqs += slctdSeqs[seqIdx][0]
        #initialize matrix to hold info for mutation pool
        y = np.zeros((uniqSeqs, seqLength+1))
        #for each sequence in matrix x
        for j, seqInfo in enumerate(x):
            #tranfer seq index to matrix y
            y[j][0] = seqInfo[0]
            #if seq count is greater than 10,000
            if seqInfo[1] > 10000:

            # Note{The following should give accurate results only if the initial
            # seq count (slctdSeqs[seqIdx][0] is sufficiently large). If not, 
            # you have to implement an alternative method that relies on
            # random numbers :( }

            # for each possible number of mutations in any seq copy (1-seqLength)
                for mutNum in xrange(seqLength):
                    #approximate the proportion of copies that will be mutated using
                    #corresponding probability p(M=mutNum)
                    y[j][mutNum+1] = int(mutNumProbs[mutNum+1]*slctdSeqs[seqIdx][0])
            # if seq count is less than 10,000
            else:
                #draw random mutNum from the mutation distribution for each seq copy
                muts = mutDist.rvs(size=slctdSeqs[seqIdx][0])
                #remove all drawn numbers equal to zero
                muts = muts[muts != 0]
                #for each non-zero mutation number
                for mutNum in muts:
                    #increment copy number to be mutated
                    y[j][mutNum+1] += 1
        #remove all mutation numbers with zero copies to be mutated
        y = y[y[:, 1] != 0]
        #for each seq to be mutated
        for mutInfo in y:
            #add to mutation pool with it's corresponding mutation info
            mutatedPool[mutInfo[0]] = mutInfo[1:][mutInfo[1:] != 0]
        print("Amplification has been carried out")
        print("Mutation has started...")
        # generate mutants and add to the amplfied sequence pool 
        amplfdSeqs = mut.generate_mutants_FAST(mutatedPool=mutatedPool, 
                                            amplfdSeqs=slctdSeqs, 
                                            aptamerSeqs=aptamerSeqs, 
                                            alphabetSet=alphabetSet)
        print("sequence amplification has been carried out")
        return amplfdSeqs

    def randomPCR_with_ErrorsAndBias_FASTv2(self, slctdSeqs, 
                                     seqLength, pcrCycleNum, 
                                     pcrYld, errorRate, 
                                     aptamerSeqs, alphabetSet, distance):
        # initialize Mutation object from class
        mut = Mutation(seqLength=seqLength, errorRate=errorRate, 
                        pcrCycleNum=pcrCycleNum, pcrYld=pcrYld)
        # count number of seqs in selected pool
        totalseqs = 0
        uniqSeqs = 0
        #compute total seq num, unique seq num, and transfer info to x
        for i, seqIdx in enumerate(slctdSeqs):
            uniqSeqs += 1
            totalseqs += int(slctdSeqs[seqIdx][0])
        #initialize matrix to hold info for amplified pool
        x = np.zeros((uniqSeqs, pcrCycleNum+4))
        for i, seqIdx in enumerate(slctdSeqs):
            x[i][0] = seqIdx
            x[i][1] = slctdSeqs[seqIdx][0]
            x[i][2] = slctdSeqs[seqIdx][1]
            x[i][3] = slctdSeqs[seqIdx][2]
        print("number of unique seqs in selected pool prior to amplification: "+str(uniqSeqs))
        print("number of seqs in selected pool prior to amplification: "+str(totalseqs))
        # calculate probabilities of different possible mutation numbers
        mutNumProbs = mut.get_mutation_probabilities_original()
        # compute a discrete distribution of mutation numbers
        mutDist = mut.get_mutation_distribution_original()
        print("Discrete Mutation Distribution has been computed")
    # PCR Amplification
        totalseqs = 0
        # initialize dictionary to keep info on seqs to be mutated
        mutatedPool = {}
        #initialize matrix to hold info for mutation pool
        y = np.zeros((uniqSeqs, seqLength+1))
        # keep track of sequence count after each pcr cycle (except last one)
        seqPop = np.zeros(pcrCycleNum)
        # compute cycle number probabilities for this seq
        cycleNumProbs = np.zeros(pcrCycleNum)
        print("Amplification has started...")
        # for each sequence in the selected pool
        for i, seqIdx in enumerate(slctdSeqs):
            # random PCR with bias using brute force        
            for n in xrange(pcrCycleNum):
                # sequence count after n cycles
                seqPop[n] = x[i][1]
                # amplify count using initial count, polymerase yield, and bias score
                x[i][1] += int(binom(x[i][1], pcrYld+x[i][3]))
            # compute cycle number probabilities
            for s, seqNum in enumerate(seqPop):
                cycleNumProbs[s] = seqNum/np.sum(seqPop)
            # transfer info to x
            for j, cycleNumProb in enumerate(cycleNumProbs):
                x[i][j+4] = cycleNumProb
            # update total num of seqs
            totalseqs += x[i][1]
            # transfer info from x to selection pool
            slctdSeqs[int(x[i][0])] = x[i][1:]
            #tranfer seq index to matrix y
            y[i][0] = x[i][0]
            #if accumulated seq count is greater than 10,000
            if np.sum(seqPop) > 10000:
            # for each possible number of mutations in any seq copy (1-seqLength)
                for mutNum in xrange(seqLength):
                    #approximate the proportion of copies that will be mutated using
                    #corresponding probability p(M=mutNum)
                    y[i][mutNum+1] = mutNumProbs[mutNum+1]*np.sum(seqPop)
            # if seq count is less than 10,000
            else:
                # draw random mutNum from the mutation distribution for each seq copy
                muts = poisson(errorRate*seqLength, int(np.sum(seqPop))) #SLOW STEP
                # remove all drawn numbers equal to zero
                muts = muts[muts != 0]
                # for each non-zero mutation number
                for mutNum in muts:
                    #increment copy number to be mutated
                    y[i][mutNum+1] += 1
        print("Amplification carried out")
        print("Sequence selection for mutation has started...")
        #remove all mutation numbers with zero copies to be mutated
        y = y[y[:, 1] != 0]
        #for each seq to be mutated
        for mutInfo in y:
            #add to mutation pool with it's corresponding mutation info
            mutatedPool[int(mutInfo[0])] = mutInfo[1:][mutInfo[1:] != 0]
        print("Mutation selection has been carried out")
        print("Mutant generation has started...")
        if(distance == "hamming"):
            # generate mutants and add to the amplfied sequence pool 
            amplfdSeqs = mut.generate_mutants_FASTv2(mutatedPool=mutatedPool, 
                                                     amplfdSeqs=slctdSeqs, 
                                                     aptamerSeqs=aptamerSeqs, 
                                                     alphabetSet=alphabetSet)
        elif(distance == "basepair"):
            amplfdSeqs = mut.generate_mutants_2D(mutatedPool=mutatedPool,
                                                 amplfdSeqs=slctdSeqs,
                                                 aptamerSeqs=aptamerSeqs,
                                                 alphabetSet=alphabetSet)
        else:
            print("argument given for distance is invalid")
            return
        print("Mutation has been carried out")
        return amplfdSeqs




#NEED TO INCLUDE GENERATION NUMBER AND ACCOUNT FOR
# PARTIAL AMPLIFICATIONS
# This method simulates PCR under stochastic effects with Mutations  amplification
# Amplifications done using brute-force approach
# Returns amplfdSeqs and mutatedPool dicts
    def randomPCR_errorProne_biased(self, slctdSeqs, seqLength, pcrCycleNum, pcrYld, errorRate):
        totalseqs = 0
        for i, seqIdx in enumerate(slctdSeqs):
            totalseqs += 1
        print totalseqs
        mutatedPool = {} #initialize pool of seqs to be mutated
        gamma = np.arange(seqLength) #initialize vector of gamma values 
        mut_m = np.arange(seqLength) #initialize vector of mut no. (0-seqLen)
        prob_m = np.zeros(seqLength) #initialize corresponding vector of probs
# Mutation Distribution #gamma vector not used
        for i in mut_m: # for each possible no. of muts, compute corresponding probability
            const = (errorRate*seqLength)**(i)/(math.factorial(i)*(1+pcrYld)**(pcrCycleNum)) #compute prob of m[i] muts
            for j in range(pcrCycleNum): # for each pcr cycle
                prob_m[i] += const*np.exp(-j*errorRate*seqLength)*math.factorial(pcrCycleNum)/(math.factorial(j)*math.factorial(pcrCycleNum - j))*(j)**(i)*pcrYld**(j) #calculate sum (SEE PAPER)
        
        mutDist = stats.rv_discrete(name='mutDist', values=(mut_m, prob_m)) #compute discrete mutation number distribution
        print mut_m
        print prob_m
        print mutDist
        print("Discrete Mutation Distribution has been computed")
        totalseqs = 0
# PCR Amplification
#EDIT MUTATIONS SUCH THAT THEY ARE CARRIED OUT RANDOMLY ON SLCTD POOL INSTEAD OF ON INDIVIDUAL SEQUENCES
        for i, seqIdx in enumerate(slctdSeqs):
    # random PCR with bias using brute force        
            for n in xrange(pcrCycleNum):
                slctdSeqs[seqIdx][0] += binom(slctdSeqs[seqIdx][0], (pcrYld + slctdSeqs[seqIdx][2])) #amplify count
                totalseqs += slctdSeqs[seqIdx][0]
        print totalseqs
        m = mutDist.rvs(size=totalseqs) #draw no. of mutations for each seq copy
        print("mutation draws taken")
        for i, seqIdx in enumerate(slctdSeqs):
            if i==0:
                muts = m[:slctdSeqs[seqIdx][0]]
                muts = muts[muts != 0]
            else:
                muts = m[slctdSeqs[oldseqIdx][0]:slctdSeqs[oldseqIdx][0]+slctdSeqs[seqIdx][0]]
                muts = muts[muts != 0]
            oldseqIdx = seqIdx
            for i, mutNum in enumerate(muts): # for each mutation instance
                randPos = np.random.randint(seqLength-1, size=mutNum) #pick random nt positions for mutation
                if seqIdx not in mutatedPool: #if seq never mutated before
                    slctdSeqs[seqIdx][0] -= 1 #decrement seq count of wild type
                    mutatedPool.setdefault(seqIdx, []).append(1) #add to mutated pool
                    mutatedPool.setdefault(seqIdx, []).append(randPos) #append mutation positions
                else: #if seq previously mutated
                    slctdSeqs[seqIdx][0] -= 1 #decrement no. of wild type seq
                    mutatedPool[seqIdx][0] += 1 #increment no. of seq to be mutated
                    mutatedPool.setdefault(seqIdx, []).append(randPos) #append mutation positions 
        amplfdSeqs = slctdSeqs
        print("sequence amplification has been carried out")
        return amplfdSeqs, mutatedPool


#pcr modelled as ideal process that doubles each sequence per cycle
#assume no bias and 100% efficiency
    def scaleTest(self, scale, dimension):
        for s in xrange(slctd_seqs):
            slctd_seqs[seq][0]*=(2**pcr_cycles) 
        amplfd_seqs = slctd_seqs
        print("sequence amplification has been carried out")
        return amplfd_seqs


#pcr modelled as ideal process that doubles each sequence per cycle
#assume no bias and 100% efficiency
    def ampIdeal(self, slctd_seqs, pcr_cycles):
        for seq in slctd_seqs:
            slctd_seqs[seq][0]*=(2**pcr_cycles) 
        amplfd_seqs = slctd_seqs
        print("sequence amplification has been carried out")
        return amplfd_seqs

##pcr modelled as bernoulli test per cycle
    def ampEfficiencyBrute(self, slctdSeqs, pcrCycleNum, pcrYld):
        for n in xrange(pcrCycleNum):
            for seq in slctdSeqs:
                slctdSeqs[seq][1] += binom(slctdSeqs[seq][0], pcrYld)
        amplfdSeqs = slctdSeqs
        print("sequence amplification has been carried out")
        return amplfdSeqs


    def gaussianTimed(self, scale):
        start_time = time.time()
        for seq in xrange(scale):
            stats.norm.rvs(5, np.sqrt(4))
        print("Random Gaussian took %s seconds" % (time.time() - start_time))
       
        return 0

##This method is designed to measure the time performance of random PCR
# using the Brute-force approach
#pcr modelled as bernoulli test per cycle
    def bruteRandomPCRtimed(self, scale, initialCount, pcrCycleNum, pcrYield):
        start_time = time.time()
        print("scale is "+str(scale))
        for seq in xrange(scale):
            currentCount = initialCount
            for n in xrange(pcrCycleNum):
                currentCount += binom(initialCount, pcrYield)
        print("sequence amplification has been carried out")
        print("brute random PCR took %s seconds" % (time.time() - start_time))
        return 0


#This method is designed to test the time performance of random PCR using 
# the Gaussian Mixture Emitting model.
    def machineRandomPCRtimed(self, machine, scale, initialCount, pcrCycleNum, pcrYield, maxGauss):
        start_time = time.time()
        conditions = np.array([initialCount, pcrCycleNum, pcrYield])
        conditions = conditions.reshape(1, -1)
        gmmParams = np.zeros(12)
        gmm = GMM(n_components=maxGauss)
        initData = np.zeros((maxGauss, 1))
        gmm.fit(initData)
        print("scale is "+str(scale))
        for i in xrange(maxGauss):
            gmm.means_[i][0] = machine[i*3].predict(conditions)[0]
            gmm.covars_[i][0] = machine[(i*3)+1].predict(conditions)[0]
            gmm.weights_[i] = machine[(i*3)+2].predict(conditions)[0]
        currentCount = 0
        currentCount = gmm.sample(scale)
        print("sequence amplification has been carried out")
        print("machine random PCR took %s seconds" % (time.time() - start_time))
        return currentCount
      
      
    def dataWriteTest(self, scale, fileName):
        start_time = time.time()
        with open(fileName, 'w') as f:
            for i in range(scale):
                f.write(str(i)+"\n")
        print("Writing data to file took %s seconds" % (time.time() - start_time))
        return 0
       

##  TEST AREA
#amp = Amplification()
#amp.machineRandomPCRtimed(svrModel, 100, 1, 25, 0.5, 4)

# This method is intended to generate sequence population data after PCR
# This uses a brute-force method that models each individual cycle as Bernoulli test
# A histogram is generated from the data to represent a probability density function
# The pdf is then fitted with a Gaussian Mixture Model 
# Training carried out using expectation-maximization algorithm
############################################################################
# PARAMS:                                                                  #
# initialCount = the initial frequence of each sequence prior to PCR       #
# pcrCycleNum = the number of cycles to be performed prior to pdf generation #
# pcrYield = the average yield of PCR                                      #
# dataPoints = the initial number of unique sequences prior to PCR         #
# gaussNum = the max number of gaussians used to fit the pdf               #
# OUTPUT:                                                                  #
# bestModel = the best performing GMM model that approximated the pdf      #
# A plot of the hist of the actual population (green)  overlayed onto the  #
# GMM model (blue) and showing the individual gaussian components (red)    #
############################################################################
    def BruteGMM(self, initialCount, pcrCycleNum, pcrYield, dataPoints, gaussNum):
        amplfdSeqs = np.zeros((dataPoints, 1)) #create vector for seqs
        amplfdSeqsHist = np.zeros(dataPoints)
        for i in range(dataPoints):
            amplfdSeqs[i] = initialCount #assign init count to each seq
            amplfdSeqsHist[i] = initialCount
            for n in range(pcrCycleNum): #model each cycle as Bernoulli test
                amplfdSeqs[i] += binom(amplfdSeqs[i], pcrYield) 

                amplfdSeqsHist[i] += binom(amplfdSeqsHist[i], pcrYield) 
      
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
        expctdSeqNum = initialCount*(1+pcrYield)**(pcrCycleNum)
        varSeqNum = -1*initialCount*(1-pcrYield)*(1+pcrYield)**(pcrCycleNum-1)*(1-(1+pcrYield)**(pcrCycleNum+1))       
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
        ax.text(0.95, 0.95, r"Initial count = "+str(initialCount)+'\n'+"No. of cycles = "+str(pcrCycleNum)+'\n'+"Yield = "+str(pcrYield), verticalalignment='top', horizontalalignment='right', transform=ax.transAxes, color='black', fontsize=10)
        ax.text(0.935, 0.65, r"$ \begin{pmatrix} %s  \end{pmatrix}$" % annot, verticalalignment='top', horizontalalignment='right', transform=ax.transAxes, color='black', fontsize=10)


        # save plot
        plt.grid()
        plt.savefig('pcrDistEst_n'+str(pcrCycleNum)+'_i'+str(initialCount)+'_y'+str(pcrYield)+'.pdf', format='pdf')
        #plt.close()
        #plt.show()
        return bestModel

# TEST AREA - TO BE DELETED
#amp = Amplification()
#amp.BruteGMM(3, 5, 0.85, 10000, 20)

    def BruteGMMlogged(self, initialCount, pcrCycleNum, pcrYield, dataPoints, gaussNum):
        amplfdSeqs = np.zeros((dataPoints, 1)) #create vector for seqs
        amplfdSeqsHist = np.zeros(dataPoints)

    
        #calculate exact solutions to the expectation and variance of the dist
        expctdSeqNum = initialCount*(1+pcrYield)**(pcrCycleNum)

        varSeqNum = -1*initialCount*(1-pcrYield)*(1+pcrYield)**(pcrCycleNum-1)*(1-(1+pcrYield)**(pcrCycleNum+1))       
 
        for i in range(dataPoints):
            amplfdSeqs[i] = initialCount #assign init count to each seq
            amplfdSeqsHist[i] = initialCount
            for n in range(pcrCycleNum): #model each cycle as Bernoulli test
                amplfdSeqs[i] += binom(amplfdSeqs[i], pcrYield) 

                amplfdSeqsHist[i] += binom(amplfdSeqsHist[i], pcrYield) 
      
        for i in range(dataPoints):
            amplfdSeqs[i] = np.log(amplfdSeqs[i])

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
        #declare sample space 
        space = np.linspace(1, (expctdSeqNum + (np.sqrt(varSeqNum)*5)), dataPoints)
        loggedSpace = np.log(space)
        space = space.reshape(-1, 1)
        loggedSpace = loggedSpace.reshape(-1, 1)
        #calculate log likelihood of sample space using trained model
        logProbs, resps = gmmModel.score_samples(loggedSpace)
        loggedSpace = np.log(np.linspace(1, (expctdSeqNum + (np.sqrt(varSeqNum)*5)), dataPoints))
        #calculate prob density func of the model
        pdf = np.exp(logProbs)/np.exp(loggedSpace)
        #calculate prob density func of each component gaussian
        individualPDFs = resps * pdf[:, np.newaxis]
     
        fig = plt.figure()
        ax = fig.add_subplot(111)


        ax.hist(amplfdSeqsHist, bins=50, normed=True, facecolor='green', alpha=0.75, label='Brute-force')
        ax.plot(space, pdf, '-k', color='b', label='GMM')
        ax.plot(space, individualPDFs, '--k', color='r', label='component')
        ax.set_xlabel('Sequence Count')
        ax.set_ylabel('p(x)')
        ax.set_title('GMM Best-fit to Population Distribution')
        # create annotation
        annot = " \mu & \sigma^{2} & \omega \\"+"\\"
        for i, mu in enumerate(gmmModel.means_):
            annot += str(np.round(gmmModel.means_[i][0], 2))+" & "+str(np.round(gmmModel.covars_[i][0], 2))+" & "+str(np.round(gmmModel.weights_[i], 2))+" \\"+"\\ "
      
        #add plot annotations
        ax.text(0.95, 0.95, r"$Y_{0} = $"+str(initialCount)+'\n'+"$N = $"+str(pcrCycleNum)+'\n'+"$ \lambda = $"+str(pcrYield), verticalalignment='top', horizontalalignment='right', transform=ax.transAxes, color='black', fontsize=10)
        ax.text(0.935, 0.65, r"$ \begin{pmatrix} %s  \end{pmatrix}$" % annot, verticalalignment='top', horizontalalignment='right', transform=ax.transAxes, color='black', fontsize=10)


        # save plot
        plt.grid()
        plt.legend(loc=2, prop={'size': 6})
        plt.savefig('logGMM_pcrDistEst_n'+str(pcrCycleNum)+'_i'+str(initialCount)+'_y'+str(pcrYield)+'.pdf', format='pdf')
        #release memory dedicated for this plot
        plt.clf()
        #plt.show()
        return modelParams

# TEST AREA - TO BE DELETED
#amp = Amplification()
#params, loggedSpace, pdf = amp.BruteGMMlogged(3, 5, 0.85, 10000, 4)







    def BruteGMMnormed(self, initialCount, pcrCycleNum, pcrYield, dataPoints, gaussNum):
        amplfdSeqs = np.zeros((dataPoints, 1)) #create vector for seqs
        amplfdSeqsHist = np.zeros(dataPoints)

    
        #calculate exact solutions to the expectation and variance of the dist
        expctdSeqNum = initialCount*(1+pcrYield)**(pcrCycleNum)

        varSeqNum = -1*initialCount*(1-pcrYield)*(1+pcrYield)**(pcrCycleNum-1)*(1-(1+pcrYield)**(pcrCycleNum+1))       
 
        for i in range(dataPoints):
            amplfdSeqs[i] = initialCount #assign init count to each seq
            amplfdSeqsHist[i] = initialCount
            for n in range(pcrCycleNum): #model each cycle as Bernoulli test
                amplfdSeqs[i] += binom(amplfdSeqs[i], pcrYield) 
                amplfdSeqsHist[i] += binom(amplfdSeqsHist[i], pcrYield) 
      
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
        #ax.text(0.95, 0.95, r"Initial count = "+str(initialCount)+'\n'+"No. of cycles = "+str(pcrCycleNum)+'\n'+"Yield = "+str(pcrYield), verticalalignment='top', horizontalalignment='right', transform=ax.transAxes, color='black', fontsize=10)
        #ax.text(0.935, 0.65, r"$ \begin{pmatrix} %s  \end{pmatrix}$" % annot, verticalalignment='top', horizontalalignment='right', transform=ax.transAxes, color='black', fontsize=10)


        # save plot
        #plt.grid()
        #plt.savefig('pcrDistEst_n'+str(pcrCycleNum)+'_i'+str(initialCount)+'_y'+str(pcrYield)+'.pdf', format='pdf')
        #plt.close()
        #plt.show()
        return modelParams

# TEST AREA - TO BE DELETED
#amp = Amplification()
#params = amp.BruteGMMnormed(3, 5, 0.85, 10000, 8)



    def GMMTest(self, initialCount, pcrCycleNum, pcrYield, dataPoints, gaussNum):
        amplfdSeqs = np.zeros((dataPoints, 1)) #create vector for seqs
        amplfdSeqsHist = np.zeros(dataPoints)
        for i in range(dataPoints):
            amplfdSeqs[i] = initialCount #assign init count to each seq
            amplfdSeqsHist[i] = initialCount
            for n in range(pcrCycleNum): #model each cycle as Bernoulli test
                amplfdSeqs[i] += binom(amplfdSeqs[i], pcrYield) 

                amplfdSeqsHist[i] += binom(amplfdSeqsHist[i], pcrYield) 
      
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
        expctdSeqNum = initialCount*(1+pcrYield)**(pcrCycleNum) #expectation of amplification
        varSeqNum = -1*initialCount*(1-pcrYield)*(1+pcrYield)**(pcrCycleNum-1)*(1-(1+pcrYield)**(pcrCycleNum+1))
        #varSeqNum = initialCount*(1-pcrYield)(1+pcrYield)**(pcrCycleNum-1)(1-(1+pcrYield)**(pcrCycleNum+1))
        #varSeqNum = (1-pcrYield)*(1+pcrYield)**(pcrCycleNum-1)*((1+pcrYield)**(pcrCycleNum)-1) #variance of amplification
     
     
        #declare sample space 
        space = np.linspace(1, (expctdSeqNum + (np.sqrt(varSeqNum)*10)), dataPoints/100)
        #generate theoretical Gaussian pdf
        approxPDF = norm.pdf(space, expctdSeqNum, np.sqrt(varSeqNum))
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


        ax.hist(amplfdSeqsHist, bins=50, normed=True, facecolor='green', alpha=0.75, label='Brute-force')
        ax.plot(space, approxPDF, '-k', color='y', label='Theoretical Gaussian')
        ax.plot(space, pdf, '-k', color='b',label='Gaussian Mixture')
        ax.plot(space, individualPDFs, '--k', color='r', label='Gaussian components')
        ax.set_xlabel('Sequence Count')
        ax.set_ylabel('p(x)')
        ax.set_title('GMM Best-fit to Population Distribution' )
        # create annotation
        annot = " \mu & \sigma^{2} & \omega \\"+"\\"
        for i, mu in enumerate(gmmModel.means_):
            annot += str(np.round(gmmModel.means_[i][0], 2))+" & "+str(np.round(gmmModel.covars_[i][0], 2))+" & "+str(np.round(gmmModel.weights_[i], 2))+" \\"+"\\ "
      
        #add plot annotations
        ax.text(0.95, 0.95, r"$E[Y_{0}] = $"+str(initialCount)+'\n'+"$N = $"+str(pcrCycleNum)+'\n'+"$\lambda = $"+str(pcrYield), verticalalignment='top', horizontalalignment='right', transform=ax.transAxes, color='black', fontsize=10)
        ax.text(0.935, 0.65, r"$ \begin{pmatrix} %s  \end{pmatrix}$" % annot, verticalalignment='top', horizontalalignment='right', transform=ax.transAxes, color='black', fontsize=10)


        # save plot
        plt.grid()
        plt.legend(loc=2, prop={'size':6})
        plt.savefig('pcrDistEst_n'+str(pcrCycleNum)+'_i'+str(initialCount)+'_y'+str(pcrYield)+'.pdf', format='pdf')
      #plt.close()
      #plt.show()
        return gmmModel

# TEST AREA - TO BE DELETED
#amp = Amplification()
#amp.GMMTest(1, 15, 0.85, 10000, 4)






    def bruteHist(self, initialCount, pcrCycleNum, pcrYield, dataPoints, binsNum):
        amplfdSeqs = np.zeros((dataPoints, 1)) #create vector for seqs
        amplfdSeqsHist = np.zeros(dataPoints)
        for i in range(dataPoints):
            amplfdSeqs[i] = initialCount #assign init count to each seq
            amplfdSeqsHist[i] = initialCount
            for n in range(pcrCycleNum): #model each cycle as Bernoulli test
                amplfdSeqs[i] += binom(amplfdSeqs[i], pcrYield) 

                amplfdSeqsHist[i] += binom(amplfdSeqsHist[i], pcrYield) 
      
        n, bins, patches = plt.hist(amplfdSeqsHist, binsNum, normed=True, facecolor='green', alpha=0.75)
        #plt.plot(space, pdf, '-k', color='b')
        #plt.plot(space, individualPDFs, '--k')
        plt.xlabel('Sequence Count')
        plt.ylabel('p(x)')
        plt.title('GMM Best-fit')
        plt.savefig('pcrDistEst_n'+str(pcrCycleNum)+'_i'+str(initialCount)+'_y'+str(pcrYield)+'.pdf', format='pdf')
        plt.show()
# TEST AREA - TO BE DELETED
#amp = Amplification()
#amp.BruteHist(2, 20, 0.65, 10000, 20)




    def ampEffBruteTest(self, initialCount, pcrCycleNum, pcrYield, dataPoints):
        amplfdSeqs = np.zeros(dataPoints)
        for i in range(dataPoints):
            amplfdSeqs[i] = initialCount
            for n in range(pcrCycleNum):
                amplfdSeqs[i] += binom(amplfdSeqs[i], pcrYield)
        n, bins, patches = plt.hist(amplfdSeqs, 10000, normed=1, facecolor='green', alpha=0.75)
        plt.xlabel('Sequence Count')
        plt.ylabel('Probability')
        plt.title('A histogram plot of '+str(dataPoints)+' sequences after '+str(pcrCycleNum)+' PCR ciycles with '+str(pcrYield)+'% yield.\n'+'Initial count = '+str(initialCount))
        plt.savefig('pcrDistEst_n'+str(pcrCycleNum)+'_i'+str(initialCount)+'_y'+str(pcrYield)+'.png')
        plt.close()
##pcr modelled as Bernoulli process after n iterations
# assume no bias or mutations
# compute expectation after n cycles for each selected sequence
# estimate amplified number as expectation of pgf 
    def ampEfficiencyDefinite(self, slctdSeqs, pcrCycleNum, pcrYld):
        for seq in slctdSeqs:
            slctdSeqs[seq][0] = slctdSeqs[seq][0]*(1+pcrYld)**(pcrCycleNum)
        amplfdSeqs = slctdSeqs
        print("sequence amplification has been carried out")
        return amplfdSeqs



##pcr modelled as Bernoulli process after n iterations
# assume no bias or mutations
# compute expectation and variance after n cycles for each selected sequence
# estimate amplified number as draw from Gaussian distribution 
# Law of Large Numbers
# Binomial distribution approximated as Gaussian due to large copy numbers
    def ampEfficiencyStochastic(self, slctdSeqs, pcrCycleNum, pcrYld):
        for seq in slctdSeqs:
            meanSeqNum = (slctdSeqs[seq][0]*(1+pcrYld))**(pcrCycleNum) #expectation
            varSeqNum = -1*slctdSeqs[seq][0]*(1-pcrYld)*(1+pcrYld)**(pcrCycleNum-1)*(1-(1+pcrYld)**(pcrCycleNum+1))       
            slctdSeqs[seq][0] = int(stats.norm.rvs(meanSeqNum, np.sqrt(varSeqNum))) #draw number from Gaussian dist  
        amplfdSeqs = slctdSeqs
        print("sequence amplification has been carried out")
        return amplfdSeqs



# Assume all mutations happened at end of PCR (i.e. mutants are not amplified till next round)
## ADD seqLength as arg to function
## EDITED variance calculation for amplification
    def ampEffMutStochastic(self, slctdSeqs, seqLength, pcrCycleNum, pcrYld, errorRate):
        mutatedPool = {} #initialize pool of seqs to be mutated
        gamma = np.arange(seqLength) #initialize vector of gamma values 
        mut_m = np.arange(seqLength) #initialize vector of mut no. (0-seqLen)
        prob_m = np.arange(seqLength) #initialize corresponding vector of probs
# Mutation Distribution
        for i in mut_m: # for each possible no. of muts, compute corresponding probability
            for j in range(pcrCycleNum): # for each pcr cycle
                gamma[i] += np.exp(-j*errorRate*seqLength)*math.factorial(pcrCycleNum)/(math.factorial(j)*math.factorial(pcrCycleNum - j))*(j)**(i)*pcrYld**(j) #calculate sum (SEE PAPER)
            prob_m[i] = (errorRate*seqLength)**(i)*gamma[i]/(math.factorial(i)*(1+pcrYld)**(pcrCycleNum)) #compute prob of m[i] muts
        mutDist = stats.rv_discrete(name='mutDist', values=(mut_m, prob_m)) #compute discrete mutation distribution
        print("Discrete Mutation Distribution has been computed")
# PCR Amplification
        for seq in slctdSeqs:
            meanSeqNum = slctdSeqs[seq][0]*(1+pcrYld)**(pcrCycleNum) #expectation of amplification
            varSeqNum = -1*slctdSeqs[seq][0]*(1-pcrYld)*(1+pcrYld)**(pcrCycleNum-1)*(1-(1+pcrYld)**(pcrCycleNum+1))       
            slctdSeqs[seq][0] = int(stats.norm.rvs(meanSeqNum, np.sqrt(varSeqNum))) #draw int from Gaussian dist
            m = mutDist.rvs(size=slctdSeqs[seq][0]) #draw no. of mutations for each seq copy
            muts = m[m != 0] #keep only copies to be mutated (i.e. m >= 1)
            for mutNum, i in enumerate(muts): # for each mutation instance
                randPos = np.random.randint(seqLength, size=mutNum) #pick random nt positions for mutation  ##add seqLen as argumant to function
            if seq not in mutatedPool: #if seq never mutated before
                slctdSeqs[seq][0] -= 1 #decrement seq count of wild type
                mutatedPool.setdefault(seq, []).append(1) #add to mutated pool
                mutatedPool.setdefault(seq, []).append(randPos) #append mutation positions
            else: #if seq previously mutated
                slctdSeqs[seq][0] -= 1 #decrement no. of wild type seq
                mutatedPool[seq][0]+=1 #increment no. of seq to be mutated
                mutatedPool.setdefault(seq, []).append(randPos) #append mutation positions
            
        amplfdSeqs = slctdSeqs
        print("sequence amplification has been carried out")
        return amplfdSeqs, mutatedPool


# Assume all mutations happened at end of PCR (i.e. mutants are not amplified till next round)
## ADD seqLength as arg to function
## EDITED variance calculation for amplification
    def ampEffMutDefinite(self, slctdSeqs, seqLength, pcrCycleNum, pcrYld, errorRate):
        mutatedPool = {} #initialize pool of seqs to be mutated
        gamma = np.arange(seqLength) #initialize vector of gamma values 
        mut_m = np.arange(seqLength) #initialize vector of mut no.s (0-seqLen)
        prob_m = np.arange(seqLength) #initialize corresponding vector of probs
      # Generate distribution for generation number        
        prob_n = np.zeros(pcrCycleNum+1)
        genNum = np.arange(pcrCycleNum+1)
        for n in range(pcrCycleNum+1): #for each generation
            prob_n[n] = math.factorial(pcrCycleNum)*(pcrYld)**(n)/(math.factorial(n)*math.factorial(pcrCycleNum - n)*(1+pcrYld)**(pcrCycleNum)) #compute prob of picking seq from gen n
        genDist = stats.rv_discrete(name='genDist', values=(genNum, prob_n)) #construct gen Distribution
        print("Discrete Generation Number Distribution has been computed")
# Mutation Distribution
        for i in mut_m: # for each possible no. of muts, compute corresponding probability
            for j in range(pcrCycleNum+1): # for each pcr generation
                gamma[i] += np.exp(-j*errorRate*seqLength)*math.factorial(pcrCycleNum)/(math.factorial(j)*math.factorial(pcrCycleNum - j))*(j)**(i)*pcrYld**(j) #calculate sum (SEE PAPER)
            prob_m[i] = (errorRate*seqLength)**(i)*gamma[i]/(math.factorial(i)*(1+pcrYld)**(pcrCycleNum)) #compute prob of m[i] muts
        mutDist = stats.rv_discrete(name='mutDist', values=(mut_m, prob_m)) #construct discrete mutation distribution
        print("Discrete Mutation Distribution has been computed")
# PCR Amplification
        for seq in slctdSeqs: #for each seq in selected pool
            slctdSeqs[seq][0] = int(slctdSeqs[seq][0]*(1+pcrYld)**(pcrCycleNum)) #expectation of amplification
            m = mutDist.rvs(size=slctdSeqs[seq][0]) #draw no. of mutations for each seq copy
            muts = m[m != 0] #keep only copies to be mutated (i.e. m >= 1)
            g = genDist.rvs(muts.shape) #draw gen no. for each mut'd seq copy
            for i, mutNum in enumerate(muts): # for each mutation instance
                randPos = np.random.randint(seqLength, size=mutNum) #pick random nt positions for mutation  ##add seqLen as argumant to function
            if seq not in mutatedPool: #if seq never mutated before
                slctdSeqs[seq][0] -= (1+pcrYld)**(pcrCycleNum-g[i]) #decrement seq count of wild type
                mutatedPool.setdefault(seq, []).append(np.zeros(muts.shape)) # add number of copies to be mutated
                mutatedPool[seq][0][i] = ((1+pcrYld)**(pcrCycleNum-g[i])) #add expected amplfd count 
                mutatedPool.setdefault(seq, []).append(randPos) #append mutation positions
            else: #if seq previously mutated
                slctdSeqs[seq][0] -= (1+pcrYld)**(pcrCycleNum-g[i]) #decrement no. of wild type seq
                mutatedPool[seq][0][i] = (1+pcrYld)**(pcrCycleNum-g[i]) #increment no. of seq to be mutated
                mutatedPool.setdefault(seq, []).append(randPos) #append mutation positions
            
        amplfdSeqs = slctdSeqs
        print("sequence amplification has been carried out")
        return amplfdSeqs, mutatedPool


    def mathApproxHist(self, initialCount, pcrCycleNum, pcrYield, dataPoints, binsNum):
        amplfdSeqs = np.zeros((dataPoints, 1)) #create vector for seqs
        amplfdSeqsHist = np.zeros(dataPoints)
        for i in range(dataPoints):
            amplfdSeqs[i] = initialCount #assign init count to each seq
            amplfdSeqsHist[i] = initialCount
            for n in range(pcrCycleNum): #model each cycle as Bernoulli test
                amplfdSeqs[i] += np.random.binomial(amplfdSeqs[i], pcrYield) 

                amplfdSeqsHist[i] += np.random.binomial(amplfdSeqsHist[i], pcrYield) 
        expctdSeqNum = initialCount*(1+pcrYield)**(pcrCycleNum) #expectation of amplification
        varSeqNum = -1*initialCount*(1-pcrYield)*(1+pcrYield)**(pcrCycleNum-1)*(1-(1+pcrYield)**(pcrCycleNum+1))
    #declare sample space 
        space = np.linspace(1, (expctdSeqNum + (np.sqrt(varSeqNum)*10)), dataPoints/100)
        approxPDF = norm.pdf(space, expctdSeqNum, np.sqrt(varSeqNum))

        fig = plt.figure()
        ax = fig.add_subplot(111)


        ax.hist(amplfdSeqsHist, bins=50, normed=True, facecolor='green', alpha=0.75, label='Brute-force')
        ax.plot(space, approxPDF, '-k', color='y', label='Theoretical Gaussian')
        ax.set_xlabel('Sequence Count')
        ax.set_ylabel('p(x)')
        ax.set_title('Theoretical Gaussian Approximation') 
        ax.text(0.95, 0.95, r"$E[Y_{0}] = $"+str(initialCount)+'\n'+"$N = $"+str(pcrCycleNum)+'\n'+"$\lambda = $"+str(pcrYield), verticalalignment='top', horizontalalignment='right', transform=ax.transAxes, color='black', fontsize=10)
        plt.legend(loc=2, prop={'size':6})
        plt.grid()
        plt.savefig('gauss_pcrDistEst_n'+str(pcrCycleNum)+'_i'+str(initialCount)+'_y'+str(pcrYield)+'.pdf', format='pdf')
        return space

##TEST AREA
#amp = Amplification()
#space = amp.mathApproxHist(1, 15, 0.85, 10000, 50)

    def BruteVBGMM(self, initialCount, pcrCycleNum, pcrYield, dataPoints, gaussNum, a, r):
        amplfdSeqs = np.zeros((dataPoints, 1)) #create vector for seqs
        amplfdSeqsHist = np.zeros(dataPoints)
        for i in range(dataPoints):
            amplfdSeqs[i] = initialCount #assign init count to each seq
            amplfdSeqsHist[i] = initialCount
            for n in range(pcrCycleNum): #model each cycle as Bernoulli test
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
      

        expctdSeqNum = initialCount*(1+pcrYield)**(pcrCycleNum)
        varSeqNum = pcrYield*(1-pcrYield)*(initialCount*(1+pcrYield)**(pcrCycleNum - 1))*((initialCount*(1+pcrYield)**(pcrCycleNum))-1)/(initialCount*(1+pcrYield)-1)


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
    def BruteDPGMM(self, initialCount, pcrCycleNum, pcrYield, dataPoints, gaussNum):
        amplfdSeqs = np.zeros((dataPoints, 1)) #create vector for seqs
        for i in range(dataPoints):
            amplfdSeqs[i] = initialCount #assign init count to each seq
            for n in range(pcrCycleNum): #model each cycle as Bernoulli test
                amplfdSeqs[i] += np.random.binomial(amplfdSeqs[i], pcrYield) 
      
        amplfdSeqs = amplfdSeqs.reshape(-1, 1)

      
        dpgmmModels = DPGMM(n_components=gaussNum, covariance_type='full').fit(amplfdSeqs)
      

        dpgmmAIC = dpgmmModels.aic(amplfdSeqs)
        dpgmmBIC = dpgmmModels.bic(amplfdSeqs)


        bestModel = dpgmmModels
      

        expctdSeqNum = initialCount*(1+pcrYield)**(pcrCycleNum)
        varSeqNum = pcrYield*(1-pcrYield)*(initialCount*(1+pcrYield)**(pcrCycleNum - 1))*((initialCount*(1+pcrYield)**(pcrCycleNum))-1)/(initialCount*(1+pcrYield)-1)


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


