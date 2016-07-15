import sys
import time
import random
import math
from scipy.stats import norm
import numpy as np
import Predictor
import matplotlib
matplotlib.use('Agg')
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt

## USER INPUT PARAMETERS
firstDatafile = str(sys.argv[1])

secondDatafile = str(sys.argv[2])

## MATRIX INITIALIZATION
firstCon = np.zeros((30, 6))
secondCon = np.zeros((30, 6))
thirdCon = np.zeros((30, 6))
fourthCon = np.zeros((30, 6))
## DATA COLLECTION
firstData = open(firstDatafile, 'r')

for j in range(30):
    line = firstData.readline()
    firstCon[j][0] = float(line.split()[0]) #initial count
    firstCon[j][2] = float(line.split()[1]) #cycle no.
    firstCon[j][1] = float(line.split()[2]) #yield
    firstCon[j][3] = float(line.split()[3]) #mean of 1st component
    firstCon[j][4] = float(line.split()[4]) #variance of 1st component
    firstCon[j][5] = float(line.split()[5]) #weight of 1st component

for j in range(30):
    line = firstData.readline()

for j in range(30):
    line = firstData.readline()
    secondCon[j][0] = float(line.split()[0]) #initial count
    secondCon[j][2] = float(line.split()[1]) #cycle no.
    secondCon[j][1] = float(line.split()[2]) #yield
    secondCon[j][3] = float(line.split()[3]) #mean of 1st component
    secondCon[j][4] = float(line.split()[4]) #variance of 1st component
    secondCon[j][5] = float(line.split()[5]) #weight of 1st component

firstData.close()
## SecondDatafile contains all data pofloat at 40% yield
secondData = open(secondDatafile, 'r')

for j in range(30):
    line = secondData.readline()
    thirdCon[j][0] = float(line.split()[0]) #initial count
    thirdCon[j][2] = float(line.split()[1]) #cycle no
    thirdCon[j][1] = float(line.split()[2]) #yield
    thirdCon[j][3] = float(line.split()[3]) #mean of 1st comp
    thirdCon[j][4] = float(line.split()[4]) #variance of 1st comp
    thirdCon[j][5] = float(line.split()[5]) #weight of 1st comp
    # initial count = 2 (SKIP)
for j in range(30):
    line = secondData.readline()
    # initial count = 3
for j in range(30):
    line = secondData.readline()
    fourthCon[j][0] = float(line.split()[0]) #initial count
    fourthCon[j][2] = float(line.split()[1]) #cycle no.
    fourthCon[j][1] = float(line.split()[2]) #yield
    fourthCon[j][3] = float(line.split()[3]) #mean of 1st comp
    fourthCon[j][4] = float(line.split()[4]) #var of 1st comp
    fourthCon[j][5] = float(line.split()[5]) #weight of 1st comp
            
secondData.close()

## DATA VISUALIZATION

# Plot 1: effect of cycle no on normalized mean of 1st GMM component

space = np.linspace(1, 30, 30)

fig1 = plt.figure()
ax = fig1.add_subplot(111)

ax.plot(space, np.log(firstCon[:,3]), "o", color='r')
ax.plot(space, np.log(secondCon[:,3]), "o", color='b')
ax.plot(space, np.log(thirdCon[:,3]), "o", color='y')
ax.plot(space, np.log(fourthCon[:,3]), "o", color='g')

ax.set_xlabel('Cycle no.')
ax.set_ylabel('log normed mean')
ax.set_title('Effect of cycle number on log normed mean of 1st GMM component')

plt.grid()
plt.savefig('pcr_logmean_g1_i1_3_y4_8.pdf', format='pdf')



class Visualizations:

    def gmmValidation(self, xTest, yTest, preds, maxGauss):
        realPDF = np.zeros((xTest.shape[0], 100))
        predPDF = np.zeros((xTest.shape[0], 100))
        individual_realPDFs = np.zeros((xTest.shape[0], maxGauss, 100))
        individual_predPDFs = np.zeros((xTest.shape[0], maxGauss, 100))
        for i, sample in enumerate(xTest):
            initialCount = int(sample[0])
            cycleNum = int(sample[1])
            yld = sample[2]
    #calculate exact solutions to the expectation and variance of the dist
            expctdSeqNum = initialCount*(1+yld)**(cycleNum)
            varSeqNum = yld*(1-yld)*(initialCount*(1+yld)**(cycleNum - 1))*((initialCount*(1+yld)**(cycleNum))-1)/(initialCount*(1+yld)-1)
            space = np.linspace(1, (expctdSeqNum + (np.sqrt(varSeqNum)*10)), 100)
            normedSpace = np.linspace((1/(2*expctdSeqNum)), ((expctdSeqNum + (np.sqrt(varSeqNum)*10))/(2*expctdSeqNum)), 100)
            for j in xrange(0, (maxGauss*3), 3):
                individual_realPDFs[i][j/3] = yTest[i][j+2]*norm.pdf(normedSpace, yTest[i][j], np.sqrt(yTest[i][j+1]))
                realPDF[i] += individual_realPDFs[i][j/3]
                individual_predPDFs[i][j/3] = preds[i][j+2]*norm.pdf(normedSpace, preds[i][j], np.sqrt(preds[i][j+1]))
                predPDF[i] += individual_predPDFs[i][j/3]

            fig = plt.figure()
            ax = fig.add_subplot(111)


            #ax.hist(amplfdSeqsHist, bins=50, normed=True, facecolor='green', alpha=0.75)
            ax.plot(space, (realPDF[i]/(2*expctdSeqNum)), '-k', color='b', label='real GMM')
            ax.plot(space, (predPDF[i]/(2*expctdSeqNum)), '-k', color='g', label='predicted GMM')
            for k in xrange(individual_realPDFs.shape[1]):
                ax.plot(space, (individual_realPDFs[i][k]/(2*expctdSeqNum)), '--k', color='r')
                ax.plot(space, (individual_predPDFs[i][k]/(2*expctdSeqNum)), '--k', color='y')
            ax.set_xlabel('Sequence Count')
            ax.set_ylabel('p(x)')
            ax.set_title('GMM of PCR Population Distribution (real vs predicted model)')
      # create annotation
      #annot = " \mu & \sigma^{2} & \omega \\"+"\\"
      #for i, mu in enumerate(gmmModel.means_):
          #annot += str(np.round(gmmModel.means_[i][0], 2))+" & "+str(np.round(gmmModel.covars_[i][0], 2))+" & "+str(np.round(gmmModel.weights_[i], 2))+" \\"+"\\ "
      
      #add plot annotations
            ax.text(0.95, 0.95, r"$E[Y_{0}] = $"+str(initialCount)+'\n'+"$N = $"+str(cycleNum)+'\n'+"$ \lambda = $"+str(yld), verticalalignment='top', horizontalalignment='right', transform=ax.transAxes, color='black', fontsize=10)
      #ax.text(0.935, 0.65, r"$ \begin{pmatrix} %s  \end{pmatrix}$" % annot, verticalalignment='top', horizontalalignment='right', transform=ax.transAxes, color='black', fontsize=10)


      # save plot
            plt.grid()
            plt.legend(loc=2, prop={'size':6})
            plt.savefig('gmmPredictor_n'+str(cycleNum)+'_i'+str(initialCount)+'_y'+str(yld)+'.pdf', format='pdf')
      #plt.close()
      #plt.show()
 
        return realPDF, predPDF

viz = Visualizations()






