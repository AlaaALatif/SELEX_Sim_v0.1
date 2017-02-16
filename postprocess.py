import sys
import numpy as np
import matplotlib
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick


def dataAnalysis(seqLength, roundNum, outputFileNames, plots, distanceMeasure):
    avgDist_per_rnd = np.zeros(roundNum)
    weighted_avgDist_per_rnd = np.zeros(roundNum)
    total_seqs_freqs = np.zeros(roundNum)
    uniq_seqs_freqs = np.zeros(roundNum)
    distFreqs = np.zeros((roundNum, seqLength+1))
    weighted_distFreqs = np.zeros((roundNum, seqLength+1))
    for rnd in xrange(roundNum):
        total_seq_num = 0
        uniq_seq_num = 0
        distance = 0
        weighted_distance = 0
        with open(outputFileNames+"_R"+str(rnd+1)) as SELEX_round:
            for line in SELEX_round:
                columns = line.split()
                distance += int(columns[2])
                weighted_distance += int(columns[1])*int(columns[2])
                total_seq_num += int(columns[1])
                uniq_seq_num += 1
                distFreqs[rnd][int(columns[2])] += 1
                weighted_distFreqs[rnd][int(columns[2])] += int(columns[1])
        total_seqs_freqs[rnd] = total_seq_num
        avgDist_per_rnd[rnd] = int(distance/uniq_seq_num)
        weighted_avgDist_per_rnd[rnd] = int(weighted_distance/total_seq_num)
        uniq_seqs_freqs[rnd] = uniq_seq_num
	for i in xrange(seqLength+1):
		distFreqs[rnd][i] /= uniq_seqs_freqs[rnd]
		weighted_distFreqs[rnd][i] /= total_seqs_freqs[rnd]
    print(weighted_avgDist_per_rnd)
    print(avgDist_per_rnd)
    with open(outputFileNames+"_processed_results", 'w') as p:
        for rnd in xrange(roundNum): 
            p.write(str(int(total_seqs_freqs[rnd]))+'\t')
            p.write(str(int(uniq_seqs_freqs[rnd]))+'\t')
            p.write(str(int(avgDist_per_rnd[rnd]))+'\t')
            p.write(str(int(weighted_avgDist_per_rnd[rnd]))+'\t')
            for l in xrange(seqLength+1):
                p.write(str(int(distFreqs[rnd][l]))+'\t')
                p.write(str(int(weighted_distFreqs[rnd][l]))+'\t')
            p.write('\n')
    if(plots==True):
        if(distanceMeasure=="hamming"):
            roundNumAxis = np.linspace(1, roundNum, roundNum)
            fig0, axes = plt.subplots(2, 2)
            plotsList = [total_seqs_freqs, uniq_seqs_freqs, weighted_avgDist_per_rnd, avgDist_per_rnd]
            for i, ax in enumerate(axes.reshape(-1)):
                ax.plot(roundNumAxis, plotsList[i])
            fig0.text(0.5, 0.04, 'Round Number', ha='center')
            fig0.text(0.04, 0.5, 'Average Distance', va='center', rotation='vertical')
            fig0.savefig(str(outputFileNames)+"_SELEX_Analytics_distance.pdf")
            fig1, axes = plt.subplots(2, 3)
            for i, ax in enumerate(axes.reshape(-1)):
                for d in range(3):
                    ax.plot(roundNumAxis, distFreqs[:,d+(3*i)+1], label='d = '+str(d+(3*i)+1))
                ax.ticklabel_format(syle='sci', axis='y', scilimits=(0,0))
                ax.legend(prop={'size':6})
            fig1.text(0.5, 0.04, 'Round Number', ha='center')
            fig1.text(0.04, 0.5, 'Fractional Frequency', va='center', rotation='vertical')
            fig1.savefig(str(outputFileNames)+"_SELEX_Analytics_distFreqs.pdf")
            # weighted fractional sequency plots
            fig2, axes = plt.subplots(2, 3)
            for i, ax in enumerate(axes.reshape(-1)):
                for d in range(3):
                    ax.plot(roundNumAxis, weighted_distFreqs[:,d+(3*i)+1], label='d = '+str(d+(3*i)+1))
                ax.ticklabel_format(syle='sci', axis='y', scilimits=(0,0))
                ax.legend(prop={'size':6})
            fig2.text(0.5, 0.04, 'Round Number', ha='center')
            fig2.text(0.04, 0.5, 'Fractional Frequency', va='center', rotation='vertical')
            fig2.savefig(str(outputFileNames)+"_SELEX_Analytics_weighted_distFreqs.pdf")
        elif(distanceMeasure=="basepair"):
            roundNumAxis = np.linspace(1, roundNum, roundNum)
            plotsList = [total_seqs_freqs, uniq_seqs_freqs, weighted_avgDist_per_rnd, avgDist_per_rnd]
            fig0, axes = plt.subplots(2, 2)
            for i, ax in enumerate(axes.reshape(-1)):
                ax.plot(roundNumAxis, plotsList[i])
            fig0.text(0.5, 0.04, 'Round Number', ha='center')
            fig0.text(0.04, 0.5, 'Average Distance', va='center', rotation='vertical')
            fig0.savefig(str(outputFileNames)+"_SELEX_Analytics_distance.pdf")
            #figures for distance analytics
            fig1, axes = plt.subplots(2, 3)
            for i, ax in enumerate(axes.reshape(-1)):
                for d in range(3):
                    ax.plot(roundNumAxis, distFreqs[:,d+(3*i)+1], label='d = '+str(d+(3*i)+1))
                ax.ticklabel_format(syle='sci', axis='y', scilimits=(0,0))
                ax.legend(prop={'size':6})
            fig1.text(0.5, 0.04, 'Round Number', ha='center')
            fig1.text(0.04, 0.5, 'Fractional Frequency', va='center', rotation='vertical')
            fig1.savefig(str(outputFileNames)+"_SELEX_Analytics_distFreqs.pdf")
            # weighted fractional sequency plots
            fig2, axes = plt.subplots(2, 3)
            for i, ax in enumerate(axes.reshape(-1)):
                for d in range(3):
                    ax.plot(roundNumAxis, weighted_distFreqs[:,d+(3*i)+1], label='d = '+str(d+(3*i)+1))
                ax.ticklabel_format(syle='sci', axis='y', scilimits=(0,0))
                ax.legend(prop={'size':6})
            fig2.text(0.5, 0.04, 'Round Number', ha='center')
            fig2.text(0.04, 0.5, 'Fractional Frequency', va='center', rotation='vertical')
            fig2.savefig(str(outputFileNames)+"_SELEX_Analytics_weighted_distFreqs.pdf")
        else:
            return
#TEST
#dataAnalysis(20, 15, "complex", True, "basepair")
