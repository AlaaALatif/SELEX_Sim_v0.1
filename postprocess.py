import sys
import numpy as np
import matplotlib
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick


def dataAnalysis(seqLength, roundNum, outputFileNames, plots):

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
        roundNumAxis = np.linspace(1, roundNum, roundNum)
        #figures for distance analytics
        cm = plt.get_cmap('gist_rainbow')    
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.set_color_cycle([cm(1.*i/seqLength+1) for i in range(seqLength+1)])
        for i in xrange(seqLength+1):
            ax.plot(roundNumAxis, distFreqs[:,i], label='d = '+str(i))
        ax.legend(prop={'size':6})
        fig.savefig(str(outputFileNames)+"_SELEX_Analytics_A7A", format='png')
    else:
        return
#TEST
dataAnalysis(20, 15, "thunder", True)
