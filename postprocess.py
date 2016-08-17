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
            fig1 = plt.figure()
            ax1 = fig1.add_subplot(321)
            ax2 = fig1.add_subplot(322)
            ax3 = fig1.add_subplot(323)
            ax4 = fig1.add_subplot(324)
            ax5 = fig1.add_subplot(325)
            ax6 = fig1.add_subplot(326)
            for i in range(3):
                ax1.plot(roundNumAxis, distFreqs[:,i+1], label='d = '+str(i+1))
                ax2.plot(roundNumAxis, distFreqs[:,i+4], label='d = '+str(i+4))
                ax3.plot(roundNumAxis, distFreqs[:,i+7], label='d = '+str(i+7))
                ax4.plot(roundNumAxis, distFreqs[:,i+11], label='d = '+str(i+11))
                ax5.plot(roundNumAxis, distFreqs[:,i+14], label='d = '+str(i+14))
                ax6.plot(roundNumAxis, distFreqs[:,i+17], label='d = '+str(i+17))
            ax3.plot(roundNumAxis, distFreqs[:,10], label='d = '+str(10))
            colormap = plt.cm.gist_ncar
            colors = [colormap(i) for i in np.linspace(0, 1, 6)]
            for i, j in enumerate(ax1.lines):
                j.set_color(colors[i])
            for i, j in enumerate(ax2.lines):
                j.set_color(colors[i])
            for i, j in enumerate(ax3.lines):
                j.set_color(colors[i])
            for i, j in enumerate(ax4.lines):
                j.set_color(colors[i])
            for i, j in enumerate(ax5.lines):
                j.set_color(colors[i])
            for i, j in enumerate(ax6.lines):
                j.set_color(colors[i])
            ax1.legend(prop={'size':4})
            ax2.legend(prop={'size':4})
            ax3.legend(prop={'size':4})
            ax4.legend(prop={'size':4})
            ax5.legend(prop={'size':4})
            ax6.legend(prop={'size':4})
            fig1.text(0.5, 0.04, 'Round Number', ha='center')
            fig1.text(0.04, 0.5, 'Fractional Sequence Count', va='center', rotation='vertical')
            fig1.savefig(str(outputFileNames)+"_SELEX_Analytics_distFreqs", format='pdf')
            fig2 = plt.figure()
            ax1 = fig2.add_subplot(321)
            ax2 = fig2.add_subplot(322)
            ax3 = fig2.add_subplot(323)
            ax4 = fig2.add_subplot(324)
            ax5 = fig2.add_subplot(325)
            ax6 = fig2.add_subplot(326)
        #cm = plt.get_cmap('gist_rainbow')
	#colormap = plt.cm.gist_ncar
	#plt.gca().set_color_cycle([colormap(i) for i in np.linspace(0, 0.9, seqLength+1)])
        #ax = fig.add_subplot(111)
        #ax.set_color_cycle([cm(1.*i/seqLength+1) for i in range(seqLength+1)])
            for i in range(3):
                ax1.plot(roundNumAxis, weighted_distFreqs[:,i+1], label='d = '+str(i+1))
                ax2.plot(roundNumAxis, weighted_distFreqs[:,i+4], label='d = '+str(i+4))
                ax3.plot(roundNumAxis, weighted_distFreqs[:,i+7], label='d = '+str(i+7))
                ax4.plot(roundNumAxis, weighted_distFreqs[:,i+11], label='d = '+str(i+11))
                ax5.plot(roundNumAxis, weighted_distFreqs[:,i+14], label='d = '+str(i+14))
                ax6.plot(roundNumAxis, weighted_distFreqs[:,i+17], label='d = '+str(i+17))
                ax3.plot(roundNumAxis, weighted_distFreqs[:,10], label='d = '+str(10))
            colormap = plt.cm.gist_ncar
            colors = [colormap(i) for i in np.linspace(0, 1, 6)]
            for i, j in enumerate(ax1.lines):
                j.set_color(colors[i])
            for i, j in enumerate(ax2.lines):
                j.set_color(colors[i])
            for i, j in enumerate(ax3.lines):
                j.set_color(colors[i])
            for i, j in enumerate(ax4.lines):
                j.set_color(colors[i])
            for i, j in enumerate(ax5.lines):
                j.set_color(colors[i])
            for i, j in enumerate(ax6.lines):
                j.set_color(colors[i])
            ax1.legend(prop={'size':4})
            ax2.legend(prop={'size':4})
            ax3.legend(prop={'size':4})
            ax4.legend(prop={'size':4})
            ax5.legend(prop={'size':4})
            ax6.legend(prop={'size':4})
            fig2.savefig(str(outputFileNames)+"_SELEX_Analytics_weighted_distFreqs", format='pdf')
        elif(distanceMeasure=="basepair"):
            roundNumAxis = np.linspace(1, roundNum, roundNum)
            #figures for distance analytics
            fig1 = plt.figure()
            ax1 = fig1.add_subplot(321)
            ax2 = fig1.add_subplot(322)
            ax3 = fig1.add_subplot(323)
            ax4 = fig1.add_subplot(324)
            ax5 = fig1.add_subplot(325)
            ax6 = fig1.add_subplot(326)
        #cm = plt.get_cmap('gist_rainbow')
	#colormap = plt.cm.gist_ncar
	#plt.gca().set_color_cycle([colormap(i) for i in np.linspace(0, 0.9, seqLength+1)])
        #ax = fig.add_subplot(111)
        #ax.set_color_cycle([cm(1.*i/seqLength+1) for i in range(seqLength+1)])
            for i in range(3):
                ax1.plot(roundNumAxis, distFreqs[:,i+1], label='d = '+str(i+1))
                ax2.plot(roundNumAxis, distFreqs[:,i+4], label='d = '+str(i+4))
                ax3.plot(roundNumAxis, distFreqs[:,i+7], label='d = '+str(i+7))
                ax4.plot(roundNumAxis, distFreqs[:,i+11], label='d = '+str(i+11))
                ax5.plot(roundNumAxis, distFreqs[:,i+14], label='d = '+str(i+14))
                ax6.plot(roundNumAxis, distFreqs[:,i+17], label='d = '+str(i+17))
            ax3.plot(roundNumAxis, distFreqs[:,10], label='d = '+str(10))
            colormap = plt.cm.gist_ncar
            colors = [colormap(i) for i in np.linspace(0, 1, 6)]
            for i, j in enumerate(ax1.lines):
                j.set_color(colors[i])
            for i, j in enumerate(ax2.lines):
                j.set_color(colors[i])
            for i, j in enumerate(ax3.lines):
                j.set_color(colors[i])
            for i, j in enumerate(ax4.lines):
                j.set_color(colors[i])
            for i, j in enumerate(ax5.lines):
                j.set_color(colors[i])
            for i, j in enumerate(ax6.lines):
                j.set_color(colors[i])
            ax1.legend(prop={'size':4})
            ax2.legend(prop={'size':4})
            ax3.legend(prop={'size':4})
            ax4.legend(prop={'size':4})
            ax5.legend(prop={'size':4})
            ax6.legend(prop={'size':4})
            fig1.savefig(str(outputFileNames)+"_SELEX_Analytics_distFreqs", format='pdf')
            # weighted fractional sequency plots
            fig2, axes = plt.figure(2, 3)
            colormap = plt.cm.gist_ncar
            colors = [colormap(i) for i in np.linspace(0, 1, 6)]
            for i, ax in enumerate(axes):
                ax = fig2.add_subplot(321+i)
                for d in range(3):
                    ax.plot(roundNumAxis, weighted_distFreqs[:,d+(3*i)+1], label='d = '+str(d+1))
            for i, ax in enumerate(axes):
                ax.legend(prop={'size':4})
                for j, line in enumerate(ax.lines):
                    line.set_color(color[i])
            fig2.savefig(str(outputFileNames)+"_SELEX_Analytics_weighted_distFreqs", format='pdf')
        else:
            return
#TEST
#dataAnalysis(20, 15, "complex", True, "basepair")
