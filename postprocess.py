import numpy as np
import matplotlib.pyplot as plt

import RNA

import Distance

D = Distance.Distance()

dpival = 300
params = {'legend.fontsize': 'medium',
          # 'figure.figsize': (15, 5),
          'axes.labelsize': 'small',
          'axes.titlesize': 'small',
          'xtick.labelsize': 'xx-small',
          'ytick.labelsize': 'xx-small'}
plt.rcParams.update(params)


# This generate the main plots from the simulation results
# The plots include changes in total and unique sequence numbers, in average distance
# and changes in the average distance of each affinity group
def dataAnalysis(seqLength, roundNum, outputFileNames, plots, distanceMeasure,
                 aptSeq=None, aptStruct=None, aptLoop=None, imgformat="pdf"):
    avgDist_per_rnd = np.zeros(roundNum)
    weighted_avgDist_per_rnd = np.zeros(roundNum)
    total_seqs_freqs = np.zeros(roundNum)
    uniq_seqs_freqs = np.zeros(roundNum)
    distFreqs = np.zeros((roundNum, seqLength+5))
    weighted_distFreqs = np.zeros((roundNum, seqLength+5))
    for rnd in range(roundNum):
        total_seq_num = 0
        uniq_seq_num = 0
        distance = 0.0
        weighted_distance = 0.0
        with open(outputFileNames + "_R{:03d}".format(rnd+1)) as SELEX_round:
            for line in SELEX_round:
                columns = line.split()
                distance += int(columns[1])
                weighted_distance += int(columns[1])*int(columns[2])
                total_seq_num += int(columns[2])
                uniq_seq_num += 1
                distFreqs[rnd][int(columns[1])] += 1
                weighted_distFreqs[rnd][int(columns[1])] += int(columns[2])
        avgDist_per_rnd[rnd] = distance/uniq_seq_num
        weighted_avgDist_per_rnd[rnd] = weighted_distance/total_seq_num
        total_seqs_freqs[rnd] = total_seq_num
        uniq_seqs_freqs[rnd] = uniq_seq_num
        for i in range(seqLength+5):
            distFreqs[rnd][i] /= uniq_seqs_freqs[rnd]
            weighted_distFreqs[rnd][i] /= total_seqs_freqs[rnd]
    print(weighted_avgDist_per_rnd)
    print(avgDist_per_rnd)
    with open(outputFileNames+"_processed_results", 'w') as p:
        for rnd in range(roundNum):
            p.write(str(int(total_seqs_freqs[rnd]))+'\t')
            p.write(str(int(uniq_seqs_freqs[rnd]))+'\t')
            p.write(str(int(avgDist_per_rnd[rnd]))+'\t')
            p.write(str(int(weighted_avgDist_per_rnd[rnd]))+'\t')
            for l in range(seqLength+1):
                p.write(str(int(distFreqs[rnd][l]))+'\t')
                p.write(str(int(weighted_distFreqs[rnd][l]))+'\t')
            p.write('\n')
    # If the user requested generating plots
    if plots:
        # 30 colors
        # co30 = plt.cm.jet(np.linspace(0,1,30))
        Ncolors = 2**7-1
        co30 = plt.cm.rainbow(np.linspace(0, 1, Ncolors+1))
        # generates indices as far apart from each other as possible
        cidx = [int("{:0>7}".format(bin(Ncolors ^ i)[2:])[::-1], base=2) for i in range(Ncolors)]
        co30 = co30[cidx]
        # linestyles = ['-', '--', '-.', ':']
        # change line below to change color palette used
        # plt.style.use("seaborn-white")
        # If Hamming distances were used
        roundNumAxis = np.linspace(1, roundNum, roundNum)
        fig0, axes = plt.subplots(2, 2, sharex=True)
        plotsList = [total_seqs_freqs, uniq_seqs_freqs, weighted_avgDist_per_rnd, avgDist_per_rnd]
        for i, ax in enumerate(axes.reshape(-1)):
            ax.plot(roundNumAxis, plotsList[i], color='C{}'.format(i))
            if i <= 1:
                ax.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))
        axes[0, 0].set_title('Total sequences')
        axes[0, 0].set_ylabel('Frequency')
        axes[1, 0].set_xlabel('Round Number')
        axes[1, 0].set_ylabel('Average Distance')
        axes[0, 1].set_title('Unique sequences')
        plt.tight_layout()
        fig0.savefig("{}_SELEX_Analytics_distance.{}".format(outputFileNames, imgformat), dpi=dpival)
        fig1, axes = plt.subplots(2, 3, sharex=True)  # , sharey=True)
        for i, ax in enumerate(axes.reshape(-1)):
            for d in range(3):
                ax.plot(roundNumAxis, distFreqs[:, d+(3*i)+1],
                        label='d = '+str(d+(3*i)+1), color=co30[i*3+d])
            ax.ticklabel_format(syle='sci', axis='y', scilimits=(0, 0))
            ax.legend(prop={'size': 6})
        axes[0, 0].set_ylim((-0.1, distFreqs.max()))
        fig1.suptitle('Unique Sequences')
        axes[0, 0].set_ylabel('Fractional Frequency')
        axes[1, 1].set_xlabel('Round Number')
        plt.tight_layout()
        plt.subplots_adjust(top=0.90)
        fig1.savefig("{}_SELEX_Analytics_distFreqs.{}".format(outputFileNames, imgformat), dpi=dpival)
        # weighted fractional sequency plots
        fig2, axes = plt.subplots(2, 3)
        for i, ax in enumerate(axes.reshape(-1)):
            for d in range(3):
                ax.plot(roundNumAxis, weighted_distFreqs[:, d+(3*i)+1],
                        label='d = '+str(d+(3*i)+1), color=co30[i*3+d])
            ax.ticklabel_format(syle='sci', axis='y', scilimits=(0, 0))
            ax.legend(prop={'size': 6})
        axes[0, 0].set_ylim((-0.1, weighted_distFreqs.max()))
        fig2.suptitle('Total Sequences')
        axes[0, 0].set_ylabel('Fractional Frequency')
        axes[1, 1].set_xlabel('Round Number')
        plt.tight_layout()
        plt.subplots_adjust(top=0.90)
        fig2.savefig("{}_SELEX_Analytics_weighted_distFreqs.{}".format(outputFileNames, imgformat), dpi=dpival)
        # probability densities for each round in the same plot, a la aptaSim
        fig3, ax = plt.subplots()
        for i, ax in enumerate(axes.reshape(-1)):
            for d in range(3):
                ax.plot(roundNumAxis, weighted_distFreqs[:, d+(3*i)+1],
                        label='d = '+str(d+(3*i)+1), color=co30[i*3+d])
            ax.ticklabel_format(syle='sci', axis='y', scilimits=(0, 0))
            ax.legend(prop={'size': 6})
        axes[0, 0].set_ylim((-0.1, weighted_distFreqs.max()))
        fig2.suptitle('Total Sequences')
        axes[0, 0].set_ylabel('Fractional Frequency')
        axes[1, 1].set_xlabel('Round Number')
        plt.tight_layout()
        plt.subplots_adjust(top=0.90)
        fig2.savefig("{}_SELEX_Analytics_weighted_distFreqs.{}".format(outputFileNames, imgformat), dpi=dpival)


# for loop only
#            if(aptSeq != None and aptStruct != None and aptLoop != None):
#                dist_matrx = np.zeros((roundNum, seqLength+1, seqLength-6))
#                uniq_dist_matrx = np.zeros((roundNum, seqLength+1, seqLength-6))
#                for rnd in xrange(roundNum):
#                    with open(outputFileNames + "_R{:03d}_dist_components_results".format(rnd+1), 'r') as p:
#                        for line in p:
#                            row = line.split()
#                            loop_dist = float(row[0])
#                            bp_dist = float(row[1])
#                            count = float(row[2])
#                            uniq = float(row[3])
#                            dist_matrx[rnd][int(loop_dist)][int(bp_dist)] += int(count)
#                            uniq_dist_matrx[rnd][int(loop_dist)][int(bp_dist)] += int(uniq)
#                fig3, axes = plt.subplots(2, 4)
#                for i, ax in enumerate(axes.reshape(-1)):
#                    if(i == 0):
#                        cax = ax.imshow(dist_matrx[i+1], interpolation='nearest',
#                                                                cmap=cm.coolwarm)
#                        cbar = fig3.colorbar(cax, ticks=[np.min(uniq_dist_matrx[i+1]),np.max(uniq_dist_matrx[i+1])], ax=ax)
#                        ax.set_title('Round '+str(i+1))
#                    else:
#                        cax = ax.imshow(dist_matrx[i*5], interpolation='nearest',
#                                                                cmap=cm.coolwarm)
#                        cbar = fig3.colorbar(cax, ticks=[np.min(uniq_dist_matrx[i*5]),np.max(uniq_dist_matrx[i*5])], ax=ax)
#                        ax.set_title('Round '+str(i*5))
#                fig3.savefig(str("pp2"+outputFileNames)+"_SELEX_Analytics_dist_heatmap.pdf", format='pdf')
#                fig3.text(0.5, 0.98, 'Total Sequences', ha='center')


def plot_histo(Nrounds, prefix, target, method, imgformat="pdf"):
    plt.style.use("seaborn-white")
    fig, axes = plt.subplots(1, Nrounds, figsize=(2.1*Nrounds, 10), sharey=True)
    plot_histo_(Nrounds, prefix, target, method, axes)
    fig.suptitle("Distribution of the distance over %d rounds" % Nrounds)
    plt.savefig("{}_SELEX_histo.{}".format(prefix, imgformat))


def plot_histo_(Nrounds, prefix, target, method, axes):
    bins = range(len(target))
    for i, ax in enumerate(axes):
        samp = list()
        wsamp = list()
        with open("{}_R{:03d}".format(prefix, i+1), 'r') as sf:
            for l in sf:
                ls = l.split()
                samp.append(ls[0])
                if len(ls) > 2:
                    wsamp.append(int(ls[2]))
                else:
                    wsamp.append(1)
        # rounds.append([D.hamming_func(target, i) for i in samp])
        if method == "hamming":
            rd = [D.hamming_func(target, i_) for i_ in samp]
        else:
            struct_target = RNA.fold(target)[0]
            rd = [RNA.bp_distance(struct_target, RNA.fold(i_)[0]) for i_ in samp]
        ax.hist(rd, bins=bins, normed=True, weights=wsamp, orientation="horizontal", label="weighted")
        # # plot unweighted graph if weights present
        # if sum(wsamp) > len(wsamp):
        #     ax.hist(rd, bins=bins, normed=True, orientation="horizontal", histtype="step", color="C1", linewidth=2, label="unweighted")
        ax.set_ylim((0, len(target)))
        # ax.set_xscale("log")
        ax.set_xticklabels([])
        ax.set_xlabel("R {:d}".format(i+1))
        #if i > 0 and (i+1) < Nrounds:
        #    ax.set_yticklabels([])

    plt.tight_layout()
    plt.subplots_adjust(top=0.90)
