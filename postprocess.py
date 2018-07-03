import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

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
    pstats = pd.DataFrame(columns=range(roundNum), index=["total", "unique", "avdist", "wavdist"])
    pdf = pd.DataFrame(columns=range(roundNum))
    wpdf = pd.DataFrame(columns=range(roundNum))
    for rnd in range(roundNum):
        data = pd.read_table("{}_R{:03d}".format(outputFileNames, rnd+1), names=["seq", "dist", "count"])
        data["wdist"] = data["dist"]*data["count"]
        bin_edges = range(data["dist"].min(), data["dist"].max())
        c, v = np.histogram(data["dist"], bins=bin_edges, density=True)
        wc, wv = np.histogram(data["dist"], bins=bin_edges, density=True, weights=data["count"])
        for vi, ci in zip(v, c):
            pdf.loc[vi, rnd] = ci
        for vi, ci in zip(wv, wc):
            wpdf.loc[vi, rnd] = ci
        pstats.loc["avdist", rnd] = data["dist"].mean()
        pstats.loc["wavdist", rnd] = data["wdist"].sum() / data["count"].sum()
        pstats.loc["total", rnd] = data["count"].sum()
        pstats.loc["unique", rnd] = data["count"].sum()
    pdf[pdf.isnull()] = 0
    pdf.sort_index(inplace=True)
    wpdf[wpdf.isnull()] = 0
    wpdf.sort_index(inplace=True)
    print(pstats)
    with open(outputFileNames+"_stats.csv", 'w') as p:
        pstats.to_csv(p)
    with open(outputFileNames+"_dists.csv", 'w') as p:
        pdf.to_csv(p)
        wpdf.to_csv(p)
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
        fig0, axes = plt.subplots(2, 2, sharex=True)
        plotsList = ["total", "unique", "wavdist", "avdist"]
        pstats = pstats.T
        for i, ax in enumerate(axes.reshape(-1)):
            ax.plot(pstats.index, pstats[plotsList[i]], color='C{}'.format(i))
            if i <= 1:
                ax.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))
        axes[0, 0].set_title('Total sequences')
        axes[0, 0].set_ylabel('Counts')
        axes[1, 0].set_xlabel('Round Number')
        axes[1, 0].set_ylabel('Average Distance')
        axes[0, 1].set_title('Unique sequences')
        plt.tight_layout()
        fig0.savefig("{}_SELEX_Analytics_distance.{}".format(outputFileNames, imgformat), dpi=dpival)
        n1 = "{}_SELEX_Analytics_distFreqs.{}".format(outputFileNames, imgformat)
        n2 = "{}_SELEX_Analytics_weighted_distFreqs.{}".format(outputFileNames, imgformat)
        for pdf_, oname, tname in zip([pdf, wpdf], [n1, n2], ["Total", "Unique"]):
            ps = pdf_.sum(axis=1).cumsum()
            idxs = ps[ps < np.percentile(ps, 80)].index
            s = len(idxs)//6
            fig1, axes = plt.subplots(2, 3, sharex=True, sharey=True)
            for i, ax in enumerate(axes.reshape(-1)):
                for si in range(s):
                    idx = idxs[s*i+si]
                    ax.plot(pdf_.columns, pdf_.loc[idx], label='d = {}'.format(idx), color=co30[si])
                ax.ticklabel_format(syle='sci', axis='y', scilimits=(0, 0))
                ax.legend(prop={'size': 6})
            axes[0, 0].set_ylim((pdf_.min().min(), pdf_.max().max()))
            fig1.suptitle(f"{tname} Sequences")
            axes[0, 0].set_ylabel('Fractional Frequency')
            axes[1, 1].set_xlabel('Round Number')
            plt.tight_layout()
            plt.subplots_adjust(top=0.90)
            fig1.savefig(oname, dpi=dpival)


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


def plot_histo(Nrounds, prefix, target, imgformat="pdf", method=None):
    plt.style.use("seaborn-white")
    fig, axes = plt.subplots(1, Nrounds, figsize=(2.1*Nrounds, 10), sharey=True)
    plot_histo_(Nrounds, prefix, target, axes, method)
    fig.suptitle("Distribution of the distance over %d rounds" % Nrounds)
    plt.savefig("{}_SELEX_histo.{}".format(prefix, imgformat))


def plot_histo_(Nrounds, prefix, target, axes, method=None):
    bins = range(len(target))
    for i, ax in enumerate(axes):
        data = pd.read_table("{}_R{:03d}".format(prefix, i+1), names=["seq", "dist", "count"])
        if method is not None:
            if method == "hamming":
                rd = [D.hamming_func(target, i_) for i_ in data["seq"]]
            else:
                struct_target = RNA.fold(target)[0]
                rd = [RNA.bp_distance(struct_target, RNA.fold(i_)[0]) for i_ in data["seq"]]
        else:
            rd = data["dist"]
        ax.hist(rd, bins=bins, normed=True, weights=data["count"], orientation="horizontal", label="weighted")
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
