import numpy as np
from scipy.interpolate import spline
import matplotlib.pyplot as plt

import Distance

d = Distance.Distance()

#This calculates the overall average bias for total and unique sequences
#for the given round
def bias_avg(seqFile, seqLength):
    bias = 0
    w_bias = 0
    totalSeqs = 0
    uniqSeqs = 0
    with open(seqFile, 'r') as f:
        for line in f:
            row = line.split()
            seq = row[0]
            bias += d.bias_func(seq, seqLength)
            w_bias += int(row[1])*d.bias_func(seq, seqLength)
            totalSeqs += int(row[1])
            uniqSeqs += 1
    avg_bias = bias/uniqSeqs
    weighted_avg_bias = w_bias/totalSeqs
    return weighted_avg_bias, avg_bias

#This calculates the average bias for each affinity group in the given round
def bias_avg_per_dist(seqFile, seqLength):
    bias_per_dist = np.zeros(seqLength+5)
    w_bias_per_dist = np.zeros(seqLength+5)
    totalSeqs_per_dist = np.zeros(seqLength+5)
    uniqSeqs_per_dist = np.zeros(seqLength+5)
    with open(seqFile, 'r') as f:
        for line in f:
            row = line.split()
            #grab sequence
            seq = row[0]
            #grab frequency
            freq = int(row[1])
            #grab distance
            dist = int(row[2])
            uniqSeqs_per_dist[dist] += 1
            totalSeqs_per_dist[dist] += freq
            #add bias to corresponding dist
            bias_per_dist[dist] += d.bias_func(seq, seqLength)
            w_bias_per_dist[dist] += freq*d.bias_func(seq, seqLength)
    #calculate averages
    for dist in range(seqLength+5):
        if(uniqSeqs_per_dist[dist] > 0):
            bias_per_dist[dist] /= uniqSeqs_per_dist[dist]
            w_bias_per_dist[dist] /= totalSeqs_per_dist[dist]
    return w_bias_per_dist[:seqLength+1], bias_per_dist[:seqLength+1]


#bias_avg("window_R14", 20)

#This plots the theoretical distribution of Hamming distances
def seq_div_hamm(seqLength, alphabetSet):
    uniqSeqNum_per_dist = np.zeros(seqLength+1)
    for h in range(seqLength+1):
        uniqSeqNum_per_dist[h] = (len(alphabetSet)-1)**(h)*comb(seqLength, h)
    hammDistAxis = np.linspace(0, seqLength, seqLength+1)
    hammDistAxis_smooth = np.linspace(0, seqLength, 200)
    uniqSeqNum_smooth = spline(hammDistAxis, uniqSeqNum_per_dist, hammDistAxis_smooth)
    fig, ax = plt.subplots(1,1)
    ax.plot(hammDistAxis_smooth, uniqSeqNum_smooth)
    ax.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    fig.text(0.5, 0.95, 'Unique Sequences', ha='center')
    fig.text(0.5, 0.04, 'Hamming Distance', ha='center')
    fig.text(0.04, 0.5, 'Frequency', va='center', rotation='vertical')
    fig.savefig("SELEX_Analytics_seqDiv_20nt.png")
    return 0
#seq_div_hamm(20, 'ACGT')

#This plots the range of distances for the Hamming, BP, and
#Loop-based metric. The scale defines the number of samples
#to use to construct the distribution 
def distance_range(scale, ref_seq, seqLength, alphabetSet):
    ref_struct = fold(ref_seq)[0]
    ref_loop = apt_loopFinder(ref_seq, ref_struct)
    hamm_dist_array = np.zeros(int(seqLength*1.5))
    bp_dist_array = np.zeros(int(seqLength*1.5))
    loop_dist_array = np.zeros(int(seqLength*1.5))
    randIdxs = random.randint(0, 4**(20)-1, size=scale)
    for i in range(scale):
        randIdx = randIdxs[i]
        randSeq = apt.pseudoAptamerGenerator(randIdx, alphabetSet, seqLength)
        randHammDist = d.hamming_func(randSeq, ref_seq)
        randbpDist = d.bp_func(ref_struct, randSeq)
        randLoopDist = d.loop_func(ref_seq, ref_struct, ref_loop, randSeq, seqLength)
        hamm_dist_array[randHammDist] += 1
        bp_dist_array[randbpDist] += 1
        loop_dist_array[randLoopDist] += 1
    for dist in range(int(seqLength*1.5)):
        hamm_dist_array[dist] /= scale
        bp_dist_array[dist] /= scale
        loop_dist_array[dist] /= scale
    fig, axis = plt.subplots(1,1)
    distAxis = np.linspace(0, int(seqLength+9), int(seqLength+10))
    distAxis_smooth = np.linspace(0, int(seqLength+9), 200)
    hamm_dist_smooth = spline(distAxis, hamm_dist_array, distAxis_smooth)
    bp_dist_smooth = spline(distAxis, bp_dist_array, distAxis_smooth)
    loop_dist_smooth = spline(distAxis, loop_dist_array, distAxis_smooth)
    axis.plot(distAxis_smooth, hamm_dist_smooth, label='Hamming')
    axis.plot(distAxis_smooth, bp_dist_smooth, label='Base-Pair')
    axis.plot(distAxis_smooth, loop_dist_smooth, label='Loop')
    axis.set_xlim([0, 25])
    axis.set_ylim([0, 0.4])
    axis.legend()
    fig.text(0.5, 0.04, 'Distance', ha='center')
    fig.text(0.04, 0.5, 'Fractional Frequency', va='center', rotation='vertical')
    fig.text(0.5, 0.95, 'Distance Distributions', ha='center')
    fig.savefig("SELEX_Analytics_distance_distributions.png")
    return hamm_dist_array

#Identify the aptamer candidate with the highest frequency in the final pool
#outputs the candidate's structure in .svg format
def aptamer_structs(fileNames, seqLength, roundNum, rounds='final'):
    if(rounds == 'final'):
        top_seq_info = [0,0]
        with open(fileNames+"_R"+str(roundNum), 'r') as f:
            for line in f:
                row = line.split()
                seq = str(row[0])
                count = int(row[1])
                if(count > top_seq_info[1]):
                    top_seq_info[0] = seq
                    top_seq_info[1] = count
        with open(fileNames+"_R"+str(roundNum)+"_topstructure_info", 'w') as f:
            seq = top_seq_info[0]
            seq_struct = fold(seq)[0]
            seq_mfe = fold(seq)[1]
            seq_count = top_seq_info[1]
            f.write(seq+'\t'+seq_struct+'\t'+str(seq_mfe)+'\t'+str(seq_count)+'\n')
        svg_rna_plot(seq, seq_struct, fileNames+"_R"+str(roundNum)+"_topstructure.svg")
        return 0
    elif(rounds == 'all'):
        top_seqs_info = []
        for rnd in range(roundNum):
            with open(fileNames+"_R"+str(rnd+1), 'r') as f:
                for line in f:
                    row = line.split()
                    seq = str(row[0])
                    count = int(row[1])
                    if(count > top_seq_info[1]):
                        top_seqs_info.append([seq,count])
        with open(fileNames+"_R"+str(roundNum)+"_topstructures_info", 'w') as f:
            for rnd in range(roundNum):
                seq = top_seqs_info[rnd][0]
                seq_struct = fold(seq)[0]
                seq_mfe = fold(seq)[1]
                seq_count = top_seqs_info[rnd][1]
                f.write(seq+'\t'+seq_struct+'\t'+str(seq_mfe)+'\t'+str(seq_count)+'\n')
                svg_rna_plot(seq, seq_struct, fileNames+"_R"+str(rnd+1)+"_topstructure.svg")
        return 0
    else:
        print("invalid option for string varible rounds. Exiting...")

#aptamer_structs('he4_hamm_small', 20, 40, 'final')


#identify the aptamer candidates with the highest affinity in the final pool
#outputs the structure of the candidate in .svg format
def aptamer_structs_aff(fileNames, seqLength, roundNum, rounds='final'):
    if(rounds == 'final'):
        top_seq_info = [0,0,np.infty]
        with open(fileNames+"_R"+str(roundNum), 'r') as f:
            for line in f:
                row = line.split()
                seq = str(row[0])
                count = int(row[1])
                dist = int(row[2])
                if(dist < top_seq_info[2]):
                    top_seq_info[0] = seq
                    top_seq_info[1] = count
                    top_seq_info[2] = dist
        with open(fileNames+"_R"+str(roundNum)+"_affstructure_info", 'w') as f:
            seq = top_seq_info[0]
            seq_struct = fold(seq)[0]
            seq_mfe = fold(seq)[1]
            seq_count = top_seq_info[1]
            seq_dist = top_seq_info[2]
            f.write(seq+'\t'+seq_struct+'\t'+str(seq_mfe)+'\t'+str(seq_count)+'\t'+str(seq_dist)+'\n')
        svg_rna_plot(seq, seq_struct, fileNames+"_R"+str(roundNum)+"_affstructure.svg")
        return 0
    elif(rounds == 'all'):
        top_seqs_info = []
        for rnd in range(roundNum):
            with open(fileNames+"_R"+str(rnd+1), 'r') as f:
                for line in f:
                    row = line.split()
                    seq = str(row[0])
                    count = int(row[1])
                    dist = int(row[2])
                    if(dist > top_seq_info[2]):
                        top_seqs_info.append([seq,count,dist])
        with open(fileNames+"_R"+str(roundNum)+"_affstructures_info", 'w') as f:
            for rnd in range(roundNum):
                seq = top_seqs_info[rnd][0]
                seq_struct = fold(seq)[0]
                seq_mfe = fold(seq)[1]
                seq_count = top_seqs_info[rnd][1]
                seq_dist = top_seqs_info[rnd][2]
                f.write(seq+'\t'+seq_struct+'\t'+str(seq_mfe)+'\t'+str(seq_count)++'\t'+str(seq_dist)+'\n')
                svg_rna_plot(seq, seq_struct, fileNames+"_R"+str(rnd+1)+"_affstructure.svg")
        return 0
    else:
        print("invalid option for string varible rounds. Exiting...")




#This generate a plot of the average bias of total
#and unique sequences during selex
def generate_bias_plot(fileNames, roundNum, seqLength):
    weighted_bias_per_rnd = np.zeros(roundNum)
    bias_per_rnd = np.zeros(roundNum)
    for rnd in range(roundNum):
        weighted_bias_per_rnd[rnd], bias_per_rnd[rnd] = bias_avg(fileNames+"_R{:03d}".format(rnd+1), seqLength)
    roundNumAxis = np.linspace(0, roundNum, roundNum)
    plotsList = [weighted_bias_per_rnd, bias_per_rnd]
    fig0, axes = plt.subplots(1, 2, sharex=True, sharey=True)
    for i, ax in enumerate(axes.reshape(-1)):
        ax.plot(roundNumAxis, plotsList[i], color="C{}".format(i))
        ax.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
        ax.set_xlabel('Round Number')
    axes[0].set_ylabel('Average Bias')
    axes[0].set_title('Total Sequences')
    axes[1].set_title('Unique Sequences')
    fig0.savefig(str(fileNames)+"_SELEX_Analytics_bias.png", dpi=300)

#generate_bias_plot('he4_loop_small', 40, 20)

#This generates a plot of the average bias of each affinity group
#during selex
def generate_bias_per_dist_plot(fileNames, roundNum, seqLength, distance):
    weighted_bias_per_rnd = np.zeros((seqLength+1, roundNum))
    bias_per_rnd = np.zeros((seqLength+1, roundNum))
    if(distance=="hamming"):
        for rnd in range(roundNum):
            weighted_bias_per_rnd[:, rnd], bias_per_rnd[:, rnd] = bias_avg_per_dist(fileNames+"_R{:03d}".format(rnd+1), seqLength)
        roundNumAxis = np.linspace(0, roundNum, roundNum)
        roundNumAxis_smooth = np.linspace(0, roundNum, 200)
        y_smooth = np.zeros(roundNum)
        plotsList = [weighted_bias_per_rnd[1:seqLength-1], bias_per_rnd[1:seqLength-1]]
        fig0, axes = plt.subplots(2, 1)
        cm = plt.cm.gist_ncar
        colors = [cm(i) for i in np.linspace(0, 0.9, seqLength-1)]
        for i, ax in enumerate(axes.reshape(-1)):
            for d in range(seqLength-2):
                y_smooth = spline(roundNumAxis, plotsList[i][d], roundNumAxis_smooth)
                ax.plot(roundNumAxis_smooth, y_smooth, color=colors[d], label=str(d+1))
                ax.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
            if(i==0):
                ax.legend(loc=2, ncol=4, prop={'size':5})
        fig0.text(0.53, 0.02, 'Round Number', ha='center')
        fig0.text(0.5, 0.96, 'Total Sequences', ha='center')
        fig0.text(0.5, 0.48, 'Unique Sequences', ha='center')
        fig0.text(0.01, 0.5, 'Average Bias', va='center', rotation='vertical')
        fig0.text(0.07, 0.98, '(a)', ha='center')
        fig0.text(0.07, 0.5, '(b)', ha='center')
        fig0.savefig(str(fileNames)+"_SELEX_Analytics_biasDist.png")
        return fig0
    elif(distance=="basepair"):
        for rnd in range(roundNum):
            weighted_bias_per_rnd[:, rnd], bias_per_rnd[:, rnd] = bias_avg_per_dist(fileNames+"_R"+str(rnd+1), seqLength)
        roundNumAxis = np.linspace(0, roundNum, roundNum)
        roundNumAxis_smooth = np.linspace(0, roundNum, 200)
        y_smooth = np.zeros(roundNum)
        plotsList = [weighted_bias_per_rnd[:seqLength-8], bias_per_rnd[:seqLength-8]]
        fig0, axes = plt.subplots(2, 1)
        cm = plt.cm.gist_ncar
        colors = [cm(i) for i in np.linspace(0, 0.9, seqLength-8)]
        for i, ax in enumerate(axes.reshape(-1)):
            for d in range(seqLength-8):
                y_smooth = spline(roundNumAxis, plotsList[i][d], roundNumAxis_smooth)
                ax.plot(roundNumAxis_smooth, y_smooth, color=colors[d], label=str(d))
                ax.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
            if(i==0):
                ax.legend(loc=2, ncol=4, prop={'size':5})
        fig0.text(0.53, 0.02, 'Round Number', ha='center')
        fig0.text(0.5, 0.96, 'Total Sequences', ha='center')
        fig0.text(0.5, 0.48, 'Unique Sequences', ha='center')
        fig0.text(0.01, 0.5, 'Average Bias', va='center', rotation='vertical')
        fig0.text(0.07, 0.98, '(a)', ha='center')
        fig0.text(0.07, 0.5, '(b)', ha='center')
        fig0.savefig(str(fileNames)+"_SELEX_Analytics_biasDist", format='pdf')
        return fig0
    elif(distance=="loop"):
        for rnd in range(roundNum):
            weighted_bias_per_rnd[:, rnd], bias_per_rnd[:, rnd] = bias_avg_per_dist(fileNames+"_R"+str(rnd+1), seqLength)
        roundNumAxis = np.linspace(0, roundNum, roundNum)
        roundNumAxis_smooth = np.linspace(0, roundNum, 200)
        y_smooth = np.zeros(roundNum)
        plotsList = [weighted_bias_per_rnd[:seqLength+1], bias_per_rnd[:seqLength+1]]
        fig0, axes = plt.subplots(2, 1)
        cm = plt.cm.gist_ncar
        colors = [cm(i) for i in np.linspace(0.1, 0.9, seqLength+1)]
        for i, ax in enumerate(axes.reshape(-1)):
            for d in range(seqLength+1):
                y_smooth = spline(roundNumAxis, plotsList[i][d], roundNumAxis_smooth)
                ax.plot(roundNumAxis_smooth, y_smooth, color=colors[d], label=str(d))
                ax.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
            if(i==0):
                ax.legend(loc=2, ncol=4, prop={'size':5})
        fig0.text(0.53, 0.02, 'Round Number', ha='center')
        fig0.text(0.5, 0.96, 'Total Sequences', ha='center')
        fig0.text(0.5, 0.48, 'Unique Sequences', ha='center')
        fig0.text(0.01, 0.5, 'Average Bias', va='center', rotation='vertical')
        fig0.text(0.07, 0.98, '(a)', ha='center')
        fig0.text(0.07, 0.5, '(b)', ha='center')
        fig0.savefig(str(fileNames)+"_SELEX_Analytics_biasDist", format='pdf')
        return fig0
    else:
        print("Invalid distance metric. Exiting...")
        return 0

#fig = generate_bias_per_dist_plot('he4_loop_small', 40, 20, "loop")
