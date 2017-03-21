#! /usr/bin/env python

import argparse

import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
import seaborn as sns
sns.set(style="white")


import sys
sys.path.append("/home/tc427/software/lib/python3.5/site-packages")
import RNA


def hf(str1, str2):
    s = 0
    for i,j in zip(str1,str2):
        if i != j:
            s +=1
    return s



if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Parse arguments.')
    parser.add_argument('-N', '--Nrounds', type=int, default=9)
    parser.add_argument('-p', '--prefix', type=str, default="hamming_samples_R")
    parser.add_argument('-t', '--target', type=str, default="GUGCCAGAUAGUAAUGGCGU")
    parser.add_argument('-m', '--method', type=str, default="hamming")
    parser.add_argument('-o', '--outfile', type=str, default=None)

    args = parser.parse_args()


    rounds = list()

    fig = plt.figure(figsize=(2.1*args.Nrounds, 10))
    G = gridspec.GridSpec(1, args.Nrounds)
    bins = range(len(args.target))
    for i, g in enumerate(G):
        print("graphing %d" % (i+1))
        samp = list()
        with open(args.prefix+"{:03d}".format(i+1), 'r') as sf:
            for l in sf:                  
                samp.append(l.split()[0])        
        #rounds.append([hf(args.target, i) for i in samp])
        if args.method == "hamming":
            rd = [hf(args.target, i_) for i_ in samp]
        else:
            struct_target = RNA.fold(args.target)[0]
            rd = [RNA.bp_distance(struct_target, RNA.fold(i_)[0]) for i_ in samp]
        ax = plt.subplot(g)
        #sns.distplot(rd, hist=True, rug=True, vertical=True, ax=ax, bins=60)#, hist_kws={"weights": r[:, 1]})#, hist_kws={'bins': "fd"})
        ax.hist(rd, bins=bins, orientation="horizontal")
        ax.set_ylim((0, len(args.target)))
        ax.set_xticklabels([])
        #ax.set_xlabel("Round {:d}".format(i+1))
        ax.set_xlabel("R {:d}".format(i+1))
        if i > 0 and (i+1) < args.Nrounds:
            ax.set_yticklabels([])


    sns.despine(bottom=True)
    fig.suptitle("Distribution of the distance over %d rounds" % args.Nrounds)

    if args.outfile is None:
        plt.show()
    else:
        plt.savefig(args.outfile)
