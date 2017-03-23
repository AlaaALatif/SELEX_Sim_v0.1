#! /usr/bin/env python

import argparse

import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
import seaborn as sns


import RNA


def hf(str1, str2):
    s = 0
    for i, j in zip(str1, str2):
        if i != j:
            s += 1
    return s


def plot_histo(Nrounds, prefix, target, method, imgformat="pdf"):
    plt.style.use("seaborn-white")
    # rounds = list()

    fig = plt.figure(figsize=(2.1*Nrounds, 10))
    G = gridspec.GridSpec(1, Nrounds)
    bins = range(len(target))
    for i, g in enumerate(G):
        print("graphing %d" % (i+1))
        samp = list()
        with open("{}_samples_R{:03d}".format(prefix, i+1), 'r') as sf:
            for l in sf:
                samp.append(l.split()[0])
        # rounds.append([hf(target, i) for i in samp])
        if method == "hamming":
            rd = [hf(target, i_) for i_ in samp]
        else:
            struct_target = RNA.fold(target)[0]
            rd = [RNA.bp_distance(struct_target, RNA.fold(i_)[0]) for i_ in samp]
        ax = plt.subplot(g)
        # sns.distplot(rd, hist=True, rug=True, vertical=True, ax=ax, bins=60)#, hist_kws={"weights": r[:, 1]})#, hist_kws={'bins': "fd"})
        ax.hist(rd, bins=bins, orientation="horizontal")
        ax.set_ylim((0, len(target)))
        # ax.set_xscale("log")
        ax.set_xticklabels([])
        ax.set_xlabel("R {:d}".format(i+1))
        if i > 0 and (i+1) < Nrounds:
            ax.set_yticklabels([])

    sns.despine(bottom=True)
    fig.suptitle("Distribution of the distance over %d rounds" % Nrounds)

    plt.tight_layout()
    plt.subplots_adjust(top=0.90)
    plt.savefig("{}_SELEX_histo.{}".format(prefix, imgformat))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Parse arguments.')
    parser.add_argument('-N', '--Nrounds', type=int, default=9)
    parser.add_argument('-p', '--prefix', type=str, default="hamming_samples_R")
    parser.add_argument('-t', '--target', type=str, default="GUGCCAGAUAGUAAUGGCGU")
    parser.add_argument('-m', '--method', type=str, default="hamming")
    parser.add_argument('-f', '--imgformat', type=str, default="png")

    args = parser.parse_args()

    plot_histo(args.Nrounds, args.prefix, args.target, args.method, args.imgformat)
