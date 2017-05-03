import argparse
import os.path
import sys
import random
from Aptamers import Aptamers
from Selection import Selection
from Amplification import Amplification
from Mutation import Mutation
import postprocess
import utils

# Fetch experiment parameters from the settings file
import configparser


def call_post_process(target):
    print("Data post-processing has started...")
    postprocess.dataAnalysis(seqLength, roundNum, outputFileNames, post_process, distanceMeasure, imgformat=img_format)
    # postprocess.dataAnalysis(seqLength, roundNum, "{}_samples".format(outputFileNames), post_process, distanceMeasure, imgformat=img_format)
    postprocess.plot_histo(roundNum, outputFileNames, target, "hamming", "png")
    postprocess.plot_histo(roundNum, "{}_samples".format(outputFileNames), target, "hamming", "png")
    print("Data post-processing is complete.")
    return


parser = argparse.ArgumentParser(description='Parse arguments.')
parser.add_argument('-p', '--postprocess', action='store_true')
parser.add_argument('-s', '--settings', type=str, default="settings.init")

args = parser.parse_args()

if not os.path.exists(args.settings):
    print("Settings file '{}' does not exists, aborting.".format(args.settings))
    sys.exit()

settings = configparser.ConfigParser({"initial_samples": "100000",
                                      "random_seed": "0",
                                      "img_format": "pdf"})
settings.read(args.settings)

aptamerType = settings.get('general', 'selex_type')
aptamerNum = settings.getint('general', 'aptamer_mode')
aptamerSeq = settings.get('general', 'reference_aptamer')
seqLength = settings.getint('general', 'sequence_length')
rng_seed = settings.getint('general', 'random_seed')
roundNum = settings.getint('general', 'number_of_rounds')
outputFileNames = settings.get('general', 'experiment_name')
# how many sampled sequence to output each round, stored in output_samples_Ri.txt
samplingSize = settings.getint('general', 'sampling_size')
post_process = settings.get('general', 'post_process')
img_format = settings.get('general', 'img_format')

# how many sequence to select each round
initialSamples = settings.getint('selectionparams', 'initial_samples')
selectionThreshold = settings.getint('selectionparams', 'scale')
distanceMeasure = settings.get('selectionparams', 'distance')
stringency = settings.getint('selectionparams', 'stringency')

pcrCycleNum = settings.getint('amplificationparams', 'number_of_pcr')
pcrYield = settings.getfloat('amplificationparams', 'pcr_efficiency')
pcrErrorRate = settings.getfloat('amplificationparams', 'pcr_error_rate')

# Instantiating classes
Apt = Aptamers()
S = Selection()
Amplify = Amplification()
Mut = Mutation()

if args.postprocess:
    call_post_process(aptamerSeq)
    sys.exit()

if distanceMeasure not in ("hamming", "basepair", "loop"):
    print("Invalid argument for distance measure")
    sys.exit()

if rng_seed == 0:
    rng_seed = random.randint(0, 2**32)
print("Random seed: {}".format(rng_seed))
random.seed(rng_seed)

assert len(aptamerSeq) == seqLength

# SELEX simulation based on random aptamer assignment, hamming-based definite selection, and
# non-ideal stochastic amplfication with no bias.
for r in range(roundNum):
    if(r == 0):
        if(aptamerType == 'DNA'):
            alphabetSet = 'ACGT'
        elif(aptamerType == 'RNA'):
            alphabetSet = 'ACGU'
        else:
            print("Error: Simulation of %.s aptamers not supported" % aptamerType)
            break
        if aptamerNum > 0:
            aptamerSeqs, initialSeqNum = Apt.optimumAptamerGenerator(aptamerNum, alphabetSet, seqLength)
        else:
            aptamerSeqs = aptamerSeq
            initialSeqNum = len(alphabetSet)**len(aptamerSeq)
        print("optimum sequences have been chosen: %s" % aptamerSeqs)
        print("SELEX Round 1 has started")
        print("total number of sequences in initial library = "+str(initialSeqNum))
        if distanceMeasure == "hamming":
            slctdSeqs = S.stochasticHammingSelection_initial(alphabetSet, seqLength, aptamerSeqs, initialSamples, initialSeqNum, samplingSize, outputFileNames, r, stringency)
        elif distanceMeasure == "basepair":
            slctdSeqs = S.stochasticBasePairSelection_initial(alphabetSet, seqLength, aptamerSeqs, initialSamples, initialSeqNum, samplingSize, outputFileNames, r, stringency)
        else:
            slctdSeqs = S.stochasticLoopSelection_initial(alphabetSet, seqLength, aptamerSeqs, initialSamples, initialSeqNum, samplingSize, outputFileNames, r, stringency)
        print("selection carried out for R1")
        amplfdSeqs = Amplify.randomPCR_with_ErrorsAndBias(slctdSeqs, seqLength, pcrCycleNum, pcrYield, pcrErrorRate, aptamerSeqs, alphabetSet, distanceMeasure)
        print("amplification carried out for R1")
        outFile = outputFileNames + "_R{:03d}".format(r+1)
        nxtRnd = open(outFile, 'w')
        print("writing R1 seqs to file")
        for seqIdx in amplfdSeqs:
            seq = Apt.pseudoAptamerGenerator(seqIdx, alphabetSet, seqLength)
            # write seqIdx, count, distance, and bias...for now
            nxtRnd.write(str(seq)+'\t'+str(int(amplfdSeqs[seqIdx][0]))+'\t'+str(int(amplfdSeqs[seqIdx][1]))+'\t'+'\n')
        nxtRnd.close()
    else:
        del(slctdSeqs)
        print("SELEX Round "+str(r+1)+" has started")
        totalSeqNum, uniqSeqNum = utils.seqNumberCounter(amplfdSeqs)
        print("total number of sequences in initial pool = "+str(totalSeqNum))
        print("total number of unique sequences in initial pool = "+str(int(uniqSeqNum)))
        if(distanceMeasure == "hamming"):
            slctdSeqs = S.stochasticHammingSelection(alphabetSet, seqLength, amplfdSeqs, selectionThreshold, uniqSeqNum, totalSeqNum, samplingSize, outputFileNames, r, stringency)
        elif distanceMeasure == "basepair":
            slctdSeqs = S.stochasticBasePairSelection(alphabetSet, seqLength, amplfdSeqs, selectionThreshold, uniqSeqNum, totalSeqNum, samplingSize, outputFileNames, r, stringency)
        else:
            slctdSeqs = S.stochasticLoopSelection(alphabetSet, seqLength, amplfdSeqs, selectionThreshold, uniqSeqNum, totalSeqNum, samplingSize, outputFileNames, r, stringency)
        print("Selection carried for R"+str(r+1))
        amplfdSeqs = Amplify.randomPCR_with_ErrorsAndBias(slctdSeqs, seqLength, pcrCycleNum, pcrYield, pcrErrorRate, aptamerSeqs, alphabetSet, distanceMeasure)
        print("Amplification carried for R"+str(r+1))
        outFile = outputFileNames + "_R{:03d}".format(r+1)
        nxtRnd = open(outFile, 'w')
        print("writing R"+str(r+1)+" seqs to file")
        for seqIdx in amplfdSeqs:
            seq = Apt.pseudoAptamerGenerator(seqIdx, alphabetSet, seqLength)
            # write idx and count for now
            nxtRnd.write(str(seq)+'\t'+str(int(amplfdSeqs[seqIdx][0]))+'\t'+str(int(amplfdSeqs[seqIdx][1]))+'\n')
        nxtRnd.close()
print("SELEX completed")

if post_process:
    call_post_process(aptamerSeqs)
    print("The simulation has ended.")
else:
    print("The simulation has ended without post-processing.")
