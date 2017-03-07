import sys
from Aptamers import Aptamers
from Selection import Selection
from Amplification import Amplification
from Mutation import Mutation
from postprocess import dataAnalysis
import utils

# Fetch experiment parameters from the settings file
import ConfigParser
settings = ConfigParser.ConfigParser()
settings.read('settings.init')

aptamerType = settings.get('general', 'selex_type')
aptamerNum = settings.getint('general', 'aptamer_mode')
aptamerSeq = settings.get('general', 'reference_aptamer')
seqLength = settings.getint('general', 'sequence_length')
roundNum = settings.getint('general', 'number_of_rounds')
outputFileNames = settings.get('general', 'experiment_name')
# how many sampled sequence to output each round, stored in output_samples_Ri.txt
samplingSize = settings.getint('general', 'sampling_size')
post_process = settings.get('general', 'post_process')

# how many sequence to select each round
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


if distanceMeasure not in ("hamming", "basepair"):
    print("Invalid argument for distance measure")
    sys.exit()

# SELEX simulation based on random aptamer assignment, hamming-based definite selection, and
# non-ideal stochastic amplfication with no bias.
for r in range(roundNum):
    if(r == 0):
        if(aptamerType == 'DNA'):
            alphabetSet = 'ACGT'
            aptamerSeqs, initialSeqNum = Apt.optimumAptamerGenerator(aptamerNum, alphabetSet, seqLength)
        elif(aptamerType == 'RNA'):
            alphabetSet = 'ACGU'
            aptamerSeqs, initialSeqNum = Apt.optimumAptamerGenerator(aptamerNum, alphabetSet, seqLength)
        else:
            print("Error: Simulation of %.s aptamers not supported" % aptamerType)
            break
        print("optimum sequences have been chosen: %s" % aptamerSeqs)
        print("SELEX Round 1 has started")
        print("total number of sequences in initial library = "+str(initialSeqNum))
        if(distanceMeasure == "hamming"):
            slctdSeqs = S.stochasticHammingSelection_initial(alphabetSet, seqLength, aptamerSeqs, selectionThreshold, initialSeqNum, samplingSize, outputFileNames, r)
        else:
            slctdSeqs = S.stochasticBasePairSelection_initial(alphabetSet, seqLength, aptamerSeqs, selectionThreshold, initialSeqNum, samplingSize, outputFileNames, r)
        print("selection carried out for R1")
        amplfdSeqs = Amplify.randomPCR_with_ErrorsAndBias(slctdSeqs, seqLength, pcrCycleNum, pcrYield, pcrErrorRate, aptamerSeqs, alphabetSet, distanceMeasure)
        print("amplification carried out for R1")
        outFile = outputFileNames + "_R{:03d}".format(r+1)
        nxtRnd = open(outFile, 'w')
        print("writing R1 seqs to file")
        for seqIdx in amplfdSeqs:
            seq = Apt.pseudoAptamerGenerator(seqIdx, alphabetSet, seqLength)
            nxtRnd.write(str(seq)+'\t'+str(int(amplfdSeqs[seqIdx][0]))+'\t'+str(int(amplfdSeqs[seqIdx][1]))+'\t'+'\n') #write seqIdx, count, distance, and bias...for now
        nxtRnd.close()
    else:
        del(slctdSeqs)
        print("SELEX Round "+str(r+1)+" has started")
        totalSeqNum, uniqSeqNum = utils.seqNumberCounter(amplfdSeqs)
        print("total number of sequences in initial pool = "+str(totalSeqNum))
        print("total number of unique sequences in initial pool = "+str(int(uniqSeqNum)))
        if(distanceMeasure == "hamming"):
            slctdSeqs = S.stochasticHammingSelection(alphabetSet, seqLength, amplfdSeqs, selectionThreshold, uniqSeqNum, totalSeqNum, samplingSize, outputFileNames, r)
        else:
            slctdSeqs = S.stochasticBasePairSelection(alphabetSet, seqLength, amplfdSeqs, selectionThreshold, uniqSeqNum, totalSeqNum, samplingSize, outputFileNames, r)
        print("Selection carried for R"+str(r+1))
        del(amplfdSeqs)
        amplfdSeqs = Amplify.randomPCR_with_ErrorsAndBias(slctdSeqs, seqLength, pcrCycleNum, pcrYield, pcrErrorRate, aptamerSeqs, alphabetSet, distanceMeasure)
        print("Amplification carried for R"+str(r+1))
        outFile = outputFileNames + "_R{:03d}".format(r+1)
        nxtRnd = open(outFile, 'w')
        print("writing R"+str(r+1)+" seqs to file")
        for seqIdx in amplfdSeqs:
            seq = Apt.pseudoAptamerGenerator(seqIdx, alphabetSet, seqLength)
            nxtRnd.write(str(seq)+'\t'+str(int(amplfdSeqs[seqIdx][0]))+'\t'+str(int(amplfdSeqs[seqIdx][1]))+'\n') #write idx and count for now
        nxtRnd.close()
print("SELEX completed")

if(post_process == True):
    print("Data post-processing has started...")
    dataAnalysis(seqLength, roundNum, outputFileNames, post_process, distanceMeasure)
    print("Data post-processing is complete")
    print("The Simulation has ended")
else:
    print("The Simulation has ended without post-processing")
