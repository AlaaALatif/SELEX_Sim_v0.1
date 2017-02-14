import sys
from Aptamers import Aptamers
from Selection import Selection
from Amplification import Amplification
from Mutation import Mutation
from postprocess import dataAnalysis
import utils

aptamerType = str(sys.argv[1])
aptamerNum = int(sys.argv[2])
seqLength = int(sys.argv[3])
selectionThreshold = long(sys.argv[4])
distanceMeasure = str(sys.argv[5])
roundNum = int(sys.argv[6])
pcrCycleNum = int(sys.argv[7])
pcrYield = float(sys.argv[8])
pcrErrorRate = float(sys.argv[9])
samplingSize = int(sys.argv[10])
outputFileNames = str(sys.argv[11])
post_process = bool(sys.argv[12])

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
    if(r==0):
        if(aptamerType == 'DNA'):
            alphabetSet = 'ACGT'
            aptamerSeqs, initialSeqNum = Apt.optimumAptamerGenerator(aptamerNum, alphabetSet, seqLength)
        elif(aptamerType == 'RNA'):
            alphabetSet = 'ACGU'
            aptamerSeqs, initialSeqNum = Apt.optimumAptamerGenerator(aptamerNum, alphabetSet, seqLength)
        else:
            print("Error: Simulation of %.s aptamers not supported" %aptamerType)
            break
        print("optimum sequences have been chosen")
        print("SELEX Round 1 has started")
        print("total number of sequences in initial library = "+str(initialSeqNum))
        if(distanceMeasure == "hamming"):
            slctdSeqs = S.stochasticHammingSelection_initial(alphabetSet, seqLength, aptamerSeqs, selectionThreshold, initialSeqNum, samplingSize, outputFileNames, r)
        else:
            slctdSeqs = S.stochasticBasePairSelection_initial(alphabetSet, seqLength, aptamerSeqs, selectionThreshold, initialSeqNum, samplingSize, outputFileNames, r)
        print("selection carried out for R1")
        amplfdSeqs = Amplify.randomPCR_with_ErrorsAndBias(slctdSeqs, seqLength, pcrCycleNum, pcrYield, pcrErrorRate, aptamerSeqs, alphabetSet, distanceMeasure)
        print("amplification carried out for R1")
        outFile = outputFileNames + "_R" + str(r+1)
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
        outFile = outputFileNames + "_R" + str(r+1)
        nxtRnd = open(outFile, 'w')
        print("writing R"+str(r+1)+" seqs to file")
        for seqIdx in amplfdSeqs:
            seq = Apt.pseudoAptamerGenerator(seqIdx, alphabetSet, seqLength)
            nxtRnd.write(str(seq)+'\t'+str(int(amplfdSeqs[seqIdx][0]))+'\t'+str(int(amplfdSeqs[seqIdx][1]))+'\n') #write idx and count for now
        nxtRnd.close()
print("SELEX completed")

if(post_process==True):
    print("Data post-processing has started...")
    dataAnalysis(seqLength, roundNum, outputFileNames, post_process, distanceMeasure)
    print("Data post-processing is complete")
    print("The Simulation has ended")
else:
    print("The Simulation has ended without post-processing")

