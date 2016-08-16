import sys
from Aptamers import Aptamers
from Selection import Selection
from Amplification import Amplification
from Mutation import Mutation
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
outputFileNames = str(sys.argv[10])

# Instantiating classes
Apt = Aptamers()
S = Selection()
Amplify = Amplification()
Mut = Mutation()

if(distanceMeasure == "hamming"):
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
            slctdSeqs = S.stochasticHammingSelection_initial(alphabetSet, seqLength, aptamerSeqs, selectionThreshold, initialSeqNum)
            print("selection carried out for R0")
            amplfdSeqs = Amplify.randomPCR_with_ErrorsAndBias_FASTv2(slctdSeqs, seqLength, pcrCycleNum, pcrYield, pcrErrorRate, aptamerSeqs, alphabetSet, distanceMeasure)
            print("amplification carried out for R0")
            outFile = outputFileNames + "_R" + str(r+1)
            nxtRnd = open(outFile, 'w')
            print("writing R1 seqs to file")
            for seqIdx in amplfdSeqs:
                nxtRnd.write(str(int(seqIdx))+'\t'+str(int(amplfdSeqs[seqIdx][0]))+'\t'+str(int(amplfdSeqs[seqIdx][1]))+'\t'+'\n') #write seqIdx, count, distance, and bias...for now
            nxtRnd.close()
        else:
            del(slctdSeqs)
            print("SELEX Round "+str(r+1)+" has started")
            totalSeqNum, uniqSeqNum = utils.seqNumberCounter(amplfdSeqs)
            print("total number of sequences in initial pool = "+str(totalSeqNum))
            print("total number of unique sequences in initial pool = "+str(int(uniqSeqNum)))
            slctdSeqs = S.stochasticHammingSelection(alphabetSet, seqLength, amplfdSeqs, selectionThreshold, uniqSeqNum, totalSeqNum)
            print("Selection carried for R"+str(r+1))
            del(amplfdSeqs)
            amplfdSeqs = Amplify.randomPCR_with_ErrorsAndBias_FASTv2(slctdSeqs, seqLength, pcrCycleNum, pcrYield, pcrErrorRate, aptamerSeqs, alphabetSet, distanceMeasure)
            print("Amplification carried for R"+str(r+1))
            outFile = outputFileNames + "_R" + str(r+1)
            nxtRnd = open(outFile, 'w')
            print("writing R"+str(r+1)+" seqs to file")
            for seqIdx in amplfdSeqs:
                nxtRnd.write(str(int(seqIdx))+'\t'+str(int(amplfdSeqs[seqIdx][0]))+'\t'+str(int(amplfdSeqs[seqIdx][1]))+'\n') #write idx and count for now
            nxtRnd.close()
    print("SELEX completed")

elif(distanceMeasure == "basepair"):
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
            slctdSeqs = S.stochasticBasePairSelection_initial(alphabetSet, seqLength, aptamerSeqs, selectionThreshold, initialSeqNum)
            print("selection carried out for R0")
            amplfdSeqs = Amplify.randomPCR_with_ErrorsAndBias_FASTv2(slctdSeqs, seqLength, pcrCycleNum, pcrYield, pcrErrorRate, aptamerSeqs, alphabetSet, distanceMeasure)
            print("amplification carried out for R0")
            outFile = outputFileNames + "_R" + str(r+1)
            nxtRnd = open(outFile, 'w')
            print("writing R1 seqs to file")
            for seqIdx in amplfdSeqs:
                nxtRnd.write(str(int(seqIdx))+'\t'+str(int(amplfdSeqs[seqIdx][0]))+'\t'+str(int(amplfdSeqs[seqIdx][1]))+'\t'+'\n') #write seqIdx, count, distance, and bias...for now
            nxtRnd.close()
        else:
            del(slctdSeqs)
            print("SELEX Round "+str(r+1)+" has started")
            totalSeqNum, uniqSeqNum = utils.seqNumberCounter(amplfdSeqs)
            print("total number of sequences in initial pool = "+str(totalSeqNum))
            print("total number of unique sequences in initial pool = "+str(int(uniqSeqNum)))
            slctdSeqs = S.stochasticBasePairSelection(alphabetSet, seqLength, amplfdSeqs, selectionThreshold, uniqSeqNum, totalSeqNum)
            print("Selection carried for R"+str(r+1))
            del(amplfdSeqs)
            amplfdSeqs = Amplify.randomPCR_with_ErrorsAndBias_FASTv2(slctdSeqs, seqLength, pcrCycleNum, pcrYield, pcrErrorRate, aptamerSeqs, alphabetSet, distanceMeasure)
            print("Amplification carried for R"+str(r+1))
            outFile = outputFileNames + "_R" + str(r+1)
            nxtRnd = open(outFile, 'w')
            print("writing R"+str(r+1)+" seqs to file")
            for seqIdx in amplfdSeqs:
                nxtRnd.write(str(int(seqIdx))+'\t'+str(int(amplfdSeqs[seqIdx][0]))+'\t'+str(int(amplfdSeqs[seqIdx][1]))+'\n') #write idx and count for now
            nxtRnd.close()
    print("SELEX completed")
else:
    print("Invalid argument for distance measure")

