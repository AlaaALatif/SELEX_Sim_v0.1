import sys
import Aptamers
import Selection
import Amplification
import Mutation
import utils

aptamerType = str(sys.argv[1])
aptamerNum = int(sys.argv[2])
seqLength = int(sys.argv[3])
selectionThreshold = long(sys.argv[4])
roundNum = int(sys.argv[5])
pcrCycleNum = int(sys.argv[6])
pcrYield = float(sys.argv[7])
pcrErrorRate = float(sys.argv[8])
outputFileNames = str(sys.argv[9])

# Instantiating classes
Apt = Aptamers.Aptamers()
S = Selection.Selection()
Amplify = Amplification.Amplification()
Mut = Mutation.Mutation()
util = utils.utils()

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

        slctdSeqs = S.stochasticHammingSelection(alphabetSet, seqLength, None, aptamerSeqs, selectionThreshold, initialSeqNum, r)
        print("selection carried out for R0")
        amplfdSeqs, mutatedPool = Amplify.randomPCR_errorProne_biased(slctdSeqs, seqLength, pcrCycleNum, pcrYield, pcrErrorRate) #change var name returned by amplification
        print("amplification carried out for R0")
        outPool = Mut.mutGen(seqLength, mutatedPool, amplfdSeqs, aptamerSeqs, alphabetSet)
        print("Mutation carried out for R0")
        outFile = outputFileNames + "_R" + str(i+1)
        nxtRnd = open(outFile, 'w')
        print("writing R1 seqs to file")
        for seqIdx in outPool:
            nxtRnd.write(str(seqIdx)+'\t'+str(outPool[seqIdx][0])+'\n') #write seqIdx and count for now
        nxtRnd.close()
    else:
        print("SELEX Round "+str(i+1)+" has started")
        totalSeqNum, uniqSeqNum = util.seqNum(outPool)
        print("totalseqs = "+str(totalSeqNum))
        slctdSeqs = S.stochasticSelection(alphabetSet, seqLength, outPool, None, selectionThreshold, totalSeqNum, r)
        print("Selection carried for R"+str(i+1))
        amplfdSeqs, mutatedPool = Amplify.ampEffMutDefinite(slctdSeqs, seqLength, pcrCycleNum, pcrYield, pcrErrorRate)
        print("Amplification carried for R"+str(i+1))
        outPool = Mut.mutGen(mutatedPool, amplfdSeqs)   
        print("Mutation carried for R"+str(i+1))
        outFile = outputFileNames + "_R" + str(i+1)
        nxtRnd = open(outFile, 'w')
        print("writing R"+str(i+1)+" seqs to file")
        for seqIdx in outPool:
            nxtRnd.write(str(seqIdx)+'\t'+str(outPool[seqIdx][0])+'\n') #write idx and count for now
        nxtRnd.close()
