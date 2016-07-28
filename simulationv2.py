import sys
import Aptamers
import Distance
import Selection
import Amplification
import Mutation

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
a = Aptamers.Aptamers()
d = Distance.Distance()
s = Selection.Selection()
amplify = Amplification.Amplification()
mut = Mutation.Mutation()

# SELEX simulation based on random aptamer assignment, hamming-based definite selection, and
# non-ideal stochastic amplfication with no bias. 
for r in range(roundNum):
    if(i==0):
        if(aptamerType = 'DNA'):
            aptamerSeqs, initialSeqNum = a.optimumAptamerGenerator(aptamerNum, 'ACGT', seqLength)
        elif(aptamerType = 'RNA'):
            aptamerSeqs, initialSeqNum = a.optimumAptamerGenerator(aptamerNum, 'ACGU', seqLength)
        else:
            print("Error: Simulation of %.s aptamers not supported" %aptamerType)
            break
 
        
        print("optimum sequences have been chosen")
        print("SELEX Round 1 has started")
        print("total number of sequences in initial library = "+str(initialSeqNum))
        slctdSeqs = s.stochasticHammingSelection(seqLength, aptamerSeqs, selectionThreshold, initialSeqNum, i)
        print("selection carried out for R0")
        amplfdSeqs, mutatedPool = amplify.randomPCR_errorProne_biased(slctdSeqs, seqLength, num_pcr_cycles, pcrYield, errorRate) #change var name returned by amplification
        print("amplification carried out for R0")
        outPool = mut.mutGen(seqLength, mutatedPool, amplfdSeqs, aptamerSeqs, alphabetSet)
        print("Mutation carried out for R0")
        outFile = initPool + "_R" + str(i+1)
        nxtRnd = open(outFile, 'w')
        print("writing R1 seqs to file")
        for seq in outPool:
            nxtRnd.write(seq[0]+'\n') #write seqIdx for now
    else:
        print("SELEX Round "+str(i+1)+" has started")
        junk, totalSeqNum = a.optimumSeqs(nseq, outFile)
        print("totalseqs = "+str(totalSeqNum))
        slctdSeqs = s.stochasticSelection(seqLength, outPool, selectionRate, totalSeqNum)
        print("Selection carried for R"+str(i+1))
        amplfdSeqs, mutatedPool = amplify.ampEffMutDefinite(slctdSeqs, seqLength, num_pcr_cycles, pcrYield, errorRate)
        print("Amplification carried for R"+str(i+1))
        del(slctdSeqs)
        outPool = mut.mutGen(mutatedPool, amplfdSeqs)   
        print("Mutation carried for R"+str(i+1))
        del(amplfdSeqs)
        del(mutatedPool)
        outFile = initPool + "_R" + str(i+1)
        nxtRnd = open(outFile, 'w')
        print("writing R"+str(i+1)+" seqs to file")
        for seq in outPool[seq][0]:
            for j in range(outPool[seq][0]):
                nxtRnd.write(seq+'\n')
