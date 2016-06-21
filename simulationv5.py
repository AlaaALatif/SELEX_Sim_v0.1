import sys
import Aptamers
import Distance
import Selection
import Amplification
import Mutation

nseq=int(sys.argv[1])
initPool=str(sys.argv[2]) 
selectionRate=float(sys.argv[4])
numRnds=int(sys.argv[3])
num_pcr_cycles = int(sys.argv[5])   
pcrYield = float(sys.argv[6])
seqLength = int(sys.argv[7])
errorRate = float(sys.argv[8])
a = Aptamers.Aptamers()
d = Distance.Distance()
s = Selection.Selection()
amplify = Amplification.Amplification()
mut = Mutation.Mutation()

# SELEX simulation based on random aptamer assignment, hamming-based definite selection, and
# non-ideal stochastic amplfication with no bias. 
for i in range(numRnds):
    if(i==0):
        print("SELEX Round 1 has started")
        optimSeqs, totalSeqNum = a.optimumSeqs(nseq, initPool)
        print("optimum seq have been chosen")
        seqs = d.seqsHamming(optimSeqs, initPool)
        print("distances have been calculated for R0")
        print("totalseqs = "+str(totalSeqNum))
        print("selection rate = "+str(selectionRate))
        slctdSeqs = s.stochasticSelection(seqLength, seqs, selectionRate, totalSeqNum)
        del(seqs)
        print("selection carried out for R0")
        print("totalseqs = "+str(totalSeqNum))
        amplfdSeqs, mutatedPool = amplify.ampEffMutDefinite(slctdSeqs, seqLength, num_pcr_cycles, pcrYield, errorRate) #change var name returned by amplification
        print("amplification carried out for R0")
        del(slctdSeqs)        
        outPool = mut.mutGen(mutatedPool, amplfdSeqs)
        print("Mutation carried out for R0")
        del(amplfdSeqs)
        del(mutatedPool)
        outFile = initPool + "_R" + str(i+1)
        nxtRnd = open(outFile, 'w')
        print("writing R1 seqs to file")
        for seq in outPool:
            for j in range(outPool[seq][0]):
                nxtRnd.write(seq+'\n')
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
