import sys
import Aptamers
import Distance
import Selection
import Amplification

nseq=int(sys.argv[1])
init_pool_file=str(sys.argv[2]) 
selection_rate=float(sys.argv[4])
num_rnds=int(sys.argv[3])
num_pcr_cycles = int(sys.argv[5])   
pcr_yield = float(sys.argv[6])
seqLength = int(sys.argv[7])
a = Aptamers.Aptamers()
d = Distance.Distance()
s = Selection.Selection()
amplify = Amplification.Amplification()

# SELEX simulation based on random aptamer assignment, hamming-based stochastic selection, and
# stochastic amplfication with no bias or mutations. 
for i in range(num_rnds):
    if(i==0):
        optimseqs, totalseqs = a.optimumSeqs(nseq, init_pool_file)
        print("optimum seq have been chosen")
        seqs = d.seqsHamming(optimseqs, init_pool_file)
        print("distances have been calculated for R0")
        print("totalseqs = "+str(totalseqs))
        print("selection rate = "+str(selection_rate))
        slctd_seqs = s.stochasticSelection(seqLength, seqs, selection_rate, totalseqs)
        print("selection carried out for R0")
        print("totalseqs = "+str(totalseqs))
        amplfd_seqs = amplify.ampEfficiencyStochastic(slctd_seqs, num_pcr_cycles, pcr_yield) #change var name returned by amplification
        print("amplification carried out for R0")
        outfile = init_pool_file + "_R" + str(i+1)
        nxt_rnd_file = open(outfile, 'w')
        print("writing R1 seqs to file")
        for seq in amplfd_seqs:
            for j in range(amplfd_seqs[seq][0]):
                nxt_rnd_file.write(seq+'\n')
    else:
        print("beginning R"+str(i+1))
        junk, totalseqs = a.optimumSeqs(nseq, outfile)
        slctd_seqs = s.stochasticSelection(seqLength, amplfd_seqs, selection_rate, totalseqs)
        amplfd_seqs = amplify.ampEfficiencyStochastic(slctd_seqs, num_pcr_cycles, pcr_yield)
        outfile = init_pool_file + "_R" + str(i+1)
        nxt_rnd_file = open(outfile, 'w')
        for seq in amplfd_seqs[seq][0]:
            for j in range(amplfd_seqs[seq][0]):
                nxt_rnd_file.write(seq+'\n')
