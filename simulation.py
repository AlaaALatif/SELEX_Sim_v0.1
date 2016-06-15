import sys
import global_levenshteinv2 as gl


nseq=int(sys.argv[1])
init_pool_file=sys.argv[2] 
selection_rate=float(sys.argv[4])                                                     
num_rnds=int(sys.argv[3])                                                           
num_pcr_cycles = int(sys.argv[5])   


for i in range(num_rnds):
    if(i==0):
        optimseqs, totalseqs = gl.optimumseqs(nseq, init_pool_file)
        print("optimum seq have been chosen")
        seqs = gl.seqs_hamming(optimseqs, init_pool_file)
        print("distances have been calculated for R0")
        print("totalseqs = "+str(totalseqs))
        print("selection rate = "+str(selection_rate))
        slctd_seqs = gl.selection(seqs, selection_rate, totalseqs)
        print("selection carried out for R0")
        print("totalseqs = "+str(totalseqs))
        amplfd_seqs = gl.amplification(slctd_seqs, num_pcr_cycles) #change var name returned by amplification
        print("amplification carried out for R0")
        outfile = init_pool_file + "_R" + str(i+1)
        nxt_rnd_file = open(outfile, 'w')
        print("writing R1 seqs to file")
        for seq in amplfd_seqs:
            for j in range(amplfd_seqs[seq][0]):
                nxt_rnd_file.write(seq+'\n')
    else:
        print("beginning R"+str(i+1))
        junk, totalseqs = gl.optimumseqs(nseq, outfile)
        seqs = gl.seqs_hamming(optimseqs, outfile)
        slctd_seqs = gl.selection(seqs, selection_rate, totalseqs)
        amplfd_seqs = gl.amplification(slctd_seqs, num_pcr_cycles)
        outfile = init_pool_file + "_R" + str(i+1)
        nxt_rnd_file = open(outfile, 'w')
        for seq in slctd_seqs[seq][0]:
            for j in range(amplfd_seqs[seq][0]):
                nxt_rnd_file.write(seq+'\n')
