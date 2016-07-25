import sys
import Aptamers
from mpi4py import MPI
rank = MPI.COMM_WORLD.Get_rank()
size = MPI.COMM_WORLD.Get_size()


ap = Aptamers.Aptamers()
seqLength = int(sys.argv[1])
scale = int(sys.argv[2])
outFile = str(sys.argv[3])


poolFrac_iter = ap.aptamerGenerator('ATCG', seqLength, 0, 1000000000, outFile+str(4))


#if rank == 0:
    #poolFrac_iter = ap.aptamerGenerator('ATCG', seqLength, (rank*scale/size), ((rank+1)*scale/size), (outFile+str(0)))

#if rank == 1:
    #poolFrac_iter = ap.aptamerGenerator('ATCG', seqLength, (rank*scale/size), ((rank+1)*scale/size), (outFile+str(1)))

#if rank == 2:
    #poolFrac_iter = ap.aptamerGenerator('ATCG', seqLength, (rank*scale/size), ((rank+1)*scale/size), (outFile+str(2)))

#if rank == 3:
    #poolFrac_iter = ap.aptamerGenerator('ATCG', seqLength, (rank*scale/size), ((rank+1)*scale/size), (outFile+str(3)))


#TEST AREA
#for core in range(size):
    #if rank == core:
        #print("Hello from rank "+str(rank))








