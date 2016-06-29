import Amplification
import sys
from mpi4py import MPI

maxInitialCount = int(sys.argv[1])
maxPCRcycles = int(sys.argv[2])
maxPCRyield = float(sys.argv[3])
coreNum = int(sys.argv[4])
rank = MPI.COMM_WORLD.Get_rank()
size = MPI.COMM_WORLD.Get_rank


amp = Amplification.Amplification()
y = float()
y = 0.0
if rank == 0:
    for i in range(maxInitialCount):
        for n in range(maxPCRcycles):
            amp.ampEffBruteTest(i+1, n+1, y+0.1, 1000000)

if rank == 1:
    for i in range(maxInitialCount):
        for n in range(maxPCRcycles):
            amp.ampEffBruteTest(i+1, n+1, y+0.2, 1000000)


if rank == 2:
    for i in range(maxInitialCount):
        for n in range(maxPCRcycles):
            amp.ampEffBruteTest(i+1, n+1, y+0.3, 1000000)
if rank == 3:
    for i in range(maxInitialCount):
        for n in range(maxPCRcycles):
                amp.ampEffBruteTest(i+1, n+1, y+0.4, 1000000)
if rank == 4:
    for i in range(maxInitialCount):
        for n in range(maxPCRcycles):
            amp.ampEffBruteTest(i+1, n+1, y+0.5, 1000000)
if rank == 5:
    for i in range(maxInitialCount):
        for n in range(maxPCRcycles):
            amp.ampEffBruteTest(i+1, n+1, y+0.6, 1000000)
if rank == 6:
    for i in range(maxInitialCount):
        for n in range(maxPCRcycles):
            amp.ampEffBruteTest(i+1, n+1, y+0.7, 1000000)
if rank == 7:
    for i in range(maxInitialCount):
        for n in range(maxPCRcycles):
            amp.ampEffBruteTest(i+1, n+1, y+0.8, 1000000)

if rank == 8:
    for i in range(maxInitialCount):
        for n in range(maxPCRcycles):
            amp.ampEffBruteTest(i+1, n+1, y+0.9, 1000000)
   
if rank == 9:
    for i in range(maxInitialCount):
        for n in range(maxPCRcycles):
            amp.ampEffBruteTest(i+1, n+1, y+1, 1000000)



