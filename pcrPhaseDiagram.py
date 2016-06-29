import Amplification
import sys

maxInitialCount = int(sys.argv[1])
maxPCRcycles = int(sys.argv[2])
maxPCRyield = float(sys.argv[3])

amp = Amplification.Amplification()
y = float()
while(y < (maxPCRyield)):
    for i in range(maxInitialCount):
        for n in range(maxPCRcycles/5):
            if n==0:
                amp.ampEffBruteTest(i+1, n+1, y+0.1, 1000000)
            else:
            
                amp.ampEffBruteTest(i+1, n*5, y+0.1, 1000000)
y += 0.1

