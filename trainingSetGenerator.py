## NOTE: put all input params and output components in one numpy
#        matrix. Use gather() to collect local objects of the 
#        matrix. 



import Amplification
import sys
from mpi4py import MPI

maxInitialCount = int(sys.argv[1])
maxPCRcycles = int(sys.argv[2])
maxYield = float(sys.argv[3])
yieldInc = float(sys.argv[4]) #yield increment
maxGaussNum = int(sys.argv[5])

rank = MPI.COMM_WORLD.Get_rank()
size = MPI.COMM_WORLD.Get_size()


amp = Amplification.Amplification()

sampleNum, outputDim = (maxInitialCount*maxPCRcycles*(maxYield/yieldInc)), maxGaussNum


#initialize input vectors

y = float() #initialize yield index
y = 0.0


if rank == 0:
    #initialize local sample index
    j = 0
    #initialize local input vectors
    yieldNumLoc1 = [0 for y in range(sampleNum/size)]
    countNumLoc1 = [0 for i in range(sampleNum/size)]
    cycleNumLoc1 = [0 for n in range(sampleNum/size)]
    #initialize local output matrices
    
    gmmMeansLoc1 = [[0 for d in range(outputDim)] for s in range(int(sampleNum/coreNum))]

    gmmCovarsLoc1 = [[0 for d in range(outputDim)] for s in range(int(sampleNum/coreNum))]

    gmmWeightsLoc1 = [[0 for d in range(outputDim)] for s in range(int(sampleNum/coreNum))]


    for i in range(maxInitialCount):
        for n in range(maxPCRcycles):
            gmmModel = amp.BruteGMM(i+1, n+1, y+yieldInc, 10000)
            yieldNumLoc1[j] = y+yieldInc
            countNumLoc1[j] = i+1
            cycleNumLoc1[j] = n+1
            for i, mu in enumerate(gmmModel.means_):
                gmmMeansLoc1[j][i] = gmmModel.means_[i][0]
                gmmCovarsLoc1[j][i] = gmmModel.covars_[i][0]
                gmmWeightsLoc1[j][i] = gmmModel.weights_[i]
            j+=1 #increment sample index

if rank == 1:
    #initialize local sample index
    j = 0
    #initialize local input vectors
    yieldNumLoc2 = [0 for y in range(sampleNum/size)]
    countNumLoc2 = [0 for i in range(sampleNum/size)]
    cycleNumLoc2 = [0 for n in range(sampleNum/size)]
    #initialize local output matrices
    
    gmmMeansLoc2 = [[0 for d in range(outputDim)] for s in range(int(sampleNum/coreNum))]

    gmmCovarsLoc2 = [[0 for d in range(outputDim)] for s in range(int(sampleNum/coreNum))]

    gmmWeightsLoc2 = [[0 for d in range(outputDim)] for s in range(int(sampleNum/coreNum))]


    for i in range(maxInitialCount):
        for n in range(maxPCRcycles):
            gmmModel = amp.BruteGMM(i+1, n+1, y+(2*yieldInc), 10000)
            yieldNumLoc2[j] = y+(2*yieldInc)
            countNumLoc2[j] = i+1
            cycleNumLoc2[j] = n+1
            for i, mu in enumerate(gmmModel.means_):
                gmmMeansLoc2[j][i] = gmmModel.means_[i][0]
                gmmCovarsLoc2[j][i] = gmmModel.covars_[i][0]
                gmmWeightsLoc2[j][i] = gmmModel.weights_[i]
            j+=1 #increment sample index

if rank == 2:
    #initialize local sample index
    j = 0
    #initialize local input vectors
    yieldNumLoc3 = [0 for y in range(sampleNum/size)]
    countNumLoc3 = [0 for i in range(sampleNum/size)]
    cycleNumLoc3 = [0 for n in range(sampleNum/size)]
    #initialize local output matrices
    
    gmmMeansLoc3 = [[0 for d in range(outputDim)] for s in range(int(sampleNum/coreNum))]

    gmmCovarsLoc3 = [[0 for d in range(outputDim)] for s in range(int(sampleNum/coreNum))]

    gmmWeightsLoc3 = [[0 for d in range(outputDim)] for s in range(int(sampleNum/coreNum))]


    for i in range(maxInitialCount):
        for n in range(maxPCRcycles):
            gmmModel = amp.BruteGMM(i+1, n+1, y+(3*yieldInc), 10000)
            yieldNumLoc3[j] = y+(3*yieldInc)
            countNumLoc3[j] = i+1
            cycleNumLoc3[j] = n+1
            for i, mu in enumerate(gmmModel.means_):
                gmmMeansLoc3[j][i] = gmmModel.means_[i][0]
                gmmCovarsLoc3[j][i] = gmmModel.covars_[i][0]
                gmmWeightsLoc3[j][i] = gmmModel.weights_[i]
            j+=1 #increment sample index

if rank == 3:
    #initialize local sample index
    j = 0
    #initialize local input vectors
    yieldNumLoc4 = [0 for y in range(sampleNum/size)]
    countNumLoc4 = [0 for i in range(sampleNum/size)]
    cycleNumLoc4 = [0 for n in range(sampleNum/size)]
    #initialize local output matrices
    
    gmmMeansLoc4 = [[0 for d in range(outputDim)] for s in range(int(sampleNum/coreNum))]

    gmmCovarsLoc4 = [[0 for d in range(outputDim)] for s in range(int(sampleNum/coreNum))]

    gmmWeightsLoc4 = [[0 for d in range(outputDim)] for s in range(int(sampleNum/coreNum))]


    for i in range(maxInitialCount):
        for n in range(maxPCRcycles):
            gmmModel = amp.BruteGMM(i+1, n+1, y+(4*yieldInc), 10000)
            yieldNumLoc4[j] = y+(4*yieldInc)
            countNumLoc4[j] = i+1
            cycleNumLoc4[j] = n+1
            for i, mu in enumerate(gmmModel.means_):
                gmmMeansLoc4[j][i] = gmmModel.means_[i][0]
                gmmCovarsLoc4[j][i] = gmmModel.covars_[i][0]
                gmmWeightsLoc4[j][i] = gmmModel.weights_[i]
            j+=1 #increment sample index

if rank == 4:
    #initialize local sample index
    j = 0
    #initialize local input vectors
    yieldNumLoc5 = [0 for y in range(sampleNum/size)]
    countNumLoc5 = [0 for i in range(sampleNum/size)]
    cycleNumLoc5 = [0 for n in range(sampleNum/size)]
    #initialize local output matrices
    
    gmmMeansLoc5 = [[0 for d in range(outputDim)] for s in range(int(sampleNum/coreNum))]

    gmmCovarsLoc5 = [[0 for d in range(outputDim)] for s in range(int(sampleNum/coreNum))]

    gmmWeightsLoc5 = [[0 for d in range(outputDim)] for s in range(int(sampleNum/coreNum))]


    for i in range(maxInitialCount):
        for n in range(maxPCRcycles):
            gmmModel = amp.BruteGMM(i+1, n+1, y+(5*yieldInc), 10000)
            yieldNumLoc5[j] = y+(5*yieldInc)
            countNumLoc5[j] = i+1
            cycleNumLoc5[j] = n+1
            for i, mu in enumerate(gmmModel.means_):
                gmmMeansLoc5[j][i] = gmmModel.means_[i][0]
                gmmCovarsLoc5[j][i] = gmmModel.covars_[i][0]
                gmmWeightsLoc5[j][i] = gmmModel.weights_[i]
            j+=1 #increment sample index

if rank == 5:
    #initialize local sample index
    j = 0
    #initialize local input vectors
    yieldNumLoc6 = [0 for y in range(sampleNum/size)]
    countNumLoc6 = [0 for i in range(sampleNum/size)]
    cycleNumLoc6 = [0 for n in range(sampleNum/size)]
    #initialize local output matrices
    
    gmmMeansLoc6 = [[0 for d in range(outputDim)] for s in range(int(sampleNum/coreNum))]

    gmmCovarsLoc6 = [[0 for d in range(outputDim)] for s in range(int(sampleNum/coreNum))]

    gmmWeightsLoc6 = [[0 for d in range(outputDim)] for s in range(int(sampleNum/coreNum))]


    for i in range(maxInitialCount):
        for n in range(maxPCRcycles):
            gmmModel = amp.BruteGMM(i+1, n+1, y+(6*yieldInc), 10000)
            yieldNumLoc6[j] = y+(6*yieldInc)
            countNumLoc6[j] = i+1
            cycleNumLoc6[j] = n+1
            for i, mu in enumerate(gmmModel.means_):
                gmmMeansLoc6[j][i] = gmmModel.means_[i][0]
                gmmCovarsLoc6[j][i] = gmmModel.covars_[i][0]
                gmmWeightsLoc6[j][i] = gmmModel.weights_[i]
            j+=1 #increment sample index

if rank == 6:
    #initialize local sample index
    j = 0
    #initialize local input vectors
    yieldNumLoc7 = [0 for y in range(sampleNum/size)]
    countNumLoc7 = [0 for i in range(sampleNum/size)]
    cycleNumLoc7 = [0 for n in range(sampleNum/size)]
    #initialize local output matrices
    
    gmmMeansLoc7 = [[0 for d in range(outputDim)] for s in range(int(sampleNum/coreNum))]

    gmmCovarsLoc7 = [[0 for d in range(outputDim)] for s in range(int(sampleNum/coreNum))]

    gmmWeightsLoc7 = [[0 for d in range(outputDim)] for s in range(int(sampleNum/coreNum))]


    for i in range(maxInitialCount):
        for n in range(maxPCRcycles):
            gmmModel = amp.BruteGMM(i+1, n+1, y+(7*yieldInc), 10000)
            yieldNumLoc7[j] = y+(7*yieldInc)
            countNumLoc7[j] = i+1
            cycleNumLoc7[j] = n+1
            for i, mu in enumerate(gmmModel.means_):
                gmmMeansLoc7[j][i] = gmmModel.means_[i][0]
                gmmCovarsLoc7[j][i] = gmmModel.covars_[i][0]
                gmmWeightsLoc7[j][i] = gmmModel.weights_[i]
            j+=1 #increment sample index

if rank == 7:
    #initialize local sample index
    j = 0
    #initialize local input vectors
    yieldNumLoc8 = [0 for y in range(sampleNum/size)]
    countNumLoc8 = [0 for i in range(sampleNum/size)]
    cycleNumLoc8 = [0 for n in range(sampleNum/size)]
    #initialize local output matrices
    
    gmmMeansLoc8 = [[0 for d in range(outputDim)] for s in range(int(sampleNum/coreNum))]

    gmmCovarsLoc8 = [[0 for d in range(outputDim)] for s in range(int(sampleNum/coreNum))]

    gmmWeightsLoc8 = [[0 for d in range(outputDim)] for s in range(int(sampleNum/coreNum))]


    for i in range(maxInitialCount):
        for n in range(maxPCRcycles):
            gmmModel = amp.BruteGMM(i+1, n+1, y+(8*yieldInc), 10000)
            yieldNumLoc8[j] = y+(8*yieldInc)
            countNumLoc8[j] = i+1
            cycleNumLoc8[j] = n+1
            for i, mu in enumerate(gmmModel.means_):
                gmmMeansLoc8[j][i] = gmmModel.means_[i][0]
                gmmCovarsLoc8[j][i] = gmmModel.covars_[i][0]
                gmmWeightsLoc8[j][i] = gmmModel.weights_[i]
            j+=1 #increment sample index

if rank == 8:
    #initialize local sample index
    j = 0
    #initialize local input vectors
    yieldNumLoc9 = [0 for y in range(sampleNum/size)]
    countNumLoc9 = [0 for i in range(sampleNum/size)]
    cycleNumLoc9 = [0 for n in range(sampleNum/size)]
    #initialize local output matrices
    
    gmmMeansLoc9 = [[0 for d in range(outputDim)] for s in range(int(sampleNum/coreNum))]

    gmmCovarsLoc9 = [[0 for d in range(outputDim)] for s in range(int(sampleNum/coreNum))]

    gmmWeightsLoc9 = [[0 for d in range(outputDim)] for s in range(int(sampleNum/coreNum))]


    for i in range(maxInitialCount):
        for n in range(maxPCRcycles):
            gmmModel = amp.BruteGMM(i+1, n+1, y+(9*yieldInc), 10000)
            yieldNumLoc9[j] = y+(9*yieldInc)
            countNumLoc9[j] = i+1
            cycleNumLoc9[j] = n+1
            for i, mu in enumerate(gmmModel.means_):
                gmmMeansLoc9[j][i] = gmmModel.means_[i][0]
                gmmCovarsLoc9[j][i] = gmmModel.covars_[i][0]
                gmmWeightsLoc9[j][i] = gmmModel.weights_[i]
            j+=1 #increment sample index

if rank == 9:
    #initialize local sample index
    j = 0
    #initialize local input vectors
    yieldNumLoc10 = [0 for y in range(sampleNum/size)]
    countNumLoc10 = [0 for i in range(sampleNum/size)]
    cycleNumLoc10 = [0 for n in range(sampleNum/size)]
    #initialize local output matrices
    
    gmmMeansLoc10 = [[0 for d in range(outputDim)] for s in range(int(sampleNum/coreNum))]

    gmmCovarsLoc10 = [[0 for d in range(outputDim)] for s in range(int(sampleNum/coreNum))]

    gmmWeightsLoc10 = [[0 for d in range(outputDim)] for s in range(int(sampleNum/coreNum))]


    for i in range(maxInitialCount):
        for n in range(maxPCRcycles):
            gmmModel = amp.BruteGMM(i+1, n+1, y+(10*yieldInc), 10000)
            yieldNumLoc10[j] = y+(10*yieldInc)
            countNumLoc10[j] = i+1
            cycleNumLoc10[j] = n+1
            for i, mu in enumerate(gmmModel.means_):
                gmmMeansLoc10[j][i] = gmmModel.means_[i][0]
                gmmCovarsLoc10[j][i] = gmmModel.covars_[i][0]
                gmmWeightsLoc10[j][i] = gmmModel.weights_[i]
            j+=1 #increment sample index


