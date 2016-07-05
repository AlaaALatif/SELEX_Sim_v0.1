## NOTE: put all input params and output components in one numpy
#        matrix. Use gather() to collect local objects of the 
#        matrix. 


import numpy as np
import Amplification
import sys
from mpi4py import MPI

maxInitialCount = int(sys.argv[1])
maxPCRcycles = int(sys.argv[2])
maxYield = float(sys.argv[3])
yieldInc = float(sys.argv[4]) #yield increment
maxGaussNum = int(sys.argv[5])
outFile = str(sys.argv[6])

rank = MPI.COMM_WORLD.Get_rank()
size = MPI.COMM_WORLD.Get_size()


amp = Amplification.Amplification()

inputDim = 3

sampleNum, outputDim = (maxInitialCount*maxPCRcycles*(maxYield/yieldInc)), (maxGaussNum*3)

dimension = outputDim + inputDim

#initialize input vectors
data = np.zeros((size+1, sampleNum, dimension))

y = float() #initialize yield index
y = 0.0

if rank == 0:
    #initialize local sample index
    j = 0

    for i in range(maxInitialCount):
        for n in range(maxPCRcycles):
            gmmModel = amp.GMMTest(i+1, n+1, 0.1, 10000, maxGaussNum)
            data[rank+1][j][2] = 0.1
            data[rank+1][j][0] = i+1
            data[rank+1][j][1] = n+1
            for k, mu in enumerate(gmmModel.means_):
                data[rank+1][j][(k+1)*3] = gmmModel.means_[k][0]
                data[rank+1][j][((k+1)*3)+1] = gmmModel.covars_[k][0]
                data[rank+1][j][((k+1)*3)+2] = gmmModel.weights_[k]
            j+=1 #increment sample index

    subOutfileName = outFile + str(rank)
    subOutfile = open(subOutfileName, 'w')
    writeString = str()

    
    for l, param in enumerate(data[rank+1]):
        subOutfile.write(str(param[0])+'\t'+str(param[1])+'\t'+str(param[2])+'\t'+str(param[3])+'\t'+str(param[4])+'\t'+str(param[5])+'\t'+str(param[6])+'\t'+str(param[7])+'\t'+str(param[8])+'\t'+str(param[9])+'\t'+str(param[10])+'\t'+str(param[11])+'\t'+str(param[12])+'\t'+str(param[13])+'\t'+str(param[14])+'\t'+str(param[15])+'\t'+str(param[16])+'\t'+str(param[17])+'\t'+str(param[18])+'\t'+str(param[19])+'\t'+str(param[20])+'\t'+str(param[21])+'\t'+str(param[22])+'\t'+str(param[23])+'\t'+str(param[24])+'\t'+str(param[25])+'\t'+str(param[26])+'\t'+str(param[27])+'\t'+str(param[28])+'\t'+str(param[29])+'\t'+str(param[30])+'\t'+str(param[31])+'\t'+str(param[32])+'\t'+'\n')

    subOutfile.close()
if rank == 1:
    #initialize local sample index
    j = 0

    for i in range(maxInitialCount):
        for n in range(maxPCRcycles):
            gmmModel = amp.GMMTest(i+1, n+1, 0.2, 10000, maxGaussNum)
            data[rank+1][j][2] = 0.2
            data[rank+1][j][0] = i+1
            data[rank+1][j][1] = n+1
            for k, mu in enumerate(gmmModel.means_):
                data[rank+1][j][(k+1)*3] = gmmModel.means_[k][0]
                data[rank+1][j][((k+1)*3)+1] = gmmModel.covars_[k][0]
                data[rank+1][j][((k+1)*3)+2] = gmmModel.weights_[k]
            j+=1 #increment sample index

    subOutfileName = outFile + str(rank)
    subOutfile = open(subOutfileName, 'w')
    writeString = str()

    
    
    for l, param in enumerate(data[rank+1]):
        subOutfile.write(str(param[0])+'\t'+str(param[1])+'\t'+str(param[2])+'\t'+str(param[3])+'\t'+str(param[4])+'\t'+str(param[5])+'\t'+str(param[6])+'\t'+str(param[7])+'\t'+str(param[8])+'\t'+str(param[9])+'\t'+str(param[10])+'\t'+str(param[11])+'\t'+str(param[12])+'\t'+str(param[13])+'\t'+str(param[14])+'\t'+str(param[15])+'\t'+str(param[16])+'\t'+str(param[17])+'\t'+str(param[18])+'\t'+str(param[19])+'\t'+str(param[20])+'\t'+str(param[21])+'\t'+str(param[22])+'\t'+str(param[23])+'\t'+str(param[24])+'\t'+str(param[25])+'\t'+str(param[26])+'\t'+str(param[27])+'\t'+str(param[28])+'\t'+str(param[29])+'\t'+str(param[30])+'\t'+str(param[31])+'\t'+str(param[32])+'\t'+'\n')

    subOutfile.close()
if rank == 2:
    #initialize local sample index
    j = 0

    for i in range(maxInitialCount):
        for n in range(maxPCRcycles):
            gmmModel = amp.GMMTest(i+1, n+1, 0.3, 10000, maxGaussNum)
            data[rank+1][j][2] = 0.3
            data[rank+1][j][0] = i+1
            data[rank+1][j][1] = n+1
            for k, mu in enumerate(gmmModel.means_):
                data[rank+1][j][(k+1)*3] = gmmModel.means_[k][0]
                data[rank+1][j][((k+1)*3)+1] = gmmModel.covars_[k][0]
                data[rank+1][j][((k+1)*3)+2] = gmmModel.weights_[k]
            j+=1 #increment sample index

    subOutfileName = outFile + str(rank)
    subOutfile = open(subOutfileName, 'w')
    writeString = str()

    
    for l, param in enumerate(data[rank+1]):
        subOutfile.write(str(param[0])+'\t'+str(param[1])+'\t'+str(param[2])+'\t'+str(param[3])+'\t'+str(param[4])+'\t'+str(param[5])+'\t'+str(param[6])+'\t'+str(param[7])+'\t'+str(param[8])+'\t'+str(param[9])+'\t'+str(param[10])+'\t'+str(param[11])+'\t'+str(param[12])+'\t'+str(param[13])+'\t'+str(param[14])+'\t'+str(param[15])+'\t'+str(param[16])+'\t'+str(param[17])+'\t'+str(param[18])+'\t'+str(param[19])+'\t'+str(param[20])+'\t'+str(param[21])+'\t'+str(param[22])+'\t'+str(param[23])+'\t'+str(param[24])+'\t'+str(param[25])+'\t'+str(param[26])+'\t'+str(param[27])+'\t'+str(param[28])+'\t'+str(param[29])+'\t'+str(param[30])+'\t'+str(param[31])+'\t'+str(param[32])+'\t'+'\n')

    subOutfile.close()
if rank == 3:
    #initialize local sample index
    j = 0

    for i in range(maxInitialCount):
        for n in range(maxPCRcycles):
            gmmModel = amp.GMMTest(i+1, n+1, 0.4, 10000, maxGaussNum)
            data[rank+1][j][2] = 0.4
            data[rank+1][j][0] = i+1
            data[rank+1][j][1] = n+1
            for k, mu in enumerate(gmmModel.means_):
                data[rank+1][j][(k+1)*3] = gmmModel.means_[k][0]
                data[rank+1][j][((k+1)*3)+1] = gmmModel.covars_[k][0]
                data[rank+1][j][((k+1)*3)+2] = gmmModel.weights_[k]
            j+=1 #increment sample index

    subOutfileName = outFile + str(rank)
    subOutfile = open(subOutfileName, 'w')
    writeString = str()

    for l, param in enumerate(data[rank+1]):
        subOutfile.write(str(param[0])+'\t'+str(param[1])+'\t'+str(param[2])+'\t'+str(param[3])+'\t'+str(param[4])+'\t'+str(param[5])+'\t'+str(param[6])+'\t'+str(param[7])+'\t'+str(param[8])+'\t'+str(param[9])+'\t'+str(param[10])+'\t'+str(param[11])+'\t'+str(param[12])+'\t'+str(param[13])+'\t'+str(param[14])+'\t'+str(param[15])+'\t'+str(param[16])+'\t'+str(param[17])+'\t'+str(param[18])+'\t'+str(param[19])+'\t'+str(param[20])+'\t'+str(param[21])+'\t'+str(param[22])+'\t'+str(param[23])+'\t'+str(param[24])+'\t'+str(param[25])+'\t'+str(param[26])+'\t'+str(param[27])+'\t'+str(param[28])+'\t'+str(param[29])+'\t'+str(param[30])+'\t'+str(param[31])+'\t'+str(param[32])+'\t'+'\n')

    subOutfile.close()
if rank == 4:
    #initialize local sample index
    j = 0

    for i in range(maxInitialCount):
        for n in range(maxPCRcycles):
            gmmModel = amp.GMMTest(i+1, n+1, 0.5, 10000, maxGaussNum)
            data[rank+1][j][2] = 0.5
            data[rank+1][j][0] = i+1
            data[rank+1][j][1] = n+1
            for k, mu in enumerate(gmmModel.means_):
                data[rank+1][j][(k+1)*3] = gmmModel.means_[k][0]
                data[rank+1][j][((k+1)*3)+1] = gmmModel.covars_[k][0]
                data[rank+1][j][((k+1)*3)+2] = gmmModel.weights_[k]
            j+=1 #increment sample index

    subOutfileName = outFile + str(rank)
    subOutfile = open(subOutfileName, 'w')
    writeString = str()

    for l, param in enumerate(data[rank+1]):
        subOutfile.write(str(param[0])+'\t'+str(param[1])+'\t'+str(param[2])+'\t'+str(param[3])+'\t'+str(param[4])+'\t'+str(param[5])+'\t'+str(param[6])+'\t'+str(param[7])+'\t'+str(param[8])+'\t'+str(param[9])+'\t'+str(param[10])+'\t'+str(param[11])+'\t'+str(param[12])+'\t'+str(param[13])+'\t'+str(param[14])+'\t'+str(param[15])+'\t'+str(param[16])+'\t'+str(param[17])+'\t'+str(param[18])+'\t'+str(param[19])+'\t'+str(param[20])+'\t'+str(param[21])+'\t'+str(param[22])+'\t'+str(param[23])+'\t'+str(param[24])+'\t'+str(param[25])+'\t'+str(param[26])+'\t'+str(param[27])+'\t'+str(param[28])+'\t'+str(param[29])+'\t'+str(param[30])+'\t'+str(param[31])+'\t'+str(param[32])+'\t'+'\n')

    subOutfile.close()
if rank == 5:
    #initialize local sample index
    j = 0

    for i in range(maxInitialCount):
        for n in range(maxPCRcycles):
            gmmModel = amp.GMMTest(i+1, n+1, 0.6, 10000, maxGaussNum)
            data[rank+1][j][2] = 0.6
            data[rank+1][j][0] = i+1
            data[rank+1][j][1] = n+1
            for k, mu in enumerate(gmmModel.means_):
                data[rank+1][j][(k+1)*3] = gmmModel.means_[k][0]
                data[rank+1][j][((k+1)*3)+1] = gmmModel.covars_[k][0]
                data[rank+1][j][((k+1)*3)+2] = gmmModel.weights_[k]
            j+=1 #increment sample index

    subOutfileName = outFile + str(rank)
    subOutfile = open(subOutfileName, 'w')
    writeString = str()


    for l, param in enumerate(data[rank+1]):
        subOutfile.write(str(param[0])+'\t'+str(param[1])+'\t'+str(param[2])+'\t'+str(param[3])+'\t'+str(param[4])+'\t'+str(param[5])+'\t'+str(param[6])+'\t'+str(param[7])+'\t'+str(param[8])+'\t'+str(param[9])+'\t'+str(param[10])+'\t'+str(param[11])+'\t'+str(param[12])+'\t'+str(param[13])+'\t'+str(param[14])+'\t'+str(param[15])+'\t'+str(param[16])+'\t'+str(param[17])+'\t'+str(param[18])+'\t'+str(param[19])+'\t'+str(param[20])+'\t'+str(param[21])+'\t'+str(param[22])+'\t'+str(param[23])+'\t'+str(param[24])+'\t'+str(param[25])+'\t'+str(param[26])+'\t'+str(param[27])+'\t'+str(param[28])+'\t'+str(param[29])+'\t'+str(param[30])+'\t'+str(param[31])+'\t'+str(param[32])+'\t'+'\n')

    subOutfile.close()
if rank == 6:
    #initialize local sample index
    j = 0

    for i in range(maxInitialCount):
        for n in range(maxPCRcycles):
            gmmModel = amp.GMMTest(i+1, n+1, 0.7, 10000, maxGaussNum)
            data[rank+1][j][2] = 0.7
            data[rank+1][j][0] = i+1
            data[rank+1][j][1] = n+1
            for k, mu in enumerate(gmmModel.means_):
                data[rank+1][j][(k+1)*3] = gmmModel.means_[i][k]
                data[rank+1][j][((k+1)*3)+1] = gmmModel.covars_[k][0]
                data[rank+1][j][((k+1)*3)+2] = gmmModel.weights_[k]
            j+=1 #increment sample index

    subOutfileName = outFile + str(rank)
    subOutfile = open(subOutfileName, 'w')
    writeString = str()


    for l, param in enumerate(data[rank+1]):
        subOutfile.write(str(param[0])+'\t'+str(param[1])+'\t'+str(param[2])+'\t'+str(param[3])+'\t'+str(param[4])+'\t'+str(param[5])+'\t'+str(param[6])+'\t'+str(param[7])+'\t'+str(param[8])+'\t'+str(param[9])+'\t'+str(param[10])+'\t'+str(param[11])+'\t'+str(param[12])+'\t'+str(param[13])+'\t'+str(param[14])+'\t'+str(param[15])+'\t'+str(param[16])+'\t'+str(param[17])+'\t'+str(param[18])+'\t'+str(param[19])+'\t'+str(param[20])+'\t'+str(param[21])+'\t'+str(param[22])+'\t'+str(param[23])+'\t'+str(param[24])+'\t'+str(param[25])+'\t'+str(param[26])+'\t'+str(param[27])+'\t'+str(param[28])+'\t'+str(param[29])+'\t'+str(param[30])+'\t'+str(param[31])+'\t'+str(param[32])+'\t'+'\n')

    subOutfile.close()
if rank == 7:
    #initialize local sample index
    j = 0

    for i in range(maxInitialCount):
        for n in range(maxPCRcycles):
            gmmModel = amp.GMMTest(i+1, n+1, 0.8, 10000, maxGaussNum)
            data[rank+1][j][2] = 0.8
            data[rank+1][j][0] = i+1
            data[rank+1][j][1] = n+1
            for k, mu in enumerate(gmmModel.means_):
                data[rank+1][j][(k+1)*3] = gmmModel.means_[k][0]
                data[rank+1][j][((k+1)*3)+1] = gmmModel.covars_[k][0]
                data[rank+1][j][((k+1)*3)+2] = gmmModel.weights_[k]
            j+=1 #increment sample index

    subOutfileName = outFile + str(rank)
    subOutfile = open(subOutfileName, 'w')
    writeString = str()


    for l, param in enumerate(data[rank+1]):
        subOutfile.write(str(param[0])+'\t'+str(param[1])+'\t'+str(param[2])+'\t'+str(param[3])+'\t'+str(param[4])+'\t'+str(param[5])+'\t'+str(param[6])+'\t'+str(param[7])+'\t'+str(param[8])+'\t'+str(param[9])+'\t'+str(param[10])+'\t'+str(param[11])+'\t'+str(param[12])+'\t'+str(param[13])+'\t'+str(param[14])+'\t'+str(param[15])+'\t'+str(param[16])+'\t'+str(param[17])+'\t'+str(param[18])+'\t'+str(param[19])+'\t'+str(param[20])+'\t'+str(param[21])+'\t'+str(param[22])+'\t'+str(param[23])+'\t'+str(param[24])+'\t'+str(param[25])+'\t'+str(param[26])+'\t'+str(param[27])+'\t'+str(param[28])+'\t'+str(param[29])+'\t'+str(param[30])+'\t'+str(param[31])+'\t'+str(param[32])+'\t'+'\n')

    subOutfile.close()
if rank == 8:
    #initialize local sample index
    j = 0

    for i in range(maxInitialCount):
        for n in range(maxPCRcycles):
            gmmModel = amp.GMMTest(i+1, n+1, 0.9, 10000, maxGaussNum)
            data[rank+1][j][2] = 0.9
            data[rank+1][j][0] = i+1
            data[rank+1][j][1] = n+1
            for k, mu in enumerate(gmmModel.means_):
                data[rank+1][j][(k+1)*3] = gmmModel.means_[k][0]
                data[rank+1][j][((k+1)*3)+1] = gmmModel.covars_[k][0]
                data[rank+1][j][((k+1)*3)+2] = gmmModel.weights_[k]
            j+=1 #increment sample index

    subOutfileName = outFile + str(rank)
    subOutfile = open(subOutfileName, 'w')
    writeString = str()

    for l, param in enumerate(data[rank+1]):
        subOutfile.write(str(param[0])+'\t'+str(param[1])+'\t'+str(param[2])+'\t'+str(param[3])+'\t'+str(param[4])+'\t'+str(param[5])+'\t'+str(param[6])+'\t'+str(param[7])+'\t'+str(param[8])+'\t'+str(param[9])+'\t'+str(param[10])+'\t'+str(param[11])+'\t'+str(param[12])+'\t'+str(param[13])+'\t'+str(param[14])+'\t'+str(param[15])+'\t'+str(param[16])+'\t'+str(param[17])+'\t'+str(param[18])+'\t'+str(param[19])+'\t'+str(param[20])+'\t'+str(param[21])+'\t'+str(param[22])+'\t'+str(param[23])+'\t'+str(param[24])+'\t'+str(param[25])+'\t'+str(param[26])+'\t'+str(param[27])+'\t'+str(param[28])+'\t'+str(param[29])+'\t'+str(param[30])+'\t'+str(param[31])+'\t'+str(param[32])+'\t'+'\n')

    subOutfile.close()
if rank == 9:
    #initialize local sample index
    j = 0

    for i in range(maxInitialCount):
        for n in range(maxPCRcycles):
            gmmModel = amp.GMMTest(i+1, n+1, 1.0, 10000, maxGaussNum)
            data[rank+1][j][2] = 1.0
            data[rank+1][j][0] = i+1
            data[rank+1][j][1] = n+1
            for k, mu in enumerate(gmmModel.means_):
                data[rank+1][j][(k+1)*3] = gmmModel.means_[k][0]
                data[rank+1][j][((k+1)*3)+1] = gmmModel.covars_[k][0]
                data[rank+1][j][((k+1)*3)+2] = gmmModel.weights_[k]
            j+=1 #increment sample index

    subOutfileName = outFile + str(rank)
    subOutfile = open(subOutfileName, 'w')
    writeString = str()


    for l, param in enumerate(data[rank+1]):
        subOutfile.write(str(param[0])+'\t'+str(param[1])+'\t'+str(param[2])+'\t'+str(param[3])+'\t'+str(param[4])+'\t'+str(param[5])+'\t'+str(param[6])+'\t'+str(param[7])+'\t'+str(param[8])+'\t'+str(param[9])+'\t'+str(param[10])+'\t'+str(param[11])+'\t'+str(param[12])+'\t'+str(param[13])+'\t'+str(param[14])+'\t'+str(param[15])+'\t'+str(param[16])+'\t'+str(param[17])+'\t'+str(param[18])+'\t'+str(param[19])+'\t'+str(param[20])+'\t'+str(param[21])+'\t'+str(param[22])+'\t'+str(param[23])+'\t'+str(param[24])+'\t'+str(param[25])+'\t'+str(param[26])+'\t'+str(param[27])+'\t'+str(param[28])+'\t'+str(param[29])+'\t'+str(param[30])+'\t'+str(param[31])+'\t'+str(param[32])+'\t'+'\n')

    subOutfile.close()
