import sys
import time
import random
import math
import numpy as np
import Predictor

Mlearn = Predictor.Predictors()

inputDimns = int(sys.argv[1])
outputDimns = int(sys.argv[2])
trainingFile = str(sys.argv[3])
testFile = str(sys.argv[4])

with open(trainingFile) as f: #open initial library
            for trainingSetSize, l in enumerate(f): #calculate total num of seqs
                pass
            trainingSetSize+=1
f.close()

with open(testFile) as t: #open initial library
            for testSetSize, l in enumerate(t): #calculate total num of seqs
                pass
            testSetSize+=1
f.close()



sampleNum = 100

# Initialize matrices for inputs and outputs, for both training and test set
xTrain = np.zeros((trainingSetSize, inputDimns))
yTrain = np.zeros((trainingSetSize, outputDimns))
xTest = np.zeros((testSetSize, inputDimns))
yTest = np.zeros((testSetSize, outputDimns))

# Fetch training set values
trainingSet = open(trainingFile, 'r')
line = trainingSet.readline()
for sample in range(trainingSetSize):
    for i in range(inputDimns):
        xTrain[sample][i] = float(line.split()[i])

    for j in range(inputDimns, outputDimns):
        yTrain[sample][j] = float(line.split()[j])
    line = trainingSet.readline()
trainingSet.close()

# Fetch test set values
testSet = open(testFile, 'r')
line = testSet.readline()
for sample in range(testSetSize):
    for i in range(inputDimns):
        xTest[sample][i] = float(line.split()[i])

    for j in range(outputDimns):
        yTest[sample][j] = float(line.split()[j])
    line = testSet.readline()
testSet.close()

svrModel = Mlearn.supportVectorRegression(xTrain, yTrain, xTest, yTest)
# Run a prediction experiment using linear regression
lrModel = Mlearn.linearRegression(xTrain, yTrain, xTest, yTest)
# Run a prediction experiment using ridge regression
rrModel = Mlearn.ridgeRegression(xTrain, yTrain, xTest, yTest)
# Run a prediction experiment using bayesian ridge regression
brrModel = Mlearn.bayesianRidgeRegression(xTrain, yTrain, xTest, yTest)

# See 'Predictors.py' program for more info on experiment design


