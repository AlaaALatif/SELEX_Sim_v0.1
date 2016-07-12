import sys
import time
import random
import math
import numpy as np
import Predictor

Mlearn = Predictor.Predictors()

inputDimns = 3
outputDimns = 12
trainingFile = "trainingSet"
testFile = "testSet"

#print("HELLOOOOO?")

print("input dimensions = "+str(inputDimns))
print("output dimensions = "+str(outputDimns))



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

    for j in range(outputDimns):
        yTrain[sample][j] = float(line.split()[j+inputDimns])
    line = trainingSet.readline()
trainingSet.close()


#print xTrain[5][2]
#print("yTrain[5][10] = "+str(yTrain[5][10]))
# Fetch test set values
testSet = open(testFile, 'r')
line = testSet.readline()

for sample in range(testSetSize):
    for i in range(inputDimns):
        xTest[sample][i] = float(line.split()[i])

    for j in range(outputDimns):
        yTest[sample][j] = float(line.split()[j+inputDimns])
    line = testSet.readline()
testSet.close()

#print xTest[5][1]
#print("yTest[5][9] = "+str(yTest[5][9]))

# Run a prediction experiment using decision tree regression
dtrModel = Mlearn.decisionTreeRegression(xTrain, yTrain, xTest, yTest)
# Run a prediction experiment using support vector regression
svrModel = Mlearn.supportVectorRegression(xTrain, yTrain, xTest, yTest)
# Run a prediction experiment using linear regression
lrModel = Mlearn.linearRegression(xTrain, yTrain, xTest, yTest)
# Run a prediction experiment using ridge regression
rrModel = Mlearn.ridgeRegression(xTrain, yTrain, xTest, yTest)
# Run a prediction experiment using bayesian ridge regression
brrModel = Mlearn.bayesianRidgeRegression(xTrain, yTrain, xTest, yTest)

# See 'Predictors.py' program for more info on experiment design


