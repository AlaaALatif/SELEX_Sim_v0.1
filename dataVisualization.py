import sys
import time
import random
import math
import numpy as np
import Predictor
import matplotlib
matplotlib.use('Agg')
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt

## USER INPUT PARAMETERS
firstDatafile = str(sys.argv[1])

secondDatafile = str(sys.argv[2])

## MATRIX INITIALIZATION
firstCon = np.zeros((30, 6))
secondCon = np.zeros((30, 6))
thirdCon = np.zeros((30, 6))
fourthCon = np.zeros((30, 6))
## DATA COLLECTION
firstData = open(firstDatafile, 'r')

for j in range(30):
    line = firstData.readline()
    firstCon[j][0] = float(line.split()[0]) #initial count
    firstCon[j][2] = float(line.split()[1]) #cycle no.
    firstCon[j][1] = float(line.split()[2]) #yield
    firstCon[j][3] = float(line.split()[24]) #mean of 1st component
    firstCon[j][4] = float(line.split()[25]) #variance of 1st component
    firstCon[j][5] = float(line.split()[26]) #weight of 1st component

for j in range(30):
    line = firstData.readline()

for j in range(30):
    line = firstData.readline()
    secondCon[j][0] = float(line.split()[0]) #initial count
    secondCon[j][2] = float(line.split()[1]) #cycle no.
    secondCon[j][1] = float(line.split()[2]) #yield
    secondCon[j][3] = float(line.split()[24]) #mean of 1st component
    secondCon[j][4] = float(line.split()[25]) #variance of 1st component
    secondCon[j][5] = float(line.split()[26]) #weight of 1st component

firstData.close()
## SecondDatafile contains all data pofloat at 40% yield
secondData = open(secondDatafile, 'r')

for j in range(30):
    line = secondData.readline()
    thirdCon[j][0] = float(line.split()[0]) #initial count
    thirdCon[j][2] = float(line.split()[1]) #cycle no
    thirdCon[j][1] = float(line.split()[2]) #yield
    thirdCon[j][3] = float(line.split()[24]) #mean of 1st comp
    thirdCon[j][4] = float(line.split()[25]) #variance of 1st comp
    thirdCon[j][5] = float(line.split()[26]) #weight of 1st comp
    # initial count = 2 (SKIP)
for j in range(30):
    line = secondData.readline()
    # initial count = 3
for j in range(30):
    line = secondData.readline()
    fourthCon[j][0] = float(line.split()[0]) #initial count
    fourthCon[j][2] = float(line.split()[1]) #cycle no.
    fourthCon[j][1] = float(line.split()[2]) #yield
    fourthCon[j][3] = float(line.split()[24]) #mean of 1st comp
    fourthCon[j][4] = float(line.split()[25]) #var of 1st comp
    fourthCon[j][5] = float(line.split()[26]) #weight of 1st comp
            
secondData.close()

## DATA VISUALIZATION

# Plot 1: effect of cycle no on normalized mean of 1st GMM component

space = np.linspace(1, 30, 30)

fig1 = plt.figure()
ax = fig1.add_subplot(111)

ax.plot(space, firstCon[:,3], "o", color='r')
ax.plot(space, secondCon[:,3], "o", color='b')
ax.plot(space, thirdCon[:,3], "o", color='y')
ax.plot(space, fourthCon[:,3], "o", color='g')

ax.set_xlabel('Cycle no.')
ax.set_ylabel('normed mean')
ax.set_title('Effect of cycle number on normalised mean of 8th GMM component')

plt.grid()
plt.savefig('pcr_g8_i1_3_y4_8.pdf', format='pdf')













