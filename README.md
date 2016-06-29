# SELEX Simulations
# This program is intended to simulate Systematic Evolution of Ligands through Exponential Enrichment (SELEX) experiments. This consists of two main steps: selection and amplification.

# The code is divided into separate classes. Each class consists of different methods that perform the same operation but under different assumptions. 

# The current main code "simulationsv5.py" runs a hamming-distance based stochastic selection followed by definite PCR amplification that assumes no bias (mutation is considered). The user must specify number of aptamers, initial pool file, number of rounds, selection rate, number of pcr cycles, pcr efficiency, mutation rate and length of sequences. 

# Example usage:
# python simulationv5 1 H1_R0.sense_seq.txt 3 0.0001 10 0.75 0.001 29
# The above command-line argumant should run a simulation using 1 'optimum' aptamer sequence, an initial pool file, 3 selection rounds, a selection rate of 0.1%, 10 pcr cycles per round, pcr yield of 75%, error rate of 0.1%, and a sequence length of 29 nucleotides. 

# It has been observed that population distributions after stochastic amplification have non-trivial properties.
# Although exact solutions of expectations and variance can be computed, the distribution cannot and varies 
# under different initial conditions. 
# The program 'pcrPhaseDiagram.py' is intended to generate histograms of the population distributions under 
# different initial conditions


# TO DO:
# Edit main to allow user to run under any of the specified methods for each class
# Outputing seqs takes too long between rounds, find a fix
# Incorporate mutations into pcr amplification [DONE]
# incorporate pcr bias into pcr amplification
# Add more edit distance options for selection
