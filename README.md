#Created by Alaa Abdel Latif on 25/08/2016. 
#Copyright © 2016 Alaa Abdel Latif. All rights reserved.

This program was developed as part of a research project at University of Cambridge under supervision of Dr. Lucy Colwell. 

This is a Numerical Simulation program of the Systematic Evolution of Ligands through eXponential enrichment (NSELEX)

This program allows the user to simulate the dynamic changes in sequence populations during SELEX experiments through probabilistic modelling. The simulation can be carried out using either options that take into account primary and/or secondary structures through energy-based affinity measures. 


INSTALLATION

Extract using:
$tar -czvf selex_simulation.tar.gz

PACKAGE DEPENDENCY

Ensure you have the following Python (>2.7) module dependencies:

	sys
	math
	itertools
	ConfigParser

	scipy
	numpy
	matplotlib
	seaborn

If any of these modules are missing are missing they can be installed using:
$sudo pip install [modulename]
If you do not have root priviledge and cannot use sudo, then you can install the missing module(s) locally using:
$pip install --user --upgrade [modulename]
Any further issues with installation should be informed to respective system administrators (or if it's a toughy just email us).

The program also requires the ViennaRNA package to be installed. This can be downloaded from:
https://www.tbi.univie.ac.at/RNA/#download
Extract the downloaded ViennaRNA package using:
$tar -czvf ViennaRNA-[version].tar.gz
This can then be installed using:
$cd ViennaRNA-[version]/
$./configure --without-perl --with-python
$make
$sudo make install
Again, if user does not have root privileges, then a local installation can be done using:
$./configure --prefix=[LOCALPATH] --without-perl --with-python
$make
$make install

(if you get an error during make, try rerunning configure with --distable-flto)

USAGE

All of the parameters for the simulation can be specified from the 'setting.init' file. Descriptions for each parameter is given inside the file. The default values in the settings file correspond to the conditions used to report the results in the corresponding thesis.
After specifying the parameters, save the settings file and then run the simulation from the command-line using:
$python sim_.py

Please note that under the default parameters, the simulation run takes almost 4 hours on an Intel(R)Core(TM) Quad CPU Q9400 machine. Using a large scale parameter or a large number of pcr cycles can result in excessive CPU time and memory use. 

Please report any issues to aaaa3@cam.ac.uk or ljc37@cam.ac.uk
