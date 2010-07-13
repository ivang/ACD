#!/usr/bin/python
""" Exports target GDD(lambda) from the file 'parameters.txt'.

The target dispersion is written to the file 'targetGDD.dat' in
the current directory. This file has 2 columns: the wavelength in nanometers
and the GDD in femtoseconds squared. """

import string
import sys
from Numeric import *

N = 101 # number of points to export (number of dots in the plot)
N1 = -1 # for the ANALYSIS RESOLUTION
lambda_min = 550.0
lambda_max = 1100.0
central_wavelength = 800.0
C = 300.0 # speed of light, nm/fs
target_dispersion = [0.0]

file = open("parameters.txt", "r")
for line in file.readlines():
	if line[0] == "#": continue
	words = string.split(string.upper(line))
	if "=" not in words: continue
	if len(words)>3 and words[0]=="TARGET" and words[1]=="DISPERSION" \
	       and words[2]=="=":
		if words[3][0]=="\"" or words[3][0]=="'":
			print "the target dispersion is specified \
by a file, not by its Taylor expansion"
		target_dispersion = []
		for i in range(3,len(words)):
			target_dispersion.append(string.atof(words[i]))
	elif len(words)>2 and words[0]=="RESOLUTION" and words[1]=="=":
		N = string.atoi(words[2])
	elif len(words)>3 and words[0]=="ANALYSIS" and \
		 words[1]=="RESOLUTION" and words[2]=="=":
		N1 = string.atoi(words[3])
	elif len(words)>3 and words[0] == "MINIMAL" and \
		 words[1] == "WAVELENGTH" and words[2]=="=":
		lambda_min = string.atof(words[3])
	elif len(words)>3 and words[0] == "MAXIMAL" and \
		 words[1] == "WAVELENGTH" and words[2]=="=":
		lambda_max = string.atof(words[3])
	elif len(words)>3 and words[0] == "CENTRAL" and \
		 words[1] == "WAVELENGTH" and words[2]=="=":
		central_wavelength = string.atof(words[3])
file.close()

omega0 = 2.0*pi*C / central_wavelength
omega_min = 2.0*pi*C / lambda_max
omega_max = 2.0*pi*C / lambda_min
step = (omega_max - omega_min) / (N - 1)
frequencies = arange(omega_min, omega_max + 0.5*step, step)
frequencies = frequencies[::-1]
wavelengths = 2.0*pi*C / frequencies
x = frequencies - omega0
target_GDD = zeros((len(x),), Float)
order = len(target_dispersion) + 1
for i in range(len(target_dispersion)-1,0,-1):
	target_GDD += target_dispersion[i]
	target_GDD *= x / (order-2.0)
	order -= 1
target_GDD += target_dispersion[0]

# write "targetGDD.dat"
file = open("targetGDD.dat", 'w')
file.write("# Group Delay Dispersion (fs^2) Versus Wavelength (nm)\n")
for i in range(len(wavelengths)):
	file.write("%g\t%g\n" % (wavelengths[i], target_GDD[i]))
file.close()
