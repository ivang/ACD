#!/usr/bin/python
""" Constructs a Bragg reflector.

The script takes two parameters from the command line: the number of
layers and the angle of incidence (in degrees). The result is printed
to the standard output. """

from sys import *
from string import *
from math import *

N = atoi(argv[1]) # the number of layers
alpha = atof(argv[2]) * pi/180.0 # the angle of incidence

material1 = "SiO2"
material2 = "Nb2O5"
n1 = 1.46
n2 = 2.22
lambda0 = 800.0
d1 = lambda0/4.0 / sqrt(1.0-(sin(alpha)/n1)**2) / n1
d2 = lambda0/4.0 / sqrt(1.0-(sin(alpha)/n2)**2) / n2

for i in range(N):
	if i % 2: print "%s\t%.1f" % (material2, d2)
	else: print"%s\t%.1f" % (material1, d1)
