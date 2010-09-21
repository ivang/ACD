#!/usr/bin/env python
# -*- coding: utf-8 -*-

from sys import argv, exit
from os.path import basename
from numpy import arange

__step = 5 
A = [1.872, 6.845, 4.833]  # coefficients for HfO2
			   # TODO: "hard coding" such parameters is
			   # considered not a good practice/evil.

def n_cauchy(coeff, x):
    """ Calculates the refractive index of a medium, using 
    the Cauchy's equation, given the wavelength x and the coefficients, 
    specific to that particular material. """
    return coeff[0] + coeff[1]/(x**2) + coeff[2]/(x**4)

if __name__ == '__main__':

    if len(argv) < 4:
	print "Usage", basename(argv[0]), "<filename> <lambda_min> <lambda_max>"
	exit(1)	

    _, filename, lambda_min, lambda_max = argv
    
    with open(filename, 'w') as f:
	for x in arange(float(lambda_min), float(lambda_max)+0.1, __step):
	    f.write('%6.2f %.8f %.8f\n' % (x, n_cauchy(A, x), float(0)))

