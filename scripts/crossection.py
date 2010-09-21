#!/usr/bin/env python

from sys import argv

if (len(argv) < 2):
	print "Usage", argv[0], "<design_name>"
	quit()

with open(argv[1]+'.E.dat', 'r') as f:
	lines = f.readlines()

tupline = [tuple(line.split()) for line in lines if line != "\n"]
ordered = sorted(tupline, key=lambda line: float(line[2]), reverse=True)

def uniq(seq, idfun=None): 
	if idfun is None:
		def idfun(x): return x
	seen = {}
	result = []
	for item in seq:
		marker = idfun(item)
		if marker in seen: continue
		seen[marker] = 1
		result.append(item)
	return result

unique = uniq(ordered, lambda x: x[1])

with open(argv[1]+'.max.dat', 'w') as f:
	for line in sorted(unique, key=lambda uq: float(uq[1])):
		_, wavelength, intensity = line
		f.write("%5.1f %s \n" % (float(wavelength), intensity))


