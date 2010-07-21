#!/usr/bin/python

from Numeric import *
from string import *
from sys import *

def ReadArray(file_name):
    """ Read an array from the data file.
    
    The function 'readArray' from 'Scientific.IO.ArrayIO' doesn't
    handle hash marks correctly. """
    data = []
    file = open(file_name, "r")
    for line in file.readlines():
        # if the first non-blank character is the hash mark, skip the line
        # parse the line
        words = split(line)
        record = []
        for w in words:
            if w[0][0] == "#": break
            record.append(atof(w))
        if len(record) > 1: data.append(record)
        elif len(record) == 1: data.append(record[0])
    file.close()
    return array(data)

if len(argv) < 2:
    print "USAGE: gp2xmgrace <analysis file>"
    exit(1)

data = ReadArray(argv[1])
for i in range(len(data)):
    print data[i][0], data[i][1], data[i][3]-data[i][1], data[i][1]-data[i][2]
