#!/usr/bin/python
# combines two hist files
# first argument is the first hist file
# second argument is the second hist file
# third arguement is the output file name

import os
import sys
sys.argv=["","",""]
sys.argv[1]="data/A.his"
sys.argv[2]="data/B.his"

seqDict = dict() # store sequences and counts
# read in first file
with open("./"+sys.argv[1]) as f:
     fileOne = f.readlines()
for line in fileOne:
    tokens = line.rstrip("\n").split("\t")
    seqDict[tokens[0]] = int(tokens[1])
    #print "stored value for "+tokens[0]+": "+ tokens[1]

# read in second file
with open("./"+sys.argv[2]) as f:
    fileTwo = f.readlines()
for line in fileTwo:
    tokens = line.rstrip("\n").split("\t")
    if tokens[0] in seqDict:
        #print "incrementing " + tokens[0] + " by"  +tokens[1] + " from " +str(seqDict[tokens[0]])
        seqDict[tokens[0]] = seqDict[tokens[0]] + int(tokens[1])
        #print "resulting in :" + str(seqDict[tokens[0]])
    else:
        seqDict[tokens[0]] = int(tokens[1])
file = open("./combined.his", 'w+')
sortedKeys = sorted(seqDict.keys())
for key in sortedKeys:
    #print key+"\t"+str(seqDict[key])
    file.write(key+"\t")
    file.write(str(seqDict[key])+"\n")
