#!/usr/bin/python
# first argument is path to first file
# second argument is path to second file
# prints sequences that appear in both files
# and the ratio of ratios of both sequenecs second/first
import os
import sys
sys.argv=["","",""]
path = "data/results/results552/"
sys.argv[1] = path + "round1.csv"
sys.argv[2] = path + "round2.csv"

fileDict1 = dict() # key is sequence, value is rest of contents of csv file line
fileDict2 = dict() # key is sequence, value is rest of contents of csv file line

with open(sys.argv[1]) as f:
     file1 = f.readlines()
for line in file1:
    tokens = line.rstrip("\n").split(",")
    sequence = tokens[0]
    contents = tokens[1:]
    #print contents
    #print sequence
    fileDict1[sequence] = contents

with open(sys.argv[2]) as f:
     file2 = f.readlines()
for line in file2:
    tokens = line.rstrip("\n").split(",")
    sequence = tokens[0]
    contents = tokens[1:]
    #print contents
    #print sequence
    fileDict2[sequence] = contents

count = 0;
for sequence in fileDict1.keys():
    if sequence in fileDict2.keys():
        result1 = sequence
        result2 = ""
        for item in fileDict1[sequence]:
            result1=result1+","+item
        for item in fileDict2[sequence]:
            result2=result2+","+item                          
        print result1
        print result2
        count = count + 1
print "Number of matches: "+str(count)
        
