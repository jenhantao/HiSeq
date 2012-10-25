#!/usr/bin/python
import os
import sys

with open("./data/"+sys.argv[1]) as f:
   fileOne = f.readlines()
results = dict()
for line in fileOne:
   tokens = line.rstrip("\n").split("\t")
   results[tokens[0]] = int(tokens[1])
with open("./data/"+sys.argv[2]) as f:
   fileTwo = f.readlines()
for line in fileTwo:
   tokens = line.rstrip("\n").split("\t")
   if tokens[0] in results.keys():
      #print str(results[tokens[0]]) +" + "+str(tokens[1]) +" = " +str(results[tokens[0]]+int(tokens[1]))
      results[tokens[0]] = results[tokens[0]] + int(tokens[1])
   else:
      results[tokens[0]] = int(tokens[1])
print "about to print"
for key in results.keys():
   print key+"\t"+results[key]

