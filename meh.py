#!/usr/bin/python
import os
import sys

with open("./data/"+sys.argv[1]) as f:
   fileOne = f.readlines()
results = dict()
for line in fileOne:
   tokens = line.rstrip("\n").split("\t")
   print line

