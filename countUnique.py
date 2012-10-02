#!/usr/bin/python

# produces a csv file that gives ordered pairs of (numberOfSequences, numberOfUniqueSequences)
# while iterating all fastq files in a directory

# first argument is the directory where data is contained
import os
import sys
import re

path = str(sys.argv[1]) 
sequenceDict = dict() # i use this like a hashset
result = []
fileList = os.listdir(path)
numberOfSequences = 0
for file in fileList:
   if "test" in file:
      with open(path + "/" +
file) as f:
         for line in f:
            if re.match("^[actgnryYRNACTG]*$", line) and len(line)>1:
               sequence = line.rstrip("\n")
               if sequence in sequenceDict:
                  result.append(len(sequenceDict.keys()));
                  numberOfSequences = numberOfSequences +1;
               else:
                  sequenceDict[sequence] = "" #don't need a value, because I just want a hashset
                  result.append(len(sequenceDict.keys()));
                  numberOfSequences = numberOfSequences +1;
for i in xrange(len(result)):
   print str(i+1) + "," + str(result[i])



