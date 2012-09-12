#!/usr/bin/python
import os
path="./data"  # insert the path to the directory of interest here

# I want to process the .his values which are just raw counts of every sequence that appears
# I want to process the .val files which are derviced from the .proc files
# formatting for .val files [.proc fileName] [enrichedCount] [originalCount] [enrichmentRatio]

dirList=os.listdir(path)
#read in configuration file and store each unique chemical
with open("./data/chempools.config") as f:
     configFile = f.readlines()
_chemDict = dict();
currentPool = "";
tempChemList = list();
for line in configFile:
     if line.startswith("#"):
          if len(currentPool)>0:
               _chemDict[currentPool] = set(currentChemicals);
          currentChemicals = list();
          currentPool = line.split(" ")[1].strip("\n")
          #print "current pool is: " + currentPool
     elif not line.startswith("#") and len(line) > 1:
          #print "adding chemical: "+line.rstrip("\n");          
          currentChemicals.append(line.rstrip("\n"))
          tempChemList.append(line.rstrip("\n"))
_chemList = set(tempChemList) #list of unique chemicals needed for initialzing entries in _chemHash


#initialize chemHash
_chemHash = dict(); 
#this dictionary will contain all of they sequences and the chemicals that enrich them
#key will be the sequence as a string, value will be a dictionary containing each of the chemicals
#and the number of times each chemical has enriched the given sequence

# this function will initialize the entry for a sequence in 
def addSeqEntry(sequence,chemical):
     if sequence in _chemHash:
          # if key already exists, update it
          _chemHash[sequence][chemical] = _chemHash[sequence][chemical]+1
     else:
          # if key doesn't exist, initialize an empty dictionary for the value, and then incrememnt the chemcial
          newDict = dict();
          for chem in _chemList:
                        newDict[chem]= 0
                        newDict["not_"+chem]=0; #also initialize an entry corresponding to the lack of a particuar chemical;
                                                #this will be used to account for confounding data
          newDict[chemical] = 1;
          _chemHash[sequence] = newDict;
          


#testArea
addSeqEntry("moo","Urea")
addSeqEntry("cow","Urea")
addSeqEntry("moo","Urea")


print(_chemHash["cow"]["Urea"])























     



files = [];
for fname in dirList:
     files.append(fname)
#for fileName in files:
     #print(fileName)
