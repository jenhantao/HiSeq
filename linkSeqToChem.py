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
_chemDict = dict(); # stores the chemicals that are particular to each pool
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
del _chemDict["Enrichment"] #Enrichment is a variable, not a pool
_chemDict[currentPool] = set(currentChemicals); #add the last entry
_chemList = set(tempChemList) #list of unique chemicals needed for initialzing entries in _chemHash
enrichmentThreshold = float(configFile[0].split(" ")[3]) #ratio required between new pool and original to be considered enriched


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
          

# read in the his files specified by the chempool.config

_ratioHash = dict() #a dictionary that will hold the enrichment ratios for all of the pools (excepting for the last one, the prephage selection
                    #key is the pool name as indicated in the config file, value is a dictionary containing the enrichment values
poolKeys = sorted(_chemDict.keys())
# last pool is presumed to be the pre-phage selection
prephageHistValues = dict(); #key is sequence, value is count
with open("./data/"+poolKeys[-1]+".his") as f:
     hisFile = f.readlines()
for line in hisFile:
     if len(line)>0:
          tokens = line.rstrip("\n").split("\t")
          #print "seq: "+tokens[0]+" count: "+tokens[1]
          prephageHistValues[tokens[0]]=tokens[1]
for i in range(len(poolKeys)-1):
     #print "reading in: "+poolKeys[i]
     with open("./data/"+poolKeys[i]+".his") as f:
          hisFile = f.readlines()
     currentRatioDict = dict()
     for line in hisFile:
          tokens = line.rstrip("\n").split("\t")
          if tokens[0] in prephageHistValues:
               ratio = float(tokens[1])/float(prephageHistValues[tokens[0]])
               if ratio >= enrichmentThreshold:
                    currentRatioDict[tokens[0]] = ratio
                    #print tokens[0]+"|"+str(currentRatioDict[tokens[0]])
     _ratioHash[poolKeys[i]]=currentRatioDict

# produce the data set that we're interested in, list of sequences and the chemicals that enrich them
for key in _ratioHash:
     print key
     
     
     

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
