#!/usr/bin/python

# first argument is enrichment ratio threshold
# second argument is confounding threshold
# third argument is path to data
# fourth argument is the number of times a chemical has to enrich a sequence
# ie, how many different pools has a chemical that enriches a sequence
# fifth argument gives how many pools a sequence in a pool has to beat
import os
import sys

#sys.argv = ["","", "", "", ""]
#sys.argv[1] = 5
#sys.argv[2] = 0
#sys.argv[3] = "./data"
#sys.argv[4] = 2

enrichmentThreshold = float(sys.argv[1]) #ratio required between new pool and original to be considered enriched
confoundingThreshold = float(sys.argv[2]) #difference between chemical and _notchemical to be considered nonconfounding
path = sys.argv[3] # path to the data and output log files
poolThreshold = float(sys.argv[4]) # in how many pools must a sequence be enriched in
majorityThreshold = float(sys.argv[5]) #how many pools must a sequence from a particular pool beat

# I want to process the .his values which are just raw counts of every sequence that appears
# I want to process the .val files which are derviced from the .proc files
# formatting for .val files [.proc fileName] [enrichedCount] [originalCount] [enrichmentRatio]

dirList=os.listdir(path)
#read in configuration file and store each unique chemical
with open(path+"/chempools2.config") as f:
     configFile = f.readlines()
_chemDict = dict() # stores the chemicals that are particular to each pool
_seqPoolDict = dict() # stores the pools that the sequences appear in. key is sequence, value is a list of pools
currentPool = "";
tempChemList = list();
for line in configFile:
     if line.startswith("#"):
          if len(currentPool)>0:
               _chemDict[currentPool] = set(currentChemicals);
          currentChemicals = list();
          currentPool = line.split(" ")[1].strip("\n").upper() #all pools are upper case
     elif not line.startswith("#") and len(line) > 1:
          currentChem=line.rstrip("\n").rstrip().lower().replace(" ","_") #chemicals are lower case, and have no spaces
          currentChemicals.append(currentChem)
          tempChemList.append(currentChem)

_chemDict[currentPool] = set(currentChemicals); #add the last entry
_chemList = set(tempChemList) #list of unique chemicals needed for initialzing entries in _chemHash
poolKeys = sorted(_chemDict.keys())

# store conflicting pools
# if two pools don't share at least one chemical, then eliminate the second pool
#_chemDict stores the chemicals that are particular to each pool
_conflictHash = dict()
for i in range(len(poolKeys)):
   _conflictHash[poolKeys[i]] = []
   for j in range(len(poolKeys)):
       if len(_chemDict[poolKeys[i]] &  _chemDict[poolKeys[j]]) < 1: # two pools must share at least one chemical
          # store this as a conflict
          print poolKeys[i] + " conflicts with "+ poolKeys[j]
          _conflictHash[poolKeys[i]].append(poolKeys[j])
for key in _conflictHash.keys():
   _conflictHash[key] = set(_conflictHash[key])


_populationHash = dict() #key is the pool, value is the total count for the pool

# read in the his files specified by the chempool.config

_ratioHash = dict() #a dictionary that will hold the enrichment ratios for all of the pools (excepting for the last one, the prephage selection
                    #key is the pool name as indicated in the config file, value is a dictionary containing the enrichment values
# calculate the total populations

for i in range(len(poolKeys)):
     totalCount = 0
     with open("./data/"+poolKeys[i]+".his") as f:
          hisFile = f.readlines()
     for line in hisFile:
          tokens = line.rstrip("\n").split("\t")
          totalCount = totalCount + float(tokens[1])
     #print("total population for "+poolKeys[i]+": "+ str(totalCount))     
     _populationHash[poolKeys[i]] = totalCount

# calculate ratios for each pool
_ratioHash = dict() # dictionary of dictionaries; key is the pool, value is dictionary containing sequence ratio pairs
for i in range(len(poolKeys)):
   _ratioHash[poolKeys[i]] = dict();
   with open("./data/"+poolKeys[i]+".his") as f:
      hisFile = f.readlines()
   for line in hisFile:
      tokens = line.rstrip("\n").split("\t")
      sequence = tokens[0]
      seqCount = float(tokens[1])
      _ratioHash[poolKeys[i]][sequence] = seqCount/_populationHash[poolKeys[i]]
      #print("seqCount: "+str(seqCount)+" total: "+str(_populationHash[poolKeys[i]]) +" ratio: "+ str(seqCount/_populationHash[poolKeys[i]]))     
droppedSeqLog = open("droppedSequencesV2.txt", "w") # log dropped sequences
# compare each pool to all of the other pools
_resultHash = dict() # key is sequence value is set of pools that it enriches in
for i in range(len(poolKeys)):
   for sequence in _ratioHash[poolKeys[i]].keys():   
      beatenCount = 0 # how many pools has this sequence in pool[i] beaten?
      for j in range(len(poolKeys)):
         if not i==j: # don't compare a sequence to itself 
            if poolKeys[i] not in _conflictHash[poolKeys[j]]: # pools conflict so don't look at them
               if sequence in _ratioHash[poolKeys[j]].keys():
                  if _ratioHash[poolKeys[i]][sequence] >= enrichmentThreshold*_ratioHash[poolKeys[j]][sequence]: # if the ratio is greater than the enrichment threshold
                     beatenCount = beatenCount + 1
               else:
                  droppedSeqLog.write(tokens[0]+"\n")
      if beatenCount >= majorityThreshold:
         if sequence in _resultHash.keys():
            _resultHash[sequence] = _resultHash[sequence] | set(poolKeys[i]) # union is the only way I know to add to a set
         else:
            _resultHash[sequence] = set(poolKeys[i])
         #print "we've got a hit; beatenCount: "+str(beatenCount)+" majorityThreshold: "+str(majorityThreshold)
# print out results

print("number of hits: "+str(len(_resultHash.keys())))










