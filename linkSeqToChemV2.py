#!/usr/bin/python
# first argument is path to data
# second argument is enrichment ratio threshold
# third argument gives how many pools a sequence in a pool has to beat
# fourth argument is pool threshold, in how many pools with the same chemical must a sequence enrich
import os
import sys

sys.argv = ["","", "", "", ""]
sys.argv[1] = "./data"
sys.argv[2] = 5
sys.argv[3] = 1
sys.argv[4] = 1

enrichmentThreshold = float(sys.argv[2]) #ratio required between new pool and original to be considered enriched
path = sys.argv[1] # path to the data and output log files
majorityThreshold = float(sys.argv[3]) #how many pools must a sequence from a particular pool beat
poolThreshold = float(sys.argv[4])
# I want to process the .his values which are just raw counts of every sequence that appears
# I want to process the .val files which are derviced from the .proc files
# formatting for .val files [.proc fileName] [enrichedCount] [originalCount] [enrichmentRatio]

dirList=os.listdir(path)
#read in configuration file and store each unique chemical
with open(path+"/chempools.config") as f:
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

_populationHash = dict() #key is the pool, value is the total count for the pool

# read in the his files specified by the chempool.config

_ratioHash = dict() #a dictionary that will hold the enrichment ratios for all of the pools (excepting for the last one, the prephage selection
                    #key is the pool name as indicated in the config file, value is a dictionary containing the enrichment values
# calculate the total populations

for i in range(len(poolKeys)):
     totalCount = 0
     with open(path+"/"+poolKeys[i]+".his") as f:
          hisFile = f.readlines()
     for line in hisFile:
          #print poolKeys[i]
          tokens = line.rstrip("\n").split("\t")
          totalCount = totalCount + float(tokens[1])
     _populationHash[poolKeys[i]] = totalCount

# calculate ratios for each pool
_ratioHash = dict() # dictionary of dictionaries; key is the pool, value is dictionary containing sequence ratio pairs
_averageHash = dict()
for i in range(len(poolKeys)):
   #count =0 
   _ratioHash[poolKeys[i]] = dict();
   with open(path+"/"+poolKeys[i]+".his") as f:
      hisFile = f.readlines()
   for line in hisFile:
      #print count
      #count = count + 1
      tokens = line.rstrip("\n").split("\t")
      if len(tokens)>1:
         sequence = tokens[0]
         seqCount = float(tokens[1])
         _ratioHash[poolKeys[i]][sequence] = seqCount/_populationHash[poolKeys[i]]
         if sequence not in _averageHash.keys():
            _averageHash[sequence] = 0

# compute the average ratio
for sequence in _averageHash.keys():
   count = 0;
   sum = 0;
   for i in range(len(poolKeys)):
      pool = poolKeys[i]
      if sequence in _ratioHash[pool].keys():
         count = count + 1
         sum = sum + _ratioHash[pool][sequence]
   if count > 0:
      _averageHash[sequence] = sum/count



droppedSeqLog = open("droppedSequencesV2.txt", "w") # log dropped sequences
# compare each pool to all of the other pools
_resultHash = dict() # key is sequence value is set of pools that it enriches in
for i in range(len(poolKeys)):
   print "################## "+(poolKeys[i])+" ##################"
   for sequence in _ratioHash[poolKeys[i]].keys():
      #print _ratioHash[poolKeys[i]].keys()
      beatenCount = 0 # how many pools has this sequence in pool[i] beaten?
      for j in range(len(poolKeys)):
         if not i==j: # don't compare a pool to itself
               if sequence in _ratioHash[poolKeys[j]].keys():
                  print str(_ratioHash[poolKeys[i]][sequence]) + " beats? " +str(enrichmentThreshold*_ratioHash[poolKeys[j]][sequence])
                  if _ratioHash[poolKeys[i]][sequence] >= enrichmentThreshold*_ratioHash[poolKeys[j]][sequence]: # if the ratio is greater than the enrichment threshold
                     beatenCount = beatenCount + 1
               
      if beatenCount >= majorityThreshold:
         if sequence in _resultHash.keys():
            _resultHash[sequence] = _resultHash[sequence] | set(poolKeys[i]) # union is the only way I know to add to a set
         else:
            _resultHash[sequence] = set(poolKeys[i])
      else:
         droppedSeqLog.write(tokens[0]+"\n")

# filter the result pools
for sequence in _resultHash.keys():
   currentPools = list(_resultHash[sequence])
   print currentPools
   for i in range(len(currentPools)):
      pool = currentPools[i]
      overlapCount = 1 # with how many pools does the current pool overlap with
      for j in range(len(currentPools)):
         comparedPool = currentPools[j]
         if not pool == comparedPool: #don't compare the same pool
            if len(_chemDict[pool] &  _chemDict[comparedPool]) > 0: 
               overlapCount = overlapCount + 1
      if not overlapCount >= poolThreshold: 
         # remove the pool from the results
         _resultHash[sequence] = _resultHash[sequence] - set(pool)
  
# print out results
for sequence in _resultHash.keys():
   resultString = sequence # + ", mean="+str(_averageHash[sequence])
   #chemicals that enrich are the intersection of the chemicals of the enriched pools
   chemicals = set()
   for pool in _resultHash[sequence]:
      resultString = resultString + "," + pool + "= " + str(_ratioHash[pool][sequence]/_averageHash[sequence])
   for pool in _resultHash[sequence]:
      if len(chemicals) > 0:
         chemicals = chemicals & _chemDict[pool]
      else:
         chemicals = _chemDict[pool]
   if len(chemicals) > 0:
      for chemical in chemicals:
         resultString = resultString+","+chemical
      print resultString
print("number of hits: "+str(len(_resultHash.keys())))










