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
path = sys.argv[3]
poolThreshold = float(sys.argv[4])
majorityThreshold = float(sys.argv[4])

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


#initialize chemHash
_chemHash = dict()
#this dictionary will contain all of they sequences and the chemicals that enrich them
#key will be the sequence as a string, value will be a dictionary containing each of the chemicals
#and the number of times each chemical has enriched the given sequence
_populationHash = dict() #key is the pool, value is the total count for the pool
# this function will initialize the entry for a sequence in 
def addSeqEntry(sequence,chemical):
     sequence = sequence.upper() # all sequences are upper case
     chemical = chemical.lower() # all chemicals are lower case
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
totalCount = 0;
for line in hisFile:
     if len(line)>0:
          tokens = line.rstrip("\n").split("\t")
          prephageHistValues[tokens[0]]=tokens[1]
          totalCount = totalCount + float(tokens[1])
_populationHash[poolKeys[-1]] = totalCount
#print("total population for "+poolKeys[-1]+": "+ str(totalCount))
for i in range(len(poolKeys)-1):
     totalCount = 0
     with open("./data/"+poolKeys[i]+".his") as f:
          hisFile = f.readlines()
     for line in hisFile:
          tokens = line.rstrip("\n").split("\t")
          totalCount = totalCount + float(tokens[1])
     #print("total population for "+poolKeys[i]+": "+ str(totalCount))     
     _populationHash[poolKeys[i]] = totalCount

               
for i in range(len(poolKeys)-1):
     with open("./data/"+poolKeys[i]+".his") as f:
          hisFile = f.readlines()
     currentRatioDict = dict()
     for line in hisFile:
          tokens = line.rstrip("\n").split("\t")
          if tokens[0] in prephageHistValues:
               ratio = (float(tokens[1])/_populationHash[poolKeys[i]])/(float(prephageHistValues[tokens[0]])/_populationHash[poolKeys[-1]])
               if ratio >= enrichmentThreshold:
                    currentRatioDict[tokens[0]] = ratio
     _ratioHash[poolKeys[i]]=currentRatioDict

# produce the data set that we're interested in, list of sequences and the chemicals that enrich them
for key in _ratioHash:
     positiveChemicals = set(_chemDict[key]) #chemicals that actually exist
     negativeChemicals = [] #chemicals that don't exist
     temp = set(_chemList)-positiveChemicals
     for chem in temp:
          negativeChemicals.append("not_"+chem)
     positiveChemicals = sorted(positiveChemicals)
     negativeChemicals = sorted(negativeChemicals)
     for sequence in _ratioHash[key]:
          for posChem in positiveChemicals:
               # increment the positive entries
               addSeqEntry(sequence, posChem)
               if sequence in _seqPoolDict.keys():
                    _seqPoolDict[sequence].append(key)
               else:
                    _seqPoolDict[sequence]=[key]
          for negChem in negativeChemicals:
               # increment the negative entries
               addSeqEntry(sequence, negChem)

# eliminate confounding chemicals and remove them from dictionary
for sequence in _chemHash.keys():
     for chemical in _chemList:
          if not _chemHash[sequence][chemical]-_chemHash[sequence]["not_"+chemical] >= _chemHash[sequence][chemical]-confoundingThreshold:
               del _chemHash[sequence][chemical]
               del _chemHash[sequence]["not_"+chemical]
          else:
               del _chemHash[sequence]["not_"+chemical] #not_chemical is not confounding or enriching, so delete it
# remove empty entries
toRemove = []
for sequence in _chemHash.keys():
     if not len(_chemHash[sequence]) > 0:
          toRemove.append(sequence)
for sequence in toRemove:
     del _chemHash[sequence]
# remove chemicals that didn't enrich in at least 2 pools
del toRemove[0:len(toRemove)]
for sequence in _chemHash.keys():
     pools = set(_seqPoolDict[sequence])
     if len(pools) < poolThreshold:
          toRemove.append(sequence)
for sequence in _chemHash.keys():
     pools = set(_seqPoolDict[sequence])
     chemicals = set(_chemDict[list(_chemDict.keys())[0]])
     for pool in pools:
          chemicals = chemicals & set(_chemDict[pool])
     if len(chemicals) < poolThreshold:
         toRemove.append(sequence)
     else:
          _chemHash[sequence] = chemicals
for sequence in toRemove:
     if sequence in _chemHash.keys():
          del _chemHash[sequence]
# print out results
for sequence in _chemHash.keys():
     string = sequence      
     for chemical in _chemHash[sequence]:
          string = string+","+ chemical
     for pool in set(_seqPoolDict[sequence]):
          string = string + "," + pool+"="+str(_ratioHash[pool][sequence])
     print(string)  
print ("number of results: "+str(len(_chemHash.keys())))






