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
for line in configFile:
	if line.startswith("#"):
		if len(currentPool)>0:
			_chemDict[currentPool] = currentChemicals;
		currentChemicals = list();
		currentPool = line.split(" ")[1].strip("\n")
		#print "current pool is: " + currentPool
	elif not line.startswith("#") and len(line) > 1:
		#print "adding chemical: "+line.rstrip("\n");		
		currentChemicals.append(line.rstrip("\n"))

##############################################
for key in _chemDict.keys():
	print "current pool: "+key
	for chem in _chemDict[key]:
		print chem


################################################

#initialize chemHash
_chemHash = dict(); 
#this dictionary will contain all of they sequences and the chemicals that enrich them
#key will be the sequence as a string, value will be a list of lists with each of the innerlists containing two elements, the name of the chemical and the number of times that chemical has enriched the given sequence

# this function will initialize the entry for a sequence in 
def addSeqEntry(seq,chem):
	if seq in _chemHash:
		# if key already exists, update it
		print "key already exists: "+seq
	else:
		print("initializing entry: "+seq)
		_chemHash[seq] = seq;


addSeqEntry("moo","")
addSeqEntry("cow","")
addSeqEntry("moo","")























	



files = [];
for fname in dirList:
	files.append(fname)
#for fileName in files:
	#print(fileName)