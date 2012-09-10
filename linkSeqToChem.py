import os
path="./data"  # insert the path to the directory of interest here
# I want to process the .his values which are just raw counts of every sequence that appears
# I want to process the .val files which are derviced from the .proc files
# formatting for .val files [.proc fileName] [enrichedCount] [originalCount] [enrichmentRatio]
dirList=os.listdir(path)
files = [];
for fname in dirList:
    files.append(fname)
for fileName in files:
    print(fileName)
