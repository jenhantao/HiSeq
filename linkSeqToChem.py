import os
path="./data"  # insert the path to the directory of interest here
dirList=os.listdir(path)
for fname in dirList:
    print fname
