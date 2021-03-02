#! /usr/bin/python
'''Command Line:
getFromFasta.py sourceFile.fasta seqID[.txt]
OR
pipeline | getFromFasta.py sourceFile.fasta - [ > outfile.fasta ]
'''
import sys
import os.path

scaffoldFile = sys.argv[1]
scaffoldID = sys.argv[2]
openInFile = open(scaffoldFile, "r")

sourceDict = {}

for line in openInFile:
	if line.startswith(">"):
		lineName = line.strip(">\n").split(" ")[0]
		sourceDict.update({lineName : []})
		try:
			seq = "".join(sourceDict[previousLineName])
			sourceDict[previousLineName] = seq
		except:
			pass
	else:
		previousLineName = lineName
		sourceDict[lineName].append(line.strip("\n"))
seq = "".join(sourceDict[previousLineName])
sourceDict[previousLineName] = seq

if os.path.isfile(scaffoldID):
	infile = open(scaffoldID, "r")
	for line in infile:
		if line.startswith("#"):
			pass
		else:
			name = line.strip(">\n").split(" ")[0]
			print(">" + name)
			print(sourceDict[name])
elif scaffoldID != "-":
	name = scaffoldID.strip(">\n").split(" ")[0]
	print(">" + name)
	print(sourceDict[name])
elif scaffoldID == "-":
	for line in sys.stdin:
		if line.startswith("#"):
			pass
		else:
			name = line.strip(">\n").split(" ")[0]
			print(">" + name)
			print(sourceDict[name])

openInFile.close()
