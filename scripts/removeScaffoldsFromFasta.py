#! /usr/bin/python
'''Command Line:removeScaffoldsFromFasta.py scaffoldFile.fasta ScaffoldID(.txt)'''
import sys

scaffoldFile = sys.argv[1]
scaffoldID = sys.argv[2]

openInFile = open(scaffoldFile, "r")
openOutFile = open("%s_%s-removed.fasta" %(scaffoldFile, scaffoldID), "w")

if ".txt" in scaffoldID:
	infile = open(scaffoldID, "r")
	scaffoldList = []
	writeOut=1
	for scaffoldName in infile:
		scaffoldList.append(scaffoldName.strip(">\n").split(" ")[0])
	for line in openInFile:
		if ">" in line and line.strip(">\n").split(" ")[0] in scaffoldList:
			writeOut = 0
		if ">" in line and line.strip(">\n").split(" ")[0] not in scaffoldList:
			writeOut = 1
			openOutFile.write(line)
		if ">" not in line:
			if writeOut == 1:
				openOutFile.write(line)
else:
	writeOut = 1
	for line in openInFile:
		if ">" in line and line.strip(">\n").split(" ")[0] == scaffoldID:
			writeOut = 0
		if ">" in line and line.strip(">\n").split(" ")[0] != scaffoldID:
			writeOut = 1
			openOutFile.write(line)
		if ">" not in line:
			if writeOut == 1:
				openOutFile.write(line)
openInFile.close()
openOutFile.close()
infile.close()
