#! /usr/bin/python

import sys


fasta_input = sys.argv[1]
inScaffolds = open(fasta_input, "r")
openOutput = open("MissingRemoved_%s" %(fasta_input), "w")
inScaffoldDict = {}
def readFasta():
	for line in inScaffolds:
		if ">" in line:
			try:
				joinLine = "".join(inScaffoldDict[scaffoldName]) # This must be before next line to merge previously parsed seq
				inScaffoldDict[scaffoldName] = joinLine
			except:
				pass
			scaffoldName = line
			inScaffoldDict[scaffoldName] = []
		if ">" not in line:
			stripLine = line.strip("\n")
			inScaffoldDict[scaffoldName].append(stripLine)
	joinLine = "".join(inScaffoldDict[scaffoldName]) # This line must be after for loop to join final sequence
	inScaffoldDict[scaffoldName] = joinLine
	inScaffolds.close()
	
def removeMissing():
	uniqueScaffoldDict = {}
	for key in inScaffoldDict.keys():
		uniqueScaffoldDict[key] = []
		for site in inScaffoldDict[key]:
			uniqueScaffoldDict[key].append(site)
		if len(set(uniqueScaffoldDict[key])) > 1:
			openOutput.write(key)
			openOutput.write("%s\n" %(inScaffoldDict[key]))
	openOutput.close()
readFasta()
removeMissing()
