#! /usr/bin/python

import sys


scaffoldFile = sys.argv[1]
outputFile = "%s_unique.fasta" %(scaffoldFile)
infileScaffolds = open(scaffoldFile, "r")
outfile = open(outputFile, "w")
fastaDict = {}

key = 0
fastaDict[key] = []
for line in infileScaffolds:
	if ">" in line:
		joinLine = "".join(fastaDict[key])
		fastaDict[key] = joinLine
		key += 1
		fastaDict[key] = []
	if ">" not in line:
			stripLine = line.strip("\n")
			fastaDict[key].append(stripLine)
joinLine = "".join(fastaDict[key])
fastaDict[key] = joinLine
key = 0
for item in sorted(set(fastaDict.values())):
	outfile.write(">Iengl_Schafran43_scaffold%d\n" %(key))
	outfile.write("%s" %(item))
	outfile.write("\n")
	key += 1

print "%d unique scaffolds out of %d total" %(len(set(fastaDict.values())), len(fastaDict))

outfile.close()
infileScaffolds.close()
