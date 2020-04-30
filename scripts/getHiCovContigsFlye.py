#! /usr/bin/python
#
# Script gets contigs from Flye assembly that have coverage greater than 3 median absolute deviations over the median -- intended to retrieve likely organellar contigs
#
# Usage: getHiCovContigsFlye.py assembly_info.txt assembly.fasta OutputfilePrefix


import sys
import numpy


assembly_info = open(sys.argv[1], "r")
assembly = open(sys.argv[2], "r")
hiCovContigs = open("%s_hiCovContigs.fasta" %(sys.argv[3]), "w")
hiMultContigs = open("%s_hiMultiplicityContigs.fasta" %(sys.argv[3]), "w")
reallyHiCovContigs = open("%s_reallyHiCovContigs.fasta" %(sys.argv[3]), "w")
reallyHiMultContigs = open("%s_reallyHiMultiplicityContigs.fasta" %(sys.argv[3]), "w")
normalContigs = open("%s_normalContigs.fasta" %(sys.argv[3]), "w")

contigDict = {}
covList = []

linecount = 0
for line in assembly_info:
	if linecount > 0:
		splitline = line.strip("\n").split("\t")
		seqName = splitline[0]
		length = int(splitline[1])
		cov = float(splitline[2])
		circular = splitline[3]
		repeat = splitline[4]
		multiplicity = float(splitline[5])
		graphPath = splitline[6]
		
		contigDict[seqName] = {"length":length}
		contigDict[seqName].update({"cov":cov})
		contigDict[seqName].update({"circular":circular})
		contigDict[seqName].update({"repeat":repeat})
		contigDict[seqName].update({"multiplicity":multiplicity})
		contigDict[seqName].update({"graphPath":graphPath})
		
		covList.append(cov)
		
		linecount += 1
	else:
		linecount += 1

covMedian = numpy.median(covList)
covMAD = numpy.median(abs(covList - covMedian))

hiCovContigList = []
hiMultiplicityContigList = []
reallyHiCovContigList = []
reallyHiMultContigList = []

for key in contigDict.keys():
	if contigDict[key]["cov"] > covMedian+(5*covMAD):
		hiCovContigList.append(key)
	if contigDict[key]["cov"] > covMedian+(10*covMAD):
		reallyHiCovContigList.append(key)
	if contigDict[key]["multiplicity"] > 2:
		hiMultiplicityContigList.append(key)
	if contigDict[key]["multiplicity"] > 20:
		reallyHiMultContigList.append(key)
print("High coverage contigs (>5 * median absolute deviation)")
print(hiCovContigList)
print("High multiplicity contigs (>2)")
print(hiMultiplicityContigList)
print("Really high coverage contigs (>10 * absolute median deviation)")
print(reallyHiCovContigList)
print("Really high multiplicity contigs (>20)")
print(reallyHiMultContigList)


writeOut = 0
for line in assembly:
	if ">" in line:
		if line.strip(">\n") in hiCovContigList:
			lineKey = line.strip(">\n")
			hiCovContigs.write(">%s_cov-%s_length-%s_multiplicity-%s\n" %(lineKey, contigDict[lineKey]["cov"], contigDict[lineKey]["length"], contigDict[lineKey]["multiplicity"]))
			writeOut = 1
		else:
			writeOut = 0
	if ">" not in line:
		if writeOut == 1:
			hiCovContigs.write(line)
		elif writeOut == 0:
			pass
assembly.seek(0)
writeOut = 0
for line in assembly:
		if ">" in line:
				if line.strip(">\n") in hiMultiplicityContigList:
						lineKey = line.strip(">\n")
						hiMultContigs.write(">%s_cov-%s_length-%s_multiplicity-%s\n" %(lineKey, contigDict[lineKey]["cov"], contigDict[lineKey]["length"], contigDict[lineKey]["multiplicity"]))
						writeOut = 1
				else:
						writeOut = 0
		if ">" not in line:
				if writeOut == 1:
						hiMultContigs.write(line)
				elif writeOut == 0:
						pass
assembly.seek(0)
writeOut = 0
for line in assembly:
		if ">" in line:
				if line.strip(">\n") in reallyHiCovContigList:
						lineKey = line.strip(">\n")
						reallyHiCovContigs.write(">%s_cov-%s_length-%s_multiplicity-%s\n" %(lineKey, contigDict[lineKey]["cov"], contigDict[lineKey]["length"], contigDict[lineKey]["multiplicity"]))
						writeOut = 1
				else:
						writeOut = 0
		if ">" not in line:
				if writeOut == 1:
						reallyHiCovContigs.write(line)
				elif writeOut == 0:
						pass
assembly.seek(0)
writeOut = 0
for line in assembly:
		if ">" in line:
				if line.strip(">\n") in reallyHiMultContigList:
						lineKey = line.strip(">\n")
						reallyHiMultContigs.write(">%s_cov-%s_length-%s_multiplicity-%s\n" %(lineKey, contigDict[lineKey]["cov"], contigDict[lineKey]["length"], contigDict[lineKey]["multiplicity"]))
						writeOut = 1
				else:
						writeOut = 0
		if ">" not in line:
				if writeOut == 1:
						reallyHiMultContigs.write(line)
				elif writeOut == 0:
						pass
assembly.seek(0)
writeOut = 0
for line in assembly:
		if ">" in line:
				if line.strip(">\n") not in hiMultiplicityContigList and line.strip(">\n") not in hiCovContigList:
						lineKey = line.strip(">\n")
						normalContigs.write(">%s_cov-%s_length-%s_multiplicity-%s\n" %(lineKey, contigDict[lineKey]["cov"], contigDict[lineKey]["length"], contigDict[lineKey]["multiplicity"]))
						writeOut = 1
				else:
						writeOut = 0
		if ">" not in line:
				if writeOut == 1:
						normalContigs.write(line)
				elif writeOut == 0:
						pass


assembly_info.close()
assembly.close()
hiCovContigs.close()
hiMultContigs.close()
reallyHiCovContigs.close()
reallyHiMultContigs.close()
normalContigs.close()
