#!/usr/bin/python
# 2020 May 13 (C) Peter W Schafran ps997@cornell.edu
#
# Extracts a range of bases out of a FASTA file. 
# Requires a tab-delimited file with columns: contig name, start position, end position
# Usage: trimFastaBases.py sequence.fasta targetFile.txt

import sys

openSeqFile = open(sys.argv[1], "r")
openOutFile = open("%s_extractedRegions.fasta" % sys.argv[1], "w")
openTargetFile = open(sys.argv[2], "r")

# check targetFile
for line in openTargetFile:
	if len(line.strip("\n").split("\t")) < 3:
		print("ERROR: missing elements in line'n")
		print(line)
		exit(1)
	else:
		if int(line.strip("\n").split("\t")[1]) > int(line.strip("\n").split("\t")[2]):
			print("ERROR: Start position comes after end position\n")
			print(line)
			exit(1)
		elif int(line.strip("\n").split("\t")[1]) <= 0 or int(line.strip("\n").split("\t")[2]) <= 0:
			print("ERROR: Start and end position can't be negative or zero\n")
			print(line)
			exit(1)

seqDict = {}
seqList = []
for seqLine in openSeqFile:
	if seqLine.startswith(">"):
		try: 
			if len(seqList) > 1:
				seq = "".join(seqList)
			else:
				seq = seqList[0]
			seqDict[seqName] = seq
		except:
			seqList = []
		seqName = seqLine.strip(">\n").split()[0]
	else:
		seqList.append(seqLine.strip("\n"))
if len(seqList) > 1:
	seq = "".join(seqList)
else:
	seq = seqList[0]
seqDict[seqName] = seq

openTargetFile.seek(0)
for targetLine in openTargetFile:
	targetContig = targetLine.strip(">\n").split("\t")[0]
	startPos = int(targetLine.strip(">\n").split("\t")[1])
	endPos = int(targetLine.strip(">\n").split("\t")[2])
	
	# check for bad positions
	if endPos > len(seqDict[targetContig]):
		print("ERROR: End position is beyond end of sequence!")
		print(targetLine)
		print(len(seqDict[targetContig]))
		exit(1)
	if startPos == 1 and endPos == len(seqDict[targetContig]):
		print("ERROR: Whole sequence selected to be trimmed!")
		print(targetLine)
		exit(1)
	
	openOutFile.write(">%s_pos%s-%s\n" %(targetContig, startPos, endPos))
	baseCount = 0
	writeCount = 0
	for base in seqDict[targetContig]:
		baseCount += 1
		if baseCount >= startPos and baseCount <= endPos:
			openOutFile.write(base)
			writeCount += 1
			if writeCount % 80 == 0:
				openOutFile.write("\n") 
	openOutFile.write("\n")

openSeqFile.close()
openOutFile.close()
openTargetFile.close()