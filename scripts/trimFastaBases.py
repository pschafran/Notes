#!/usr/bin/python
# 2020 April 30 (C) Peter W Schafran ps997@cornell.edu
#
# Trims a range of bases out of a single sequence in FASTA format. Base range is inclusive (start and end positions are removed).
# Usage: trimFastaBases.py singleSequence.fasta startBasePosition endBasePosition

import sys

openSeqFile = open(sys.argv[1], "r")
openOutFile = open("%s_trimmed.fasta" % sys.argv[1], "w")
startPos = int(sys.argv[2])
endPos = int(sys.argv[3])

if startPos > endPos:
	print("ERROR: Start position is after end position!")
	exit(1)

seqList = []

for line in openSeqFile:
	if line.startswith(">"):
		seqName = line.strip("\n")
	else:
		seqList.append(line.strip("\n"))
if len(seqList) > 1:
	seq = "".join(seqList)
else:
	seq = seqList[0]

if endPos > len(seq):
	print("ERROR: End position is beyond end of sequence!")
	exit(1)
if startPos == 1 and endPos == len(seq):
	print("ERROR: Whole sequence selected to be trimmed!")
	exit(1)
baseCount = 0
writeCount = 0
openOutFile.write("%s\n" % seqName)
for base in seq:
	baseCount += 1
	if baseCount < startPos or baseCount > endPos:
		openOutFile.write(base)
		writeCount += 1
		if writeCount % 80 == 0:
			openOutFile.write("\n") 
openOutFile.write("\n")
openSeqFile.close()
openOutFile.close()

