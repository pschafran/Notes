#! /home/ps997/miniconda3/bin/python

# Usage: centrifugeExtractor.py centrifuge.results taxIDsToReject.txt reads.fastq

import sys
from Bio import SeqIO

centResultFile = sys.argv[1]
taxIDsFile = sys.argv[2]
readsFile = sys.argv[3]

taxIDsList = []
keepReadsList = []
rejectReadsList = []

openTaxIDsFile = open(taxIDsFile, "r")
for line in openTaxIDsFile:
	taxIDsList.append(line.strip("\n"))
openTaxIDsFile.close()


openCentResultFile = open(centResultFile , "r")
for line in openCentResultFile:
	splitline = line.strip("\n").split("\t")
	if splitline[2] in taxIDsList:
		rejectReadsList.append(splitline[0])
	elif splitline[2] not in taxIDsList:
		keepReadsList.append(splitline[0])
openCentResultFile.close()

outputKeepFile = open("%s.keep.fastq" % readsFile, "w")
outputRejectFile = open("%s.reject.fastq" % readsFile, "w")


input_seq_dict = SeqIO.index(readsFile, "fastq")
print(input_seq_dict)
#for record in input_seq_dict:
#	 if record in keepReadsList:
#		SeqIO.write(record, outputKeepFile, "fastq")
#	else:
#		SeqIO.write(record, outputRejectFile, "fastq")