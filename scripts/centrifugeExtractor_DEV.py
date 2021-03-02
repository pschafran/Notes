#! /home/ps997/miniconda3/bin/python

# Usage: centrifugeExtractor.py centrifuge.results taxIDsToReject.txt readsfile.fq

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

keepReadsList = set(keepReadsList)
rejectReadsList = set(rejectReadsList)

outputKeepReadsList = open("%s_keepReadsList.txt" % readsFile, "w")
for readID in keepReadsList:
	outputKeepReadsList.write("%s\n" % readID)
outputKeepReadsList.close()
print("Keep-reads file created")

outputKeepFile = open("%s.keep.DEV.fastq" % readsFile, "w")
outputRejectFile = open("%s.reject.DEV.fastq" % readsFile, "w")

input_seq_dict = SeqIO.to_dict(SeqIO.parse(readsFile, "fastq"))
for record in input_seq_dict:
	if record in keepReadsList:
		SeqIO.write(input_seq_dict[record], outputKeepFile, "fastq")
	elif record in rejectReadsList:
		SeqIO.write(input_seq_dict[record], outputRejectFile, "fastq")
	else:
		print("%s not present in centrifuge results! Did you use the wrong fastq?" % record)
outputKeepFile.close()
outputRejectFile.close()
