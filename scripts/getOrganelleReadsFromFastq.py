#! /usr/bin/python
#
# 2020 January 14. Peter Schafran ps997@cornell.edu
#
# Usage: getOrganelleReadsFromFastq.py AllReads.fastq OrganelleBLASTuniqueReadIDs.txt OutputPrefix
# Files must be called in order!

import sys
from Bio import SeqIO

if len(sys.argv) != 4:
	print('''ERROR: Incorrect number of files specified! Need three files in this order.
Usage: getOrganelleReadsFromFastq.py AllReads.fastq OrganelleBLASTuniqueReadIDs.txt OutputPrefix
Files must be in this order!
''')
	exit(1)
else:
	readFile = sys.argv[1]
	hitFile = sys.argv[2]
	outFile = sys.argv[3]

	hitList = []
	infile = open(hitFile, "r")
	for line in infile:
		hitList.append(line.strip("\n"))
		
	organelleFile = open("%s_organelle.fastq" %(outFile), "w")
	nuclearFile = open("%s_nuclear.fastq" %(outFile), "w")

	for record in SeqIO.parse(readFile, "fastq"):
		if record.name in hitList:
			SeqIO.write(record, organelleFile, "fastq")
		else:
			SeqIO.write(record, nuclearFile, "fastq")

infile.close()
organelleFile.close()
nuclearFile.close()
