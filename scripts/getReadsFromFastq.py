#! /usr/bin/python
#
# 2020 January 13. Peter Schafran ps997@cornell.edu
#
# Usage: getReadsFromFastq.py AllReads.fastq ReadIDs.txt OutputReads.fastq
# Files must be called in order!

import sys
from Bio import SeqIO

if len(sys.argv) != 4:
	print('''ERROR: Incorrect number of files specified! Need three files in this order.
Usage: getReadsFromFastq.py AllReads.fastq ReadIDs.txt OutputReads.fastq
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

	input_seq_iterator = SeqIO.parse(readFile, "fastq")
	hit_seq_iterator = (record for record in input_seq_iterator \
					   if record.name in hitList)
	SeqIO.write(hit_seq_iterator, outFile, "fastq")