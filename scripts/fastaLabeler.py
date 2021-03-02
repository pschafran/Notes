#! /usr/bin/python
# (C) Peter Schafran 24 March 2017

#Command Line: fastalabeler.py LabelToAppend FastaFile

import sys


append = sys.argv[1]
file = sys.argv[2]

infile = open(file, "r")
outfile = open("%s_labeled.fasta" %(file), "w")

for line in infile:
	line = line.strip('\n')
	if ">" in line:
		outfile.write("%s_%s\n" %(line,append))
	else:
		outfile.write("%s\n" %(line))
		
infile.close()
outfile.close()