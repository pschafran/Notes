#! /usr/bin/python

import sys

file = sys.argv[1]
fileName = ".".join(file.split(".")[:-1])

infile = open(file, "r")
outfile = open("%s_noGaps.fasta" %(fileName), "w")
for line in infile:
	if ">" in line:
		outfile.write( "\n%s" %(line)) #This makes the first line of the file empty and needs to be manually removed
	else:
		stripline = line.strip("\n")
		for item in stripline:
			if item != "-":
				outfile.write(item)
infile.close()
outfile.close()