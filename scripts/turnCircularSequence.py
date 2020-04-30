#! /usr/bin/python

import sys

infile = sys.argv[1]
outname = infile.split(".fasta")[0]
openfile = open(infile, "r")

seq = []
for line in openfile:
	if ">" in line:
		name = line.strip("\n")
	else:
		stripline = line.strip("\n")
		for base in stripline:
			seq.append(base)

openfile.close()
halfway = int(len(seq)/2)
outfile = open("%s_turned.fasta" %(outname), "w")
outfile.write("%s_turned\n" %(name))
outfile.write("".join(seq[halfway:]))
outfile.write("".join(seq[0:halfway]))
outfile.write("\n")
outfile.close()



