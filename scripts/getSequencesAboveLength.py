#! /usr/bin/python

# Usage: getSeqsAboveLength.py seqs.fasta/q length

import sys
from Bio import SeqIO
 
length = int(sys.argv[2])
filetype = sys.argv[1].split(".")[-1]
filename = ".".join(sys.argv[1].split(".")[0:-1])
outfile = open("%s_%s.%s" %(filename, length, filetype), "w")

#print(filetype)
#print(filename)
#print(outfile)

input_seq_dict = SeqIO.index(sys.argv[1], filetype)
for key in input_seq_dict.keys():
	#print(key)
	#print(input_seq_dict[key].seq)
	seqLen = (len(input_seq_dict[key].seq))
	if seqLen >= length:
		SeqIO.write(input_seq_dict[key], outfile, filetype)
outfile.close()
