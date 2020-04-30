#! /usr/bin/python

# Usage: getSeqsAboveLength.py seqs.fasta/q length

import sys
from Bio import SeqIO
 
length = sys.argv[2]
fileType = sys.argv[1].split(".")[-1]

input_seq_iterator = SeqIO.parse(sys.argv, fileType)
hit_seq_iterator = (record for record in input_seq_iterator if len(record.seq) >= length)
SeqIO.write(hit_seq_iterator, outFile, fileType)

