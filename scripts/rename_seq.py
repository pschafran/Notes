#!/usr/bin/env python
# Modified by PWS 12 May 2020 to accept infile on commandline

from Bio import SeqIO
import re
import sys

inputfile = sys.argv[1]

def parse_fasta(infile):
	AllSeq = SeqIO.parse(infile, 'fasta')
	return [i for i in AllSeq]

input = parse_fasta('%s' % inputfile)
output = open('%s.renamed.fa' % inputfile, 'w')

counter = 1
for seq in input:
	output.write('>' + str(counter) + '_' + str(seq.id).replace('|','_') + '\n' + str(seq.seq) + '\n')
	counter = counter + 1
