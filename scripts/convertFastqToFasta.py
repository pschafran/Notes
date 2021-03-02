#! /usr/bin/python
#
# 2020 January 15. Peter Schafran ps997@cornell.edu

# Usage: convertFastqToFasta.py File.fastq

from Bio import SeqIO
import sys

filename = sys.argv[1]
try:
	if (filename.split(".")[-1]) != "fastq":
		print("ERROR: NOT A FASTQ FILE")
		exit(1)
	else:
		fileprefix = ".".join(filename.split(".")[0:-1])
		commandline = "%s,fastq, %s.fasta, fasta" %(filename, fileprefix)
		SeqIO.convert(filename,"fastq","%s.fasta" %(fileprefix),"fasta")
except:
	print("ERROR: Something went wrong...dumping variables")
	print("filename = %s") % filename
	print("filename-split = %s") % filename.split(".")
	print("fileprefix = %s") % ".".join(filename.split(".")[0:-1])
	print("SeqIO called with: SeqIO.convert(%s)") % commandline 
	exit(1)

