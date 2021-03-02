#! /usr/bin/python

# Usage: renameFasta.py augustus.hits.gtf contig_conversion_table.tsv

# Name table is 2-column TSV with old name in col 1 and new name in col 2

import sys


fasta = sys.argv[1]
names = sys.argv[2]

filename = ".".join(fasta.split(".")[:-1])
fileext = fasta.split(".")[-1]

seqNameDict = {}

with open(names, "r") as inNames:
	for line in inNames:
		oldName = line.strip("\n").split("\t")[0]
		newName = line.strip("\n").split("\t")[1]
		seqNameDict[oldName] = newName
outfile = open("%s_renamed_contigs.%s" % (filename, fileext), "w")
with open(fasta, "r") as inSeqs:
	for line in inSeqs:
		if line.startswith("#"):
			outfile.write(line)
		else:	
			oldName = line.strip("\n").split("\t")[0]
			restOfLine = "\t".join(line.strip("\n").split("\t")[1:])
			try:
				outfile.write("%s\t%s\n" % (seqNameDict[oldName], restOfLine) )
			except KeyError:
				print("Name missing from conversion table: %s" % oldName)
outfile.close()
