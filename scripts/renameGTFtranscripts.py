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
outfile = open("%s_renamed_genes.%s" % (filename, fileext), "w")
with open(fasta, "r") as inSeqs:
	for line in inSeqs:
		if line.startswith("#"):
			outfile.write(line)
		else:	
			commentField = line.strip("\n").split("\t")[8]
			comments = commentField.split(";")
			if len(comments) == 1:
				oldName = comments[0]
				restOfLine = "\t".join(line.strip("\n").split("\t")[0:8])
				try:
                                	outfile.write("%s\t%s\n" % (restOfLine, seqNameDict[oldName]) )
                        	except KeyError:
                                	print("Name missing from conversion table: %s" % oldName)
			else:
				transcriptField = comments[0]
				geneField = comments[1]
				transcriptID = comments[0].split("\"")[1]
                        	transcriptNum = transcriptID.split(".")[1]
                        	geneID = comments[1].split("\"")[1]				
				restOfLine = "\t".join(line.strip("\n").split("\t")[0:8])

				try:
					outfile.write('''%s\ttranscript_id "%s.%s"; gene_id "%s"\n''' % (restOfLine, seqNameDict[oldName], transcriptNum, seqNameDict[oldName]))
				except KeyError:
					print("Name missing from conversion table: %s" % oldName)
outfile.close()
