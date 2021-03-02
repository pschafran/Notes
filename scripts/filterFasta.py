#! /usr/bin/python
'''Command Line:getScaffoldsFromFasta.py SeqFile.fasta FileWithSeqNamesToKeep.txt'''
import sys

seqFile = sys.argv[1]
filterFile = sys.argv[2]


openSeqFile = open(seqFile, "r")
openFilterFile = open(filterFile, "r")
openOutFile = open("%s_%s_filtered.fasta" %(filterFile.split(".")[:-1][0], seqFile.split(".")[:-1][0]), "w")

filterList = []
for line in openFilterFile:
	filterList.append(line.strip("\n"))
openFilterFile.close()
writeOut = 0
for line in openSeqFile:
	splitline = line.strip("\n").split(">")
	try:
		if splitline[1] in filterList:
			writeOut = 1
			openOutFile.write(line)
		if ">" in line and splitline[1] not in filterList:
			writeOut = 0
	except:
		if ">" not in line:
			if writeOut == 1:
				openOutFile.write(line)
		
openSeqFile.close()
openOutFile.close()