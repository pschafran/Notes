#! /usr/bin/python

import sys

seqfile = sys.argv[1]
filename = ".".join(seqfile.split(".")[:-1])
fileext = seqfile.split(".")[-1]

openfile = open(seqfile, "r")
seqDict = {}

for line in openfile:
	if line.startswith(">"):
		seqID = line.strip(">\n").split(" ")[0]
		geneID = seqID.split(".")[0]
		transcriptID = seqID.split(".")[1]
		if geneID in seqDict.keys():
			seqDict[geneID].update({transcriptID : []})
		else:
			seqDict[geneID] = {transcriptID : []}
	else:
		seqDict[geneID][transcriptID].append(line.strip("\n"))
	seqDict[geneID][transcriptID] = ["".join(seqDict[geneID][transcriptID])]

openfile.close()
outfile = open("%s_primary_transcripts.%s" %(filename, fileext), "w")
for gene in sorted(seqDict.keys()):
	maxlen = 0
	maxTranscript = "t1"
	for transcript in seqDict[gene].keys():
		if maxlen < len(seqDict[gene][transcript][0]):
			maxlen = len(seqDict[gene][transcript][0])
			maxTranscript = transcript
	outfile.write(">%s.%s\n" %(gene, maxTranscript))
	outfile.write("%s\n" % seqDict[gene][maxTranscript][0])
outfile.close()
