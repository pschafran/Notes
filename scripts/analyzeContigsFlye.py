#! /usr/bin/python
#
#
# Usage: analyzeContigsFlye.py assembly_info.txt assembly.fasta output-prefix-from-taxonomicParser outputfilename.csv


import sys
import numpy
import glob
import re
from Bio import SeqIO
from Bio.SeqUtils import GC
from Bio.SeqUtils import lcc


assembly_info = open(sys.argv[1], "r")
assembly = sys.argv[2]
taxonomyFilePrefix = sys.argv[3]
outfile = open(sys.argv[4], "w")

inputSeqDict = SeqIO.index(assembly, "fasta")

contigDict = {}
covList = []
linecount = 0
for line in assembly_info:
	if linecount > 0:
		splitline = line.strip("\n").split("\t")
		seqName = splitline[0]
		length = int(splitline[1])
		cov = float(splitline[2])
		covList.append(cov)
		if splitline[3] == "+":
			circular = 1
		else:
			circular = 0
		if splitline[4] == "+":
			repeat = 1
		else:
			repeat = 0
		multiplicity = float(splitline[5])
		gcContent = GC(inputSeqDict[seqName].seq)
		lccList = lcc.lcc_mult(inputSeqDict[seqName].seq, 1000)
		lccMean = numpy.mean(lccList)
		lccMedian = numpy.median(lccList)
		lccMin = numpy.min(lccList)
		lccMax = numpy.max(lccList)
		lccStdev = numpy.std(lccList)
		
		contigDict[seqName] = {"length":length}
		contigDict[seqName].update({"cov":cov})
		contigDict[seqName].update({"circular":circular})
		contigDict[seqName].update({"repeat":repeat})
		contigDict[seqName].update({"multiplicity":multiplicity})
		contigDict[seqName].update({"gcContent":gcContent})
		contigDict[seqName].update({"lccMin":lccMin})
		contigDict[seqName].update({"lccMean":lccMean})
		contigDict[seqName].update({"lccMedian":lccMedian})
		contigDict[seqName].update({"lccMax":lccMax})
		contigDict[seqName].update({"lccStdev":lccStdev})
		
		linecount += 1
	else:
		linecount += 1

covMedian = numpy.median(covList)
covMAD = numpy.median(abs(covList - covMedian))

for taxonomyFile in glob.glob("%s*.fasta" %(taxonomyFilePrefix)):
	taxon = re.split("_|.fasta", taxonomyFile)[-2]
	openTaxFile = open(taxonomyFile, "r")
	for taxLine in openTaxFile:
		if taxLine.startswith(">"):
			contig = taxLine.strip(">\n")
			contigDict[contig].update({"taxon":taxon})
	openTaxFile.close()

outfile.write("seqID,length,cov,relative_cov,circular,repeat,multiplicity,GC,minLCC,meanLcc,medianLCC,maxLCC,stdevLCC,taxon\n")
for key in contigDict.keys():
	try:
		outfile.write("%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n" %(key, contigDict[key]["length"], contigDict[key]["cov"], float(((contigDict[key]["cov"] - covMAD) / covMAD)*100), contigDict[key]["circular"], contigDict[key]["repeat"], contigDict[key]["multiplicity"], contigDict[key]["gcContent"], contigDict[key]["lccMin"], contigDict[key]["lccMean"], contigDict[key]["lccMedian"], contigDict[key]["lccMax"], contigDict[key]["lccStdev"], contigDict[key]["taxon"]))
	except KeyError:
		outfile.write("%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,unclassified\n" %(key, contigDict[key]["length"], contigDict[key]["cov"], float(((contigDict[key]["cov"] - covMAD) / covMAD)*100), contigDict[key]["circular"], contigDict[key]["repeat"], contigDict[key]["multiplicity"], contigDict[key]["gcContent"], contigDict[key]["lccMin"], contigDict[key]["lccMean"], contigDict[key]["lccMedian"], contigDict[key]["lccMax"], contigDict[key]["lccStdev"]))
assembly_info.close()
outfile.close()
