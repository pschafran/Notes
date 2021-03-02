#! /usr/bin/python

import sys
import gzip
from datetime import datetime
from Bio import SeqIO
from collections import OrderedDict


# This script parses a BLAST output file with TaxIDs, then retrieves and creates new files with the reads/contigs corresponding to each taxonomic group at a specified rank


if "-h" in sys.argv or len(sys.argv) == 1:
	print('''

Help
----
This script parses a BLAST or Kraken2 output file with TaxIDs, then retrieves and creates new files with the reads/contigs corresponding to each taxonomic group at a specified rank

Required
-hits	A BLAST tabular output file (outfmt 6) with query seq id (qseqid) in column 1,  the target seq id (staxids) in column 2, and bitscore in column 3. (Standard Kraken2 output w/ -kraken flag or Centrifuge results w/ -centrifuge flag)
-i	File with reads/contigs used to create BLAST file (FASTA/Q format)
-o	Output file prefix 
-rank	Taxonomic rank to group sequences (must match one of those in nodes DB column 2)

Optional:
-nodesDB	The file containing the taxonomic database. Consists of 4 tab-separated columns: TaxID of node; Taxonomic rank of node; Name of node; TaxID of parent node (default location: /home/ps997/bin/blobtools/nodesDB.txt)
-kraken	Kraken2 formatted hits file
-centrifuge	Centrifuge formatted results file

Dependencies:
Biopython (https://biopython.org)

2020 February 5. Peter Schafran ps997@cornell.edu

''') 
	exit(0)

if "-hits" not in sys.argv or "-i" not in sys.argv or "-o" not in sys.argv or "-rank" not in sys.argv:
	print('''ERROR: One or more parameters missing. Recheck command line syntax. All parameters are required. Run script with -h flag for add'l info.
Usage: taxonomicParser.py -i Reads.fasta/q -rank taxonomic-rank -hits BLAST-hits -o OutputFilePrefix''')
	exit(1)
if "-nodesDB" not in sys.argv:
	nodesDB = "/home/ps997/bin/blobtools/data/nodesDB.txt"
krakenMode = 0
centrifugeMode = 0
print("-"*100)
print("START: %s" %(datetime.now()))
for item in sys.argv:
	if "-i" == item:
		readFile = sys.argv[sys.argv.index(item)+1]
		print("Input file: %s" %(readFile))
		try:
			if readFile.split(".")[-1] == "gz":
				compressed = 1
				fileExt = readFile.split(".")[-2]
				if fileExt == "fastq" or fileExt == "fq":
					fileType = "fastq"
				elif fileExt == "fasta" or fileExt == "fa" or fileExt == "fna" or fileExt == "faa":
					fileType = "fasta"
			else:
				compressed = 0
				fileExt = readFile.split(".")[-1]
				if fileExt == "fastq" or fileExt == "fq":
					fileType = "fastq"
				elif fileExt == "fasta" or fileExt == "fa" or fileExt == "fna" or fileExt == "faa":
					fileType = "fasta"
		except:
			print("ERROR: Input file missing or name formatted improperly!")
			exit(1)
	if "-o" == item:
		outputFilePrefix = sys.argv[sys.argv.index(item)+1]
		print("Output File Prefix: %s" %(outputFilePrefix))
	if "-nodesDB" == item:
		nodesDB = sys.argv[sys.argv.index(item)+1]
	if "-hits" == item:
		blastHits = sys.argv[sys.argv.index(item)+1]
		print("Hits: %s" %(blastHits))
	if "-rank" == item:
		taxonGrouping = sys.argv[sys.argv.index(item)+1]
		print("Taxonomic Grouping: %s" %(taxonGrouping))
	if "-kraken" == item:
		print("Hits file in Kraken2 format")
		krakenMode = 1
	if "-centrifuge" == item:
		print("Hits file in Centrifuge results format")
		centrifugeMode = 1
if krakenMode == 1 and centrifugeMode == 1:
	print("ERROR: Can't select Kraken and Centrifuge formats together!")
	exit(1)
print("-"*100)  

# Parsing nodes DB file
print("Parsing nodes DB file...")
lineCount = 0
nodesDict = {}
openNodesDB = open(nodesDB, "r")
for line in openNodesDB:
	lineCount += 1
	splitline = line.strip("\n").split("\t")
	if lineCount > 1:
		try:
			nodesDict[splitline[0]] = [splitline[1],splitline[2],splitline[3]]
		except:
			print("ERROR on line %s of %s: Incorrectly formatted,ust have 4 tab-separated columns" %(lineCount, nodesDB))
			exit(1)
openNodesDB.close()


# Parsing hits file
blastDict = {} # blastDict[queryID] = {taxID : score}
openBlastHits = open(blastHits, "r")
if krakenMode == 1:
	print("Parsing Kraken hits file...")
	for line in openBlastHits:
		splitline = line.strip("\n").split("\t")
		blastDict[splitline[1]] = {splitline[2] : 100}
elif centrifugeMode == 1:
	centrifugeLineCount = 0
	print("Parsing Centrifuge results file...")
	for line in openBlastHits:
		if centrifugeLineCount > 0:
			splitline = line.strip("\n").split("\t")
			if splitline[0] in blastDict:
				if splitline[2] in blastDict[splitline[0]].keys():
					if blastDict[splitline[0]][splitline[2]] <= float(splitline[3]):
						blastDict[splitline[0]][splitline[2]] = float(splitline[3])
					else:
						pass
				else:
					blastDict[splitline[0]].update({splitline[2] : float(splitline[3])})
			else:
				blastDict[splitline[0]] = {splitline[2] : float(splitline[3])}
			centrifugeLineCount += 1
		else:
			centrifugeLineCount += 1
else:
	print("Parsing BLAST hits file...")
	for line in openBlastHits:
		splitline = line.strip("\n").split("\t")
		if splitline[0] in blastDict:
			if splitline[1] in blastDict[splitline[0]].keys():
				if blastDict[splitline[0]][splitline[1]] <= float(splitline[2]):
					blastDict[splitline[0]][splitline[1]] = float(splitline[2])
				else:
					pass
			else:
				blastDict[splitline[0]].update({splitline[1] : float(splitline[2])})
		else:
			blastDict[splitline[0]] = {splitline[1] : float(splitline[2])}

openBlastHits.close()

# Parsing BLAST dict for greatest bitscore for each hit and enter into taxonomyDict (taxID is key to list of contigs)
print("Parsing results for best hit...")
taxonomyDict = {}
taxIDtotal = 0
taxIDmissing = 0
blastDictKeys = blastDict.keys()
for queryKey in blastDictKeys:
	taxIDtotal += 1
	if krakenMode == 1:
		bestStaxHit = blastDict[queryKey].keys()[0]
		if bestStaxHit == "0":
			bestStaxHit = "1"
	else:	
		sortedStaxDict = OrderedDict(sorted(blastDict[queryKey].items(), key = lambda item: item[1], reverse=True))
		bestStaxHit = sortedStaxDict.keys()[0]
	try:
		currentRank = nodesDict[bestStaxHit][0]
		currentNode = bestStaxHit
		parentNode = nodesDict[bestStaxHit][2]
	except KeyError:
		print("%s taxID not found in NCBI taxonomy" %(queryKey))
		taxIDmissing += 1
	while currentRank != taxonGrouping:
		try:
			currentNode = parentNode
			currentRank = nodesDict[parentNode][0]
			parentNode = nodesDict[parentNode][2]
			if currentNode == 1:
				#print("Reached root without finding corrent taxonomic rank" %(currentNode, parentNode))
				currentRank = taxonGrouping
			elif currentNode > 1 and currentNode == parentNode:
				#print("Error: current node and parent node are the same: %s\t%s" %(currentNode, parentNode))
				currentRank = taxonGrouping
		except KeyError:
			print("Desired taxonomic level not listed for %s" %(queryKey))
	
	if currentNode in taxonomyDict:
		taxonomyDict[currentNode].append(queryKey)
	else:
		taxonomyDict[currentNode] = [queryKey]


# Read in input file
print("Reading input sequences...")
if compressed == 1:
	with gzip.open(readFile, "rt") as readFileGZ:
		inputSeqDict = SeqIO.to_dict(SeqIO.parse(readFileGZ, fileType))
elif compressed == 0:
	inputSeqDict = SeqIO.index(readFile, fileType)


print("Writing output...")
taxonomyDictKeys = taxonomyDict.keys()
for key in taxonomyDictKeys:
	print("%s (TaxID %s):\t%s seqs" %(nodesDict[key][1], key, len(taxonomyDict[key])))
	outfile = open("%s_%s_DEV.%s" %(outputFilePrefix, nodesDict[key][1], fileType), "w")
	for targetSeq in taxonomyDict[key]:
		try:
			record = inputSeqDict[targetSeq]
			SeqIO.write(record, outfile, fileType)
		except KeyError:
			print("%s not found in input seqeunce file" % targetSeq)
	outfile.close()
percTaxIDmissing = (float(taxIDmissing)/float(taxIDtotal))*100
print("taxID not found in NCBI taxonomy for %d out of %d total sequences (%.1f\%)" %(taxIDmissing, taxIDtotal, percTaxIDmissing))




































