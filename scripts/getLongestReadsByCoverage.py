#! /usr/bin/python
#

from Bio import SeqIO
from collections import OrderedDict
from datetime import datetime
import sys

if "-h" in sys.argv or len(sys.argv) == 1:
	print('''

Help
----
getLongestReadsByCoverage.py retrieves the longest reads (from a FASTQ file) that will give a certain amount of coverage for the genome, and puts those reads in a new FASTQ file.

Required Parameters:
-c	The desired coverage for the genome
-g	The genome size for the organism. Can use shorthand (e.g. 150k, 200m, 3.2g)
-i	Input file (fastq format)
-o	Output file (fastq format, must have different name than input)

Dependencies:
Biopython (https://biopython.org)

2020 January 17. Peter Schafran ps997@cornell.edu

''') 
	exit(0)

if "-i" not in sys.argv or "-g" not in sys.argv or "-c" not in sys.argv or "-o" not in sys.argv:
	print('''ERROR: One or more parameters missing. Recheck command line syntax. All parameters are required.
Usage: getLongestReadsByCoverage.py -i Reads.fastq -g GenomeSize -c DesiredCoverage -o OutputFilename.fastq''')
	exit(1)

print("-"*100)
print("START: %s" %(datetime.now()))
for item in sys.argv:
	if "-i" == item:
		readFile = sys.argv[sys.argv.index(item)+1]
		try:
			testOpen = open(readFile, "r")
			testOpen.close()
			fileExt = readFile.split(".")[-1]
			if fileExt != "fastq" and fileExt != "fq":
				print("ERROR: Input file doesn't appear to be a FASTQ file!")
				print("File extension seen: %s " % fileExt)
				exit(1)
			else:
				print("Input file: %s" %(readFile))
		except:
			print("ERROR: Input file missing or name formatted improperly!")
			exit(1)
	if "-c" == item:
		coverage = float(sys.argv[sys.argv.index(item)+1])
		print("Coverage: %5.1f" %(coverage))
	if "-g" == item:
		genomeSize = sys.argv[sys.argv.index(item)+1]
		print("Genome Size: %s" %(genomeSize))
		if "k" in genomeSize:
			genomeSize = 1000*float(genomeSize.split("k")[0])
		elif "K" in genomeSize:
			genomeSize = 1000*float(genomeSize.split("K")[0])
		elif "m" in genomeSize:
			genomeSize = 1000000*float(genomeSize.split("m")[0])
		elif "M" in genomeSize:
			genomeSize = 1000000*float(genomeSize.split("M")[0])
		elif "g" in genomeSize:
			genomeSize = 1000000000*float(genomeSize.split("g")[0])
		elif "G" in genomeSize:
			genomeSize = 1000000000*float(genomeSize.split("G")[0])
		else:
			genomeSize = float(genomeSize)
		
	if "-o" == item:
		outputFile = sys.argv[sys.argv.index(item)+1]
		try:
			fileExt = readFile.split(".")[-1]
			if fileExt != "fastq" and fileExt != "fq":
				print("ERROR: Output file doesn't appear to be a FASTQ file!")
				exit(1)
			else:
				print("Output file: %s" %(outputFile))
		except:
			print("ERROR: Output file missing or name formatted improperly!")
			exit(1)
print("-"*100)

if readFile == outputFile:
	print("ERROR: Read file and output file have the same name!")
	exit(1)

print("Parsing input file...this may take a long time!")

seqLengthDict = {}
recordCounter = 0

inputFastqDict = SeqIO.index(readFile, "fastq")
inputFastqDictKeys = inputFastqDict.keys()
inputFastqDictLen = float(len(inputFastqDict))
for key in inputFastqDictKeys:
	seqLengthDict[key] = int(len(inputFastqDict[key].seq))
	# Progress Bar
	recordCounter += 1
	completionPerc = int(float(100*recordCounter)/inputFastqDictLen)
	sys.stdout.write('\r')
	sys.stdout.write("[%-100s] %d%%" % ('='*completionPerc, completionPerc))
	sys.stdout.flush()
print("\n")




sortedSeqDict = OrderedDict(sorted(seqLengthDict.items(), key = lambda item: item[1], reverse=True))
print("Selecting longest reads to get target coverage...")
bpSum = 0
targetSeqList = []
sortedSeqDictKeys = sortedSeqDict.keys()
for seq in sortedSeqDictKeys:
	currentCoverage = float(bpSum/genomeSize)
	if currentCoverage < coverage:
		targetSeqList.append(seq)
		bpSum = bpSum + sortedSeqDict[seq]
		currentCoverage = float(bpSum/genomeSize)
		#print("\tCurrent coverage: %5.1f" %(currentCoverage))
	else:
		pass
print("Found %d longest seqs provide %5.1f X coverage" %(len(targetSeqList),currentCoverage))

print("Writing sequences to output...this may take a long time!")

baseCounter = 0
recordCounter = 0
completionPerc = 0

outfile = open(outputFile, "w")
targetSeqListLen = float(len(targetSeqList))
for targetSeq in targetSeqList:
	record = inputFastqDict[targetSeq]
	SeqIO.write(record, outfile, "fastq")
	
	#Progress Bar
	recordCounter += 1
	completionPerc = int(float(100*recordCounter)/targetSeqListLen)
	sys.stdout.write('\r')
	sys.stdout.write("[%-100s] %d%%" % ('='*completionPerc, completionPerc))
	sys.stdout.flush()
print("\n")
print("END: %s" %(datetime.now()))
outfile.close()

