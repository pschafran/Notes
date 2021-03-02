#! /usr/bin/python

import sys

#buscoFile = "full_table_BUSCO.tsv"
#scaffoldFile = "scaffolds.fasta"

#singleCopyBUSCOs = "Itaiw_BUSCOs.fna"

#infileBUSCO = open(buscoFile, "r")
#infileScaffolds = open(scaffoldFile, "r")
#infileSCBUSCOs = open(singleCopyBUSCOs, "r")


commandLine = sys.argv[1:]
argList = ["-h", "--help", "-a", "--align", "-i", "--inter"]
filelist = []
for item in commandLine:
	if item not in argList:
		filelist.append(item)

def help():
	if "-h" in commandLine or "--help" in commandLine:
		print '''
Purpose:
This script combines the same benchmarking universal single-copy orthologs (BUSCOs) from the single_copy_busco_sequence directories output by BUSCO (busco.ezlab.org).
Requires that FASTA names in BUSCO output files have sample identifier in sample name (e.g. >E0902352SLX:Sample1:scaffold46145). Easiest to do by naming contig/scaffold file before running BUSCO.

Usage:
getCompleteBUSCOs.py [flags] ./path/to/files.fna

Dependencies:
MAFFT (https://mafft.cbrc.jp/) installed and in PATH. 
	On ODU HPC, enable with:
	enable_lmod
	module load mafft/7

Flags:
-h, --help	Display this help menu
-a, --align	Align combined BUSCOs
-i, --inter	Save unaligned sequence files. Assumes -a
-s, --stats	Generate missing data statistics. Assumes each subfolder represents a different sample.


Written by Peter Schafran	peterwschafran.com	November 2018
'''
		sys.exit()

def dictCreate():
	counter = 0
	for line in infileBUSCO:
		if "Complete" in line:
			splitline = line.strip("\n").split("\t")
			buscoList.append(splitline[0])
			contigList.append(splitline[2])
			buscoDict[splitline[0]] = splitline[2]
			counter += 1
	print "Found %d complete BUSCOs" %(counter)

def scaffoldSearch():
	writeOut = 0
	for line in infileScaffolds:
		scaffoldID = line.strip("\n").split(">")
		if ">" in line and scaffoldID[1] in buscoDict.values():
			#print "Found scaffold %s" %(scaffoldID[1])
			writeOut = 1
			for key in buscoDict.keys():
				if scaffoldID[1] in buscoDict[key]:
					outfile.write(">%s=%s\n" %(key,scaffoldID[1]))
		if ">" in line and scaffoldID[1] not in buscoDict.values():
			writeOut = 0
		if ">" not in line:
			if writeOut == 1:
				outfile.write(line)

def checkFile():
	counter = 0
	for key in sorted(buscoDict.keys()):
		print key
		for line in infileScaffolds:
			if key in line:
				counter += 1
				print "Found %d of %d matches" %(counter, len(buscoDict.keys()))
		infileScaffolds.seek(0)

def checkFileFormats():
	filenameList = []
	for filename in filelist:
		splitFilename = filename.split(".")
		filenameList.append(splitFilename[-1])
	if len(set(filenameList)) > 1:
		print'''
WARNING: Multiple file extensions in file list!
'''
	if len(set(filenameList)) == 0:
		print'''
ERROR: No input files specified!
'''

def combineBUSCOs():
	for file in filelist:
		openFile = open(file, "r")
		for line in openFile:
			if ">" in line:
				splitline = line.strip("\n").split(":")
				if splitline[0] not in scbuscoDict.keys():
					scbuscoDict[splitline[0]] = {}
		openFile.close()

def BUSCOscaffoldSearch():
	writeOut = 0
	for file in filelist:
		splitFile = file.split(".")
		fileExt = splitFile[-1]
		openFile = open(file, "r")
		for line in openFile:
			scaffoldID = line.strip("\n").split(":")
			if ">" in line and scaffoldID[0] in scbuscoDict.keys():
				#print "Found scaffold %s" %(scaffoldID[0])
				buscoID = scaffoldID[0]
				fileID = scaffoldID[1] + "/" + scaffoldID[2]
				writeOut = 1
				scbuscoDict[buscoID][fileID] = {}
				statsDict[sampleList].append(scaffoldID[1])
			if ">" in line and scaffoldID[0] not in scbuscoDict.keys():
				writeOut = 0
			if ">" not in line:
				if writeOut == 1:
					scbuscoDict[buscoID][fileID] = line
		openFile.close()
	statsDict[numSamples] = len(set(statsDict[sampleList]))
	for key in scbuscoDict.keys():
		numSeqs = len(scbuscoDict[key])
		if numSeqs > 0:
			splitKey = key.split(">")
			statsDict[missingDataList].append(100-((float(len(scbuscoDict[key]))/float(statsDict[numSamples]))*100))
			outFileName = "%s_combined_%sseqs.%s" %(splitKey[1],numSeqs,fileExt)
			outFile = open(outFileName, "w")
			for sample in scbuscoDict[key]:
				sampleName = sample.split("/")
				outFile.write(">%s\n" %(sampleName[1]))
				outFile.write(scbuscoDict[key][sample])
			outFile.close()
		
def calcStats():
	statsDict[minMissingData] = min(statsDict[missingDataList])
	statsDict[meanMissingData] = sum(statsDict[missingDataList])/len(statsDict[missingDataList])
	statsDict[maxMissingData] = max(statsDict[missingDataList])
	if len(statsDict[missingDataList]) % 2 == 0:
		medianIndex = len(statsDict[missingDataList])/2
		sortedMissingData = sorted(statsDict[missingDataList])
		statsDict[medianMissingData] = sortedMissingData[medianIndex]
	elif len(statsDict[missingDataList]) % 2 != 0:
		medianIndex = int(len(statsDict[missingDataList])/2)
		sortedMissingData = sorted(statsDict[missingDataList])
		statsDict[medianMissingData] = sortedMissingData[medianIndex]
	print "Number of BUSCOs:\t%d" %(len(scbuscoDict.keys()))
	print "Minimum Missing Data:\t%f" %(statsDict[minMissingData])
	print "Mean Missing Data:\t%f" %(statsDict[meanMissingData])
	print "Median Missing Data:\t%f" %(statsDict[medianMissingData])
	print "Maximum Missing Data:\t%f" %(statsDict[maxMissingData])
help()
buscoDict={}
scbuscoDict={}
statsDict={}
sampleList = "sampleList"
numSamples = "numSamples"
minMissingData = "minMissingData"
maxMissingData = "maxMissingData"
meanMissingData = "meanMissingData"
medianMissingData = "medianMissingData"
missingDataList = "missingDataList"
statsDict.keys().append(sampleList)
statsDict.keys().append(numSamples)
statsDict.keys().append(minMissingData)
statsDict.keys().append(maxMissingData)
statsDict.keys().append(meanMissingData)
statsDict.keys().append(medianMissingData)
statsDict[sampleList] = []
statsDict[missingDataList] = []
contigList = []
buscoList = []
scbuscoList=[]
#dictCreate()
#scaffoldSearch()
#checkFile()
checkFileFormats()
combineBUSCOs()
BUSCOscaffoldSearch()
calcStats()
