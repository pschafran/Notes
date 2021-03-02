#! /opt/local/bin/python2.7
# On my MBP, I had to point to another version of python because the default at /usr/bin/python would not allow updating scipy and all files are locked.

## 2019 February 26
## Last updated 2019 September 23
##Peter W. Schafran
## This script takes a matrix of distances from a tree (only tested on Geneious exported files) and an input list of diploids, 
## then annotates all samples in the matrix with the nearest diploid.

## Input: Requires minimum of two files: 
##		Diploid_List.txt  = A list of diploids in the matrix, and an abbreviation for the diploid that will be appended to sample names. One species per line.
##			E.g. Isoetes_engelmannii, engl
##				 Isoetes_echinospora, echi
##		Distance_Matrix.csv = A matrix of patristic distances from a phylogeny. Only tested on matrices exported from Geneious.
##		Seqs_Or_Alignment_For_Matrix.fasta = A fasta file containing the sequences/alignment used to construct the corresponding phylogeny for the distance matrix. (Optional)

## Output: Generates a list of all samples in the matrix with the nearest diploid species abbreviation appended.
##		E.g. Isoetes_butleri_CiafreSN1_Cluster0;size=29; ==> Isoetes_butleri_CiafreSN1_Cluster0_butl
##		If given a fasta file, the script will also output the fasta with sequence names changed as above.
##		Outputs "Summary" files: one for each unique taxon or taxa combination

## Usage: PURC_assignClustersToDiploids.py Diploid_List.txt Distance_Matrix.csv Seqs_or_Alignment_From_Matrix.fasta (optional)

import sys, time, numpy, scipy
from scipy.optimize import linear_sum_assignment


usage='''
This script requires at least 2 files in this order:
	1. Diploid_List.txt  = A list of diploids in the matrix, and an abbreviation for the diploid that will be appended to sample names. One species per line.
		E.g. 
		Isoetes_engelmannii, engl
		Isoetes_echinospora, echi
	2. Distance_Matrix.csv = A matrix of patristic distances from a phylogeny in comma-separated value format. Only accepts full matrix matrix (header and 1st column with labels, both sides of diagonal with values). 
	3. Seqs_Or_Alignment_For_Matrix.fasta = A fasta file containing the sequences/alignment used to construct the distance matrix (Optional)

Usage: PURC_assignClustersToDiploids.py Diploid_List.txt Distance_Matrix.csv Seqs_or_Alignment_From_Matrix.fasta (optional)
'''



filelist = sys.argv[1:]
try:
	diploidFile = filelist[0]
except:
	print "!" * 10
	print "ERROR: No files provided"
	print "!" * 10
	sys.exit(usage)
matrixFiles = []
fastaFiles = []



#Check the input files for type, correct format
def checkFiles():
	for file in filelist:
		if "csv" in file.split(".")[-1]:
			matrixFiles.append(file)
		elif "fasta" in file.split(".")[-1] or "fa" in file.split(".")[-1]:
			fastaFiles.append(file)
	if len(matrixFiles) == 0:
		print "!" * 10
		print "ERROR: No matrix files provided, or not in CSV format"
		print "!" * 10
		sys.exit(usage)
	if len(fastaFiles) == 0:
			print "-" * 10
			print "WARNING: No fasta files included to be renamed"
			print "-" * 10
			time.sleep(5)
	openfile = open(diploidFile, "r")
	diploidLineCount = 0
	fileCount = 0
	for line in openfile:
		if ">" in line or "\t" in line or "," not in line:
			print "ERROR: Diploid file improperly formatted! Needs to be CSV format"
			print "Error line: %s" %(line)
			time.sleep(5)
		diploidLineCount += 1
	openfile.close()
	print "%s diploid species to be assigned" %(diploidLineCount)
	time.sleep(1)
	lineCount = 0
	for file in matrixFiles:
		fileCount += 1
		openfile = open(file, "r")
		for line in openfile:
			lineCount += 1
		openfile.close()
		lineCount = lineCount - 1
	print "%s matrix file(s) with approximately %s clusters to assign" %(fileCount, lineCount)
	time.sleep(1)

checkFiles()

def assignClusters():
	for file in matrixFiles:
		print "Working on %s..." %(file)
		time.sleep(1)
		matrixOpenFile = open(file, "r")
		diploidOpenFile = open(diploidFile)
		diploidList = []
		diploidAbbrevDict = {}
		indexDict = {}
		distanceDict = {}
		
		outputFileName = file.split(".")[0]
		outfile = open("%s_renamed.txt" %(outputFileName), "w")
		
		for line in diploidOpenFile:
			splitline = line.strip("\n").split(",")
			diploidList.append(splitline[0])
			diploidAbbrevDict[splitline[0]] = splitline[1]

		lineCount = 0
		for line in matrixOpenFile:
			if lineCount == 0:
				indexSplitLine = line.strip("\n").split(",")
				indexCount = 0
				for item in indexSplitLine:
					if item in diploidList:
						indexDict[indexCount] = item
					indexCount += 1
			else:
				distanceList = []
				matrixSplitLine = line.strip("\n").split(",")
				for key in indexDict.keys():
					if matrixSplitLine[0] == indexDict[key]:
						pass
					else:
						if matrixSplitLine[0] in distanceDict.keys():
							distanceDict[matrixSplitLine[0]].update({matrixSplitLine[key] : indexDict[key]})
						else:
							distanceDict[matrixSplitLine[0]] = {matrixSplitLine[key] : indexDict[key]}
			lineCount += 1

		newNameDict = {}
		allCombosDict = {}
		polyploidCombosDict = {}
		diploidCombosDict = {}

		for key in distanceDict.keys():
			if key not in diploidList:
				minDistance = sorted(distanceDict[key].keys())[0]
				closestDiploid = distanceDict[key][minDistance]
				splitKey = key.split(";")[0]
				newName = "%s_%s\n" %(splitKey, diploidAbbrevDict[closestDiploid])
				outfile.write(newName)
				newNameDict[key] = newName
				sampleName = key.split("_Cluster")[0]
				try:
					allCombosDict[sampleName].append(diploidAbbrevDict[closestDiploid])
				except:
					allCombosDict[sampleName] = [diploidAbbrevDict[closestDiploid]]
# Separate allCombos dict into dict with just diploids (1 diploid seq) and polyploids (>1 polyploid seq)
		for comboKey in allCombosDict.keys():
			if len(allCombosDict[comboKey]) == 1:
				diploidCombosDict[comboKey] = "_".join(sorted(allCombosDict[comboKey]))
			elif len(allCombosDict[comboKey]) > 1:
				polyploidCombosDict[comboKey] = "_".join(sorted(allCombosDict[comboKey]))
			allCombosDict[comboKey] = "_".join(sorted(allCombosDict[comboKey]))

# Create "Summary" output files and write headers
		summaryOutfile = open("Summary_%s" %(file), "w")
		summaryOutfile.write("Diploid OTUs,Samples\n")
		summaryPolyploidFile = open("Summary_Polyploids_%s" %(file), "w")
		summaryPolyploidFile.write("Diploid OTUs, Samples\n")
		summaryDiploidFile = open("Summary_Diploids_%s" %(file), "w")
		summaryDiploidFile.write("Diploid OTU, Samples\n")

# Loop through dict of all taxa in dict allCombosDict[sample-name] = [diploid1_diploid2...diploid_n]
		for taxon in sorted(set(allCombosDict.values())):
# Write out txt file with all samples for each unique diploid combination
			taxonOutfile = open("%s.txt" %(taxon), "w")
			taxonOutfile.write("%s\n" %(taxon))
			summaryOutfile.write("%s" %(taxon))
			for comboKey in sorted(allCombosDict.keys(), key = str):
				if allCombosDict[comboKey] == taxon:
					taxonOutfile.write("%s\n" %(comboKey))
					summaryOutfile.write(",%s" %(comboKey))
			taxonOutfile.close()
			summaryOutfile.write("\n")
		summaryOutfile.close()

# Write out txt file for samples with each unique polyploid combo
		for polyploidTaxon in sorted(set(polyploidCombosDict.values())):
			summaryPolyploidFile.write(polyploidTaxon)
			for polyploidComboKey in sorted(polyploidCombosDict.keys(), key = str):
				if polyploidCombosDict[polyploidComboKey] == polyploidTaxon:
					summaryPolyploidFile.write(",%s" %(polyploidComboKey))
			summaryPolyploidFile.write("\n")
		summaryPolyploidFile.close()

# Write out txt file for samples with each unique diploid
		for diploidTaxon in sorted(set(diploidCombosDict.values())):
			summaryDiploidFile.write(diploidTaxon)
			for diploidComboKey in sorted(diploidCombosDict.keys(), key = str):
				if diploidCombosDict[diploidComboKey] == diploidTaxon:
					summaryDiploidFile.write(",%s" %(diploidComboKey))
			summaryDiploidFile.write("\n")
		summaryDiploidFile.close()

# Write out csv file for UpSet plot (https://cran.r-project.org/web/packages/UpSetR/vignettes/basic.usage.html)
# File is a diploid X sample matrix with presence/absence values for each diploid seq from a given sample
		upsetOut = open("UpSet_Polyploids%s" %(file), "w")
		upsetOut.write("Sample,")
		for diploid in sorted(set(diploidCombosDict.values())):
			upsetOut.write("%s," %(diploid))
		upsetOut.write("\n")
		for sample in sorted(polyploidCombosDict.keys()):
			upsetOut.write("%s," %(sample))
			for diploid in sorted(set(diploidCombosDict.values())):
				if diploid in polyploidCombosDict[sample]:
					upsetOut.write("1,")
				elif diploid not in polyploidCombosDict[sample]:
					upsetOut.write("0,")
				else:
					print "ERROR around line 195: This shouldn't happen. Something was formatted wrong in input files"
					exit(1)
			upsetOut.write("\n")
		upsetOut.close()

# Create dictionary of dictionaries to hold each sample with distances to each other sample
		matrixDistanceDict = {} # matrixDistanceDict[sampleID] = {sampleID1: dist, sampleID2: dist, ...}
		matrixOpenFile.seek(0)
		matrixList = [] # distance matrix stored as a list of lists with structure matrixList[row1][column1]
		matrixLine = 0
		for line in matrixOpenFile:
			if matrixLine == 0:
				headerLine = line.strip("\n").split(",")
				matrixList.append(headerLine)
				matrixLine += 1
			else:
				sampleLine = line.strip("\n").split(",")
				matrixList.append(sampleLine)
				matrixDistanceDict[sampleLine[0]] = {}
				count = 0
				for item in sampleLine:
					matrixDistanceDict[sampleLine[0]].update({headerLine[count]:sampleLine[count]})
					count += 1
		
# Loop through distance dictionary to extract branch lengths for sets of sequences that belong to a pair of samples. Write minimum OTU distance to all-by-all sample matrix

		polyploidDistOutfile = open("Polyploid_Sample_Distance_%s" %(file), "w")
		polyploidDistOutfile.write(" ,%s\n" %(sorted(set(polyploidCombosDict.keys()))))
		matrixPositionDict = {}

## Loop through list of samples and create new dict to hold sample-by-sample pairwise OTU distances
		for sample1 in sorted(set(polyploidCombosDict.keys())):
			matrixPositionDict[sample1] = {}
			otuMatrixColumnPositions = []

### Loop through OTUs belonging to sample
			for otu in headerLine:
				if sample1 in otu: # Match determined by sample name occurring in OTU name so no partial matches can be present in sample names
					otuMatrixColumnPositions.append(headerLine.index(otu))
			matrixPositionDict[sample1] = {"ColumnPositions" : otuMatrixColumnPositions}

### Loop though`list of samples again to do all-by-all comparison
			for sample2 in sorted(set(polyploidCombosDict.keys())):
				otuMatrixRowPositions = []
				if sample2 != sample1:

#### Loop through OTUs belonging to sample 2
					for otu2 in headerLine:
						if sample2 in otu2: # Match determined by sample name occurring in OTU name so no partial matches can be present in sample names
							otuMatrixRowPositions.append(headerLine.index(otu2))
				matrixPositionDict[sample1][sample2] = otuMatrixRowPositions

## Loop through sample lists again to get pairs of samples, then extract matrix positions from matrixPositionDict
		for sample1 in sorted(set(polyploidCombosDict.keys())):
			polyploidDistOutfile.write("%s," %(sample1))
			for sample2 in sorted(set(polyploidCombosDict.keys())):
				if sample2 in sample1 or sample1 in sample2:
					polyploidDistOutfile.write(",")
				else:
					tempColumns = matrixPositionDict[sample1]["ColumnPositions"]
					tempRows = matrixPositionDict[sample1][sample2]
					tempMatrix = []
### Matrix positions used to extract values from matrixList
					for x in tempColumns:
						tempRowDist = []
						for y in tempRows:
							tempRowDist.append(matrixList[x][y])
						tempMatrix.append(tempRowDist)
### List of lists with matrix values converted to array
					arrayName = numpy.array(tempMatrix, dtype = float)
### Hungarian Method applied to matrix to get shortest possible branch length sum and results written to file
					row_ind, col_ind = linear_sum_assignment(arrayName)
					if arrayName[row_ind, col_ind].sum() < 0:
						print "WARNING: Negative distance between %s and %s" %(sample1, sample2)
					polyploidDistOutfile.write("%s," %(arrayName[row_ind, col_ind].sum()))
			polyploidDistOutfile.write("\n")
		polyploidDistOutfile.close()
		
# Write out fasta with diploid abbreviations appended to sequence names
		if len(fastaFiles) > 0:
			try:
				inFasta = open("%s.fasta" %(outputFileName), "r")
				outFasta = open("%s_renamed.fasta" %(outputFileName), "w")
			except:
				inFasta = open("%s.fa" %(outputFileName), "r")
				outFasta = open("%s_renamed.fa" %(outputFileName), "w")
			writeOut = 0
			for line in inFasta:
				if ">" in line:
					writeOut = 0
					stripline = line.strip("\n").split(">")[1]
					if stripline in newNameDict.keys():
						writeOut = 1
						outFasta.write(">%s" %(newNameDict[stripline]))
					else:
						print "%s missing from newNameDict" %(stripline)
				else:
					if writeOut == 1:
						outFasta.write(line)
			inFasta.close()
		try:
			outFasta.close()
		except:
			pass
		outfile.close()

assignClusters()

