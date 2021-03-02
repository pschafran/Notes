#! /usr/bin/python
##2019 February 26
##Peter W. Schafran

### From PURC output files, this removes clusters that represent <10% of each sample's total reads (hardcoded around line 68), and creates CSVs with percentage of each cluster size.
### From root PURC directory, to generate from all different analyses run: PURC_cleanup.py ./purc_out_*/purc_run*4*

import sys, itertools, os


filelist = sys.argv[1:]
clusterFileList = []
for file in filelist:
	openfile = open(file, 'r')
	filename = file.split("/")
	fileout = "ClusterStats_%s.csv" %(filename[-1])
	clusterFileList.append(fileout)
	outfile = open(fileout, 'w')

	sampleDict = {}
	clusterDict = {}
	percentDict = {}

	##Loop through fasta cluster file ("purc_run_...4_clustered_reconsensus.fa") to extract cluster sizes from cluster names. 
	for line in openfile:
		if ">" in line:
			clusterDict = {}
			splitline = line.strip("\n").split(";")
			splitline2 = splitline[0].split("_")
			sampleName = "_".join(splitline2[0:-1])
			clusterName = splitline2[-1:][0]
			clusterSize = splitline[1].split("=")[1]
			clusterDict = {clusterName : clusterSize}
			if sampleName in sampleDict.keys():
				sampleDict[sampleName][clusterName] = clusterSize
			
			else:
				sampleDict[sampleName] = clusterDict

	##Calculate the max number of clusters for the dataset
	maxClusters = 1
	for sampleName in sampleDict.keys():
		if len((sampleDict[sampleName])) > maxClusters:
			maxClusters = len(sampleDict[sampleName].keys())

	##Calculate percentages for each cluster
	for sampleName in sampleDict.keys():
		readTotal = 0 
		numClusters = 0
		for cluster in sampleDict[sampleName]:
			readTotal += int(sampleDict[sampleName][cluster])
			numClusters += 1
		percentDict[sampleName] = {"NumClusters" : numClusters}
		percentDict[sampleName].update({"ReadTotal" : readTotal})
		for cluster in sampleDict[sampleName]:
			clusterPercent = round((float(sampleDict[sampleName][cluster])/float(readTotal))*100, 2)
			percentDict[sampleName].update({cluster : clusterPercent})
	
	##Write new .fasta with just clusters with >10% of reads (non-contaminants)
	fastaFileout = "CleanClusters_%s" %(filename[-1])
	cleanFasta = open(fastaFileout , "w")
	openfile.close()
	openfile = open(file, "r")
	for line in openfile:
		if ">" in line:
			writeFasta = 0
			splitline = line.strip("\n").split(";")
			splitline2 = splitline[0].split("_")
			sampleName = "_".join(splitline2[0:-1])
			clusterName = splitline2[-1:][0]
#### Hardcoded cutoff of 10% for threshold to save sequences
			if percentDict[sampleName][clusterName] > 10:
				writeFasta = 1
				cleanFasta.write(line)
		elif ">" not in line and writeFasta == 1:
			cleanFasta.write(line)
	cleanFasta.close()

	##Write header line to output file
	header = []
	header.append("SampleName")
	header.append("NumClusters")
	for x in range(maxClusters):
		header.append("Cluster%s" %x)

	header.append("TotalReads")
	for x in range(maxClusters):
		header.append("Cluster%s_Percent" %x)
	header.append("\n")
	outfile.write(",".join(header))

	##Write sample lines to output file
	for sampleName in sampleDict.keys():
		count = 1
		writeline = []
		writeline.append(sampleName)
		writeline.append(str(percentDict[sampleName]["NumClusters"]))
		for cluster in sampleDict[sampleName].keys():
			writeline.append(sampleDict[sampleName][cluster])
			count += 1
	##While loop needed to add extra space to make TotalReads and Percent columns to line up
		while count <= maxClusters:
			writeline.append(" ")
			count += 1
		writeline.append(str(percentDict[sampleName]["ReadTotal"]))
		for cluster in sampleDict[sampleName].keys():
			writeline.append(str(percentDict[sampleName][cluster]))
		writeline.append("\n")
		outfile.write(",".join(writeline))
	openfile.close()
	outfile.close()
	
	
##Start combining individual ClusterStats files in current directory
masterSampleDict = {}
uniqueSampleList = []
outfile = open("ClusterStats_allTotalClusters.csv", "w")
for file in clusterFileList:
	openfile = open(file, 'r')
	for line in openfile:
		if line.split(",")[0] not in uniqueSampleList and line.split(",")[0] != "SampleName":
			uniqueSampleList.append(line.split(",")[0])
for sample in uniqueSampleList:
	masterSampleDict[sample] = {}
	for file in clusterFileList:
		counter = 0
		openfile = open(file, "r")
		for line in openfile:
			if sample in line:
				line = line.split(",")
				masterSampleDict[sample].update({file: line[1]})
				counter += 1
		if counter == 0:
			masterSampleDict[sample].update({file: "0"})

header = []
header.append("SampleName")
for file in sorted(clusterFileList):
	header.append(file)
header.append("\n")
outfile.write(",".join(header))
for sample in sorted(masterSampleDict.keys()):
	writeline = []
	writeline.append(sample)
	for item in sorted(masterSampleDict[sample]):
		writeline.append(str(masterSampleDict[sample][item]))
	writeline.append("\n")
	outfile.write(",".join(writeline))

openfile.close()
outfile.close()






