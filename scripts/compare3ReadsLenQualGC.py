#! /home/ps997/miniconda3/bin/python
# Needs to run with python 3 on our server!!!

# 2020 January 18
# Peter W. Schafran ps997@cornell.edu

from Bio import SeqIO
from Bio.SeqUtils import GC
import pkg_resources
pkg_resources.require("numpy==1.17.3")  # modified to use specific numpy
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import sys
from matplotlib import colors
from numpy import log

### File parsing ###

fileDict = {}

for i in range(len(sys.argv)):
	if i == 0:
		pass
	elif i > 0:
		filePrefix = sys.argv[i].split(".fastq")[0]
		if "/" in filePrefix:
			filePrefix = filePrefix.split("/")[-1]
		fileDict[sys.argv[i]] = {"filePrefix":filePrefix}
		fileDict[sys.argv[i]].update({"lengthList":[]})
		fileDict[sys.argv[i]].update({"logLengthList":[]})
		fileDict[sys.argv[i]].update({"qualList":[]})
		fileDict[sys.argv[i]].update({"gcList":[]})

fileCounter = 1
for key in fileDict.keys():
	print("Parsing file %d of %d...%s" %(fileCounter, len(fileDict.keys()), key))
	fileIndex = SeqIO.index(key, "fastq")
	recordCounter = 0
	for read in fileIndex.keys():
		record = fileIndex[read]
		fileDict[key]["lengthList"].append(len(record.seq))
		fileDict[key]["logLengthList"].append(log(len(record.seq)))
		fileDict[key]["qualList"].append(sum(record.letter_annotations["phred_quality"])/len(record.letter_annotations["phred_quality"]))
		fileDict[key]["gcList"].append(GC(record.seq))
	
		#Progress bar
		recordCounter += 1
		completionPerc = int(float(100*recordCounter)/float(len(fileIndex)))
		sys.stdout.write('\r')
		sys.stdout.write("[%-100s] %d%%" % ('='*completionPerc, completionPerc))
		sys.stdout.flush()
	fileCounter += 1
	print("\n")
	
	
#file1 = sys.argv[1]
#file1Prefix = file1.split(".fastq")[0]
#if "/" in file1Prefix:
#	file1Prefix = file1Prefix.split("/")[-1]
#file2 = sys.argv[2]
#file2Prefix = file2.split(".fastq")[0]
#if "/" in file2Prefix:
#	file2Prefix = file2Prefix.split("/")[-1]
#file3 = sys.argv[3]
#file3Prefix = file3.split(".fastq")[0]
#if "/" in file3Prefix:
#	file3Prefix = file3Prefix.split("/")[-1]


#run1LengthList = []
#run2LengthList = []
#run3LengthList = []
#run1QualList = []
#run2QualList = []
#run3QualList = []
#run1GCList = []
#run2GCList = []
#run3GCList = []


	
#print("Parsing file 1...")
#run1Index = SeqIO.index(file1, "fastq")
#recordCounter = 0
#for key in run1Index.keys():
#	record = run1Index[key]
#	run1LengthList.append(len(record.seq))
#	run1QualList.append(sum(record.letter_annotations["phred_quality"])/len(record.letter_annotations["phred_quality"]))
#	run1GCList.append(GC(record.seq))
#	
#	#Progress bar
#	recordCounter += 1
#	completionPerc = int(float(100*recordCounter)/float(len(run1Index)))
#	sys.stdout.write('\r')
#	sys.stdout.write("[%-100s] %d%%" % ('='*completionPerc, completionPerc))
#	sys.stdout.flush()
#print("\n")

#print("Parsing file 2...")
#run2Index = SeqIO.index(file2, "fastq") 
#recordCounter = 0
#for key in run2Index.keys():
#	record = run2Index[key]
#	run2LengthList.append(len(record.seq))
#	run2QualList.append(sum(record.letter_annotations["phred_quality"])/len(record.letter_annotations["phred_quality"]))
#	run2GCList.append(GC(record.seq))
#	
#	#Progress bar
#	recordCounter += 1
#	completionPerc = int(float(100*recordCounter)/float(len(run2Index)))
#	sys.stdout.write('\r')
#	sys.stdout.write("[%-100s] %d%%" % ('='*completionPerc, completionPerc))
#	sys.stdout.flush()
#print("\n")

#print("Parsing file 3...")
#run3Index = SeqIO.index(file3, "fastq") 
#recordCounter = 0
#for key in run3Index.keys():
#	record = run3Index[key]
#	run3LengthList.append(len(record.seq))
#	run3QualList.append(sum(record.letter_annotations["phred_quality"])/len(record.letter_annotations["phred_quality"]))
#	run3GCList.append(GC(record.seq))
#	
#	#Progress bar
#	recordCounter += 1
#	completionPerc = int(float(100*recordCounter)/float(len(run3Index)))
#	sys.stdout.write('\r')
#	sys.stdout.write("[%-100s] %d%%" % ('='*completionPerc, completionPerc))
#	sys.stdout.flush()
#print("\n")

# Log transform lengths
#run1LogLength = [log(record) for record in run1LengthList]
#run2LogLength = [log(record) for record in run2LengthList] 
#run3LogLength = [log(record) for record in run3LengthList] 


### Plotting stuff ###
print("Plotting stuff...")
numFiles = len(fileDict.keys())
fontScale = 1/numFiles
plt.tick_params(axis='x', labelsize=15*fontScale)


# Histograms of lengths for both runs
fig, axs = plt.subplots(1, numFiles, sharex=True, sharey=True, tight_layout=False)
axsIndex = 0
for key in sorted(fileDict.keys()):
	axs[axsIndex].hist(fileDict[key]["lengthList"], bins = 100)
	axsIndex += 1
axsIndex = 0
for key in sorted(fileDict.keys()):
	if axsIndex == 0:
		axs[axsIndex].set(xlabel='Read Length', ylabel='No. of Reads')
		axs[axsIndex].set_title(fileDict[key]["filePrefix"], size=20*fontScale)
	else:
		axs[axsIndex].set(xlabel='Read Length')
		axs[axsIndex].set_title(fileDict[key]["filePrefix"], size=20*fontScale)
	axsIndex += 1
plt.savefig("Sequence_lengths.pdf" , format = "pdf")
plt.close()


# Histograms of log-lengths for both runs
fig, axs = plt.subplots(1, numFiles, sharex=True, sharey=True, tight_layout=False)
axsIndex = 0
for key in sorted(fileDict.keys()):
	axs[axsIndex].hist(fileDict[key]["logLengthList"], bins = 100)
	axsIndex += 1
axsIndex = 0
for key in sorted(fileDict.keys()):
	if axsIndex == 0:
		axs[axsIndex].set(xlabel='Log-transformed Read Length', ylabel='No. of Reads')
		axs[axsIndex].set_title(fileDict[key]["filePrefix"], fontsize = 20*fontScale)
	else:
		axs[axsIndex].set(xlabel='Log-transformed Read Length')
		axs[axsIndex].set_title(fileDict[key]["filePrefix"], fontsize = 20*fontScale)
	axsIndex += 1
plt.savefig("Log_sequence_lengths.pdf" , format = "pdf")
plt.close()
#fig, axs = plt.subplots(1, 3, sharex=True, sharey=True, tight_layout=False)
#axs[0].hist(run1LogLength, bins = 100)
#axs[1].hist(run2LogLength, bins = 100)
#axs[2].hist(run3LogLength, bins = 100)
#axs[0].set(xlabel='Log-transformed Read Length', ylabel='No. of Reads', title=file1)
#axs[1].set(xlabel='Log-transformed Read Length', title=file2)
#axs[2].set(xlabel='Log-transformed Read Length', title=file3)
#plt.savefig("%s_%s_%s_loglength.pdf" %(file1Prefix, file2Prefix, file3Prefix), format = "pdf")
#plt.close()

# Histograms of quality for both runs
fig, axs = plt.subplots(1, numFiles, sharex=True, sharey=True, tight_layout=False)
axsIndex = 0
for key in sorted(fileDict.keys()):
	axs[axsIndex].hist(fileDict[key]["qualList"], bins = 100)
	axsIndex += 1
axsIndex = 0
for key in sorted(fileDict.keys()):
	if axsIndex == 0:
		axs[axsIndex].set(xlabel='Read Avg. Quality Score', ylabel='No. of Reads', title=fileDict[key]["filePrefix"], fontsize = 15*fontScale)
	else:
		axs[axsIndex].set(xlabel='Read Avg. Quality Score', title=fileDict[key]["filePrefix"], fontsize = 15*fontSize)
	axsIndex += 1
plt.savefig("Sequence_qualities.pdf" , format = "pdf")
plt.close()

#fig, axs = plt.subplots(1, 3, sharex=True, sharey=True, tight_layout=False)
#axs[0].hist(run1QualList, bins = 100)
#axs[1].hist(run2QualList, bins = 100)
#axs[2].hist(run3QualList, bins = 100)
#axs[0].set(xlabel='Read Avg. Quality Score', ylabel='No. of Reads', title=file1)
#axs[1].set(xlabel='Read Avg. Quality Score', title=file2)
#axs[2].set(xlabel='Read Avg. Quality Score', title=file3)
#plt.savefig("%s_%s_%s_quality.pdf" %(file1Prefix, file2Prefix, file3Prefix), format = "pdf")
#plt.close()

# Histograms of GC% for both runs
fig, axs = plt.subplots(1, 3, sharex=True, sharey=True, tight_layout=False)
axs[0].hist(run1GCList, bins = 100)
axs[1].hist(run2GCList, bins = 100)
axs[2].hist(run3GCList, bins = 100)
axs[0].set(xlabel='Read GC%', ylabel='No. of Reads', title=file1)
axs[1].set(xlabel='Read GC%', title=file2)
axs[2].set(xlabel='Read GC%', title=file3)
plt.savefig("%s_%s_%s_GCperc.pdf" %(file1Prefix, file2Prefix, file3Prefix), format = "pdf")
plt.close()

# 2D Histograms of length vs. quality for both runs

fig, axs = plt.subplots(1, 3, sharey=True, tight_layout=False)
axs[0].hist2d(run1LengthList, run1QualList, bins=100, norm = colors.LogNorm())
axs[1].hist2d(run2LengthList, run2QualList, bins=100, norm = colors.LogNorm())
h = axs[2].hist2d(run3LengthList, run3QualList, bins=100, norm = colors.LogNorm())
axs[0].set(xlabel='Read Length', ylabel='Read Avg. Quality Score', title=file1)
axs[1].set(xlabel='Read Length', title=file2)
axs[2].set(xlabel='Read Length', title=file3)
cbar = plt.colorbar(h[3], ax=axs[2], pad = 0.01, aspect = 50)
cbar.set_label('# of reads', rotation=90, fontsize = 10)
cbar.ax.tick_params(labelsize=5) 
plt.savefig("%s_%s_%s_length_vs_qual.pdf" %(file1Prefix, file2Prefix, file3Prefix), format = "pdf")
plt.close()

# 2D Histogram of log-length vs. quality for both runs
fig, axs = plt.subplots(1, 3, sharex=True, sharey=True, tight_layout=False)
axs[0].hist2d(run1LogLength, run1QualList, bins=100)
axs[1].hist2d(run2LogLength, run2QualList, bins=100)
h = axs[2].hist2d(run3LogLength, run3QualList, bins=100)
axs[0].set(xlabel='Log-transformed Read Length', ylabel='Read Avg. Quality Score', title=file1)
axs[1].set(xlabel='Log-transformed Read Length', title=file2)
axs[2].set(xlabel='Log-transformed Read Length', title=file3)
cbar = plt.colorbar(h[3], ax=axs[2], pad = 0.01, aspect = 50)
cbar.set_label('# of reads', rotation=90, fontsize = 10)
cbar.ax.tick_params(labelsize=5) 
plt.savefig("%s_%s_%s_loglength_vs_qual.pdf" %(file1Prefix, file2Prefix, file3Prefix), format = "pdf")
plt.close()

# 2D Histograms of length vs. GC% for both runs
fig, axs = plt.subplots(1, 3, sharey=True, tight_layout=False)
axs[0].hist2d(run1LengthList, run1GCList, bins=100, norm = colors.LogNorm())
axs[1].hist2d(run2LengthList, run2GCList, bins=100, norm = colors.LogNorm())
h = axs[2].hist2d(run3LengthList, run3GCList, bins=100, norm = colors.LogNorm())
axs[0].set(xlabel='Read Length', ylabel='Read GC%', title=file1)
axs[1].set(xlabel='Read Length', title=file2)
axs[2].set(xlabel='Read Length', title=file3)
cbar = plt.colorbar(h[3], ax=axs[2], pad = 0.01, aspect = 50)
cbar.set_label('# of reads', rotation=90, fontsize = 10)
cbar.ax.tick_params(labelsize=5) 
plt.savefig("%s_%s_%s_length_vs_GCperc.pdf" %(file1Prefix, file2Prefix, file3Prefix), format = "pdf")
plt.close()

# 2D Histogram of log-length vs. GC% for both runs
fig, axs = plt.subplots(1, 3, sharex=True, sharey=True, tight_layout=False)
axs[0].hist2d(run1LogLength, run1GCList, bins=100, norm = colors.LogNorm())
axs[1].hist2d(run2LogLength, run2GCList, bins=100, norm = colors.LogNorm())
h = axs[2].hist2d(run3LogLength, run3GCList, bins=100, norm = colors.LogNorm())
axs[0].set(xlabel='Log-transformed Read Length', ylabel='Read GC%', title=file1)
axs[1].set(xlabel='Log-transformed Read Length', title=file2)
axs[2].set(xlabel='Log-transformed Read Length', title=file3)
cbar = plt.colorbar(h[3], ax=axs[2], pad = 0.01, aspect = 50)
cbar.set_label('# of reads', rotation=90, fontsize = 10)
cbar.ax.tick_params(labelsize=5) 
plt.savefig("%s_%s_%s_loglength_vs_GCperc.pdf" %(file1Prefix, file2Prefix, file3Prefix), format = "pdf")
plt.close()

# 2D Histograms of quality vs. GC% for both runs
fig, axs = plt.subplots(1, 3, sharey=True, tight_layout=False)
axs[0].hist2d(run1QualList, run1GCList, bins=100, norm = colors.LogNorm())
axs[1].hist2d(run2QualList, run2GCList, bins=100, norm = colors.LogNorm())
h = axs[2].hist2d(run3QualList, run3GCList, bins=100, norm = colors.LogNorm())
axs[0].set(xlabel='Read Avg. Quality Score', ylabel='Read GC%', title=file1)
axs[1].set(xlabel='Read Avg. Quality Score', title=file2)
axs[2].set(xlabel='Read Avg. Quality Score', title=file3)
cbar = plt.colorbar(h[3], ax=axs[2], pad = 0.01, aspect = 50)
cbar.set_label('# of reads', rotation=90, fontsize = 10)
cbar.ax.tick_params(labelsize=5) 
plt.savefig("%s_%s_%s_quality_vs_GCperc.pdf" %(file1Prefix, file2Prefix, file3Prefix), format = "pdf")
plt.close()

