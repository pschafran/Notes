#! /home/ps997/miniconda3/bin/python
# Needs to run with python 3 on our server!!!

# 2020 January 18
# Peter W. Schafran ps997@cornell.edu

from Bio import SeqIO
from Bio.SeqUtils import GC
import pkg_resources
#pkg_resources.require("numpy==1.17.3")  # modified to use specific numpy
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import sys
from matplotlib import colors
from numpy import log
from datetime import datetime

print(datetime.now())

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
	fileIndexKeys = fileIndex.keys()
	for read in fileIndexKeys:
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
	

### Plotting stuff ###
print("Plotting stuff...")
numFiles = len(fileDict.keys())
fontScale = 1/numFiles

if numFiles > 1:
	# Histograms of lengths for both runs
	fig, axs = plt.subplots(1, numFiles, sharex=True, sharey=True, tight_layout=False)
	axsIndex = 0
	for key in sorted(fileDict.keys()):
		axs[axsIndex].hist(fileDict[key]["lengthList"], bins = 100)
		axsIndex += 1
	axsIndex = 0
	for key in sorted(fileDict.keys()):
		if axsIndex == 0:
			axs[axsIndex].set_ylabel('No. of Reads')
			axs[axsIndex].set_xlabel('Read Length', size=20*fontScale, wrap=True)
			axs[axsIndex].tick_params(axis='x', labelsize=20*fontScale)
			axs[axsIndex].set_title(fileDict[key]["filePrefix"], size=20*fontScale, wrap=True)
		else:
			axs[axsIndex].set_xlabel('Read Length', size=20*fontScale, wrap=True)
			axs[axsIndex].tick_params(axis='x', labelsize=20*fontScale)
			axs[axsIndex].set_title(fileDict[key]["filePrefix"], size=20*fontScale, wrap=True)
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
			axs[axsIndex].set_ylabel('No. of Reads')
			axs[axsIndex].tick_params(axis='x', labelsize=20*fontScale)
			axs[axsIndex].set_xlabel('Log-transformed Read Length', size=20*fontScale, wrap = True)
			axs[axsIndex].set_title(fileDict[key]["filePrefix"], size = 20*fontScale, wrap = True)
		else:
			axs[axsIndex].set_xlabel('Log-transformed Read Length', size=20*fontScale, wrap = True)
			axs[axsIndex].tick_params(axis='x', labelsize=20*fontScale)
			axs[axsIndex].set_title(fileDict[key]["filePrefix"], size = 20*fontScale, wrap = True)
		axsIndex += 1
	plt.savefig("Log_sequence_lengths.pdf" , format = "pdf")
	plt.close()

	# Histograms of quality for both runs
	fig, axs = plt.subplots(1, numFiles, sharex=True, sharey=True, tight_layout=False)
	axsIndex = 0
	for key in sorted(fileDict.keys()):
		axs[axsIndex].hist(fileDict[key]["qualList"], bins = 100)
		axsIndex += 1
	axsIndex = 0
	for key in sorted(fileDict.keys()):
		if axsIndex == 0:
			axs[axsIndex].set_ylabel('No. of Reads')
			axs[axsIndex].tick_params(axis='x', labelsize=20*fontScale)
			axs[axsIndex].set_xlabel('Read Avg. Quality Score', size=20*fontScale, wrap = True)
			axs[axsIndex].set_title(fileDict[key]["filePrefix"], size = 20*fontScale, wrap = True)
		else:
			axs[axsIndex].set_xlabel('Read Avg. Quality Score', size=20*fontScale, wrap = True)
			axs[axsIndex].tick_params(axis='x', labelsize=20*fontScale)
			axs[axsIndex].set_title(fileDict[key]["filePrefix"], fontsize = 20*fontScale, wrap = True)
		axsIndex += 1
	plt.savefig("Sequence_qualities.pdf" , format = "pdf")
	plt.close()


	# Histograms of GC% for both runs
	fig, axs = plt.subplots(1, numFiles, sharex=True, sharey=True, tight_layout=False)

	axsIndex = 0
	for key in sorted(fileDict.keys()):
		axs[axsIndex].hist(fileDict[key]["gcList"], bins = 100)
		axsIndex += 1
	axsIndex = 0
	for key in sorted(fileDict.keys()):
		if axsIndex == 0:
			axs[axsIndex].set_ylabel('No. of Reads')
			axs[axsIndex].tick_params(axis='x', labelsize=20*fontScale)
			axs[axsIndex].set_xlabel('GC%', size=20*fontScale, wrap = True)
			axs[axsIndex].set_title(fileDict[key]["filePrefix"], size = 20*fontScale, wrap = True)
		else:
			axs[axsIndex].set_xlabel('GC%', size=20*fontScale, wrap = True)
			axs[axsIndex].tick_params(axis='x', labelsize=20*fontScale)
			axs[axsIndex].set_title(fileDict[key]["filePrefix"], fontsize = 20*fontScale, wrap = True)
		axsIndex += 1
	plt.savefig("GC_content.pdf" , format = "pdf")
	plt.close()


	# 2D Histograms of length vs. quality for both runs

	fig, axs = plt.subplots(1, numFiles, sharex=True, sharey=True, tight_layout=False)
	axsIndex = 0
	for key in sorted(fileDict.keys()):
		h = axs[axsIndex].hist2d(fileDict[key]["lengthList"], fileDict[key]["qualList"], bins = 100, norm = colors.LogNorm())
		axsIndex += 1
	axsIndex = 0
	for key in sorted(fileDict.keys()):
		if axsIndex == 0:
			axs[axsIndex].set_ylabel('Read Avg. Quality Score')
			axs[axsIndex].tick_params(axis='x', labelsize=20*fontScale)
			axs[axsIndex].set_xlabel('Read Length', size=20*fontScale, wrap = True)
			axs[axsIndex].set_title(fileDict[key]["filePrefix"], size = 20*fontScale, wrap = True)
		else:
			axs[axsIndex].set_xlabel('Read Length', size=20*fontScale, wrap = True)
			axs[axsIndex].tick_params(axis='x', labelsize=20*fontScale)
			axs[axsIndex].set_title(fileDict[key]["filePrefix"], fontsize = 20*fontScale, wrap = True)
		axsIndex += 1
	cbar = plt.colorbar(h[3], ax=axs[axsIndex-1], pad = 0.01, aspect = 50)
	cbar.set_label('# of reads', rotation=90, fontsize = 10)
	cbar.ax.tick_params(labelsize=5) 
	plt.savefig("Length_vs_Quality.pdf" , format = "pdf")
	plt.close()
	
	
	# 2D Histogram of log-length vs. quality for both runs
	fig, axs = plt.subplots(1, numFiles, sharex=True, sharey=True, tight_layout=False)
	axsIndex = 0
	for key in sorted(fileDict.keys()):
		h = axs[axsIndex].hist2d(fileDict[key]["logLengthList"], fileDict[key]["qualList"], bins = 100, norm = colors.LogNorm())
		axsIndex += 1
	axsIndex = 0
	for key in sorted(fileDict.keys()):
		if axsIndex == 0:
			axs[axsIndex].set_ylabel('Read Avg. Quality Score')
			axs[axsIndex].tick_params(axis='x', labelsize=20*fontScale)
			axs[axsIndex].set_xlabel('Log-transformed Read Length', size=20*fontScale, wrap = True)
			axs[axsIndex].set_title(fileDict[key]["filePrefix"], size = 20*fontScale, wrap = True)
		else:
			axs[axsIndex].set_xlabel('Log-transformed Read Length', size=20*fontScale, wrap = True)
			axs[axsIndex].tick_params(axis='x', labelsize=20*fontScale)
			axs[axsIndex].set_title(fileDict[key]["filePrefix"], fontsize = 20*fontScale, wrap = True)
		axsIndex += 1
	cbar = plt.colorbar(h[3], ax=axs[axsIndex-1], pad = 0.01, aspect = 50)
	cbar.set_label('# of reads', rotation=90, fontsize = 10)
	cbar.ax.tick_params(labelsize=5) 
	plt.savefig("LogLength_vs_Quality.pdf" , format = "pdf")
	plt.close()
	
	
	# 2D Histograms of length vs. GC% for both runs
	fig, axs = plt.subplots(1, numFiles, sharex=True, sharey=True, tight_layout=False)
	axsIndex = 0
	for key in sorted(fileDict.keys()):
		h = axs[axsIndex].hist2d(fileDict[key]["lengthList"], fileDict[key]["gcList"], bins = 100, norm = colors.LogNorm())
		axsIndex += 1
	axsIndex = 0
	for key in sorted(fileDict.keys()):
		if axsIndex == 0:
			axs[axsIndex].set_ylabel('GC%')
			axs[axsIndex].tick_params(axis='x', labelsize=20*fontScale)
			axs[axsIndex].set_xlabel('Read Length', size=20*fontScale, wrap = True)
			axs[axsIndex].set_title(fileDict[key]["filePrefix"], size = 20*fontScale, wrap = True)
		else:
			axs[axsIndex].set_xlabel('Read Length', size=20*fontScale, wrap = True)
			axs[axsIndex].tick_params(axis='x', labelsize=20*fontScale)
			axs[axsIndex].set_title(fileDict[key]["filePrefix"], fontsize = 20*fontScale, wrap = True)
		axsIndex += 1
	cbar = plt.colorbar(h[3], ax=axs[axsIndex-1], pad = 0.01, aspect = 50)
	cbar.set_label('# of reads', rotation=90, fontsize = 10)
	cbar.ax.tick_params(labelsize=5) 
	plt.savefig("Length_vs_GC.pdf" , format = "pdf")
	plt.close()
	
	
	# 2D Histogram of log-length vs. GC% for both runs
	fig, axs = plt.subplots(1, numFiles, sharex=True, sharey=True, tight_layout=False)
	axsIndex = 0
	for key in sorted(fileDict.keys()):
		h = axs[axsIndex].hist2d(fileDict[key]["logLengthList"], fileDict[key]["gcList"], bins = 100, norm = colors.LogNorm())
		axsIndex += 1
	axsIndex = 0
	for key in sorted(fileDict.keys()):
		if axsIndex == 0:
			axs[axsIndex].set_ylabel('GC%')
			axs[axsIndex].tick_params(axis='x', labelsize=20*fontScale)
			axs[axsIndex].set_xlabel('Log-transformed Read Length', size=20*fontScale, wrap = True)
			axs[axsIndex].set_title(fileDict[key]["filePrefix"], size = 20*fontScale, wrap = True)
		else:
			axs[axsIndex].set_xlabel('Log-transformed Read Length', size=20*fontScale, wrap = True)
			axs[axsIndex].tick_params(axis='x', labelsize=20*fontScale)
			axs[axsIndex].set_title(fileDict[key]["filePrefix"], fontsize = 20*fontScale, wrap = True)
		axsIndex += 1
	cbar = plt.colorbar(h[3], ax=axs[axsIndex-1], pad = 0.01, aspect = 50)
	cbar.set_label('# of reads', rotation=90, fontsize = 10)
	cbar.ax.tick_params(labelsize=5) 
	plt.savefig("LogLength_vs_GC.pdf" , format = "pdf")
	plt.close()
	
	
	# 2D Histograms of quality vs. GC% for both runs
	fig, axs = plt.subplots(1, numFiles, sharex=True, sharey=True, tight_layout=False)
	axsIndex = 0
	for key in sorted(fileDict.keys()):
		h = axs[axsIndex].hist2d(fileDict[key]["qualList"], fileDict[key]["gcList"], bins = 100, norm = colors.LogNorm())
		axsIndex += 1
	axsIndex = 0
	for key in sorted(fileDict.keys()):
		if axsIndex == 0:
			axs[axsIndex].set_ylabel('GC%')
			axs[axsIndex].tick_params(axis='x', labelsize=20*fontScale)
			axs[axsIndex].set_xlabel('Read Avg. Quality Score', size=20*fontScale, wrap = True)
			axs[axsIndex].set_title(fileDict[key]["filePrefix"], size = 20*fontScale, wrap = True)
		else:
			axs[axsIndex].set_xlabel('Read Avg. Quality Score', size=20*fontScale, wrap = True)
			axs[axsIndex].tick_params(axis='x', labelsize=20*fontScale)
			axs[axsIndex].set_title(fileDict[key]["filePrefix"], fontsize = 20*fontScale, wrap = True)
		axsIndex += 1
	cbar = plt.colorbar(h[3], ax=axs[axsIndex-1], pad = 0.01, aspect = 50)
	cbar.set_label('# of reads', rotation=90, fontsize = 10)
	cbar.ax.tick_params(labelsize=5) 
	plt.savefig("Quality_vs_GC.pdf" , format = "pdf")
	plt.close()
	
	print(datetime.now())
	
elif numFiles == 1:
	# Histogram of lengths
	fig, axs = plt.subplots(1, numFiles, sharex=True, sharey=True, tight_layout=False)
	for key in sorted(fileDict.keys()):
		axs.hist(fileDict[key]["lengthList"], bins = 100)
	for key in sorted(fileDict.keys()):
		axs.set_ylabel('No. of Reads')
		axs.set_xlabel('Read Length', size=10, wrap=True)
		axs.tick_params(axis='x', labelsize=10)
		axs.set_title(fileDict[key]["filePrefix"], size=10, wrap=True)
	plt.savefig("Sequence_lengths.pdf" , format = "pdf")
	plt.close()


	# Histograms of log-lengths for both runs
	fig, axs = plt.subplots(1, numFiles, sharex=True, sharey=True, tight_layout=False)
	for key in sorted(fileDict.keys()):
		axs.hist(fileDict[key]["logLengthList"], bins = 100)
	for key in sorted(fileDict.keys()):
		axs.set_ylabel('No. of Reads')
		axs.tick_params(axis='x', labelsize=10)
		axs.set_xlabel('Log-transformed Read Length', size=10, wrap = True)
		axs.set_title(fileDict[key]["filePrefix"], size = 10, wrap = True)
	plt.savefig("Log_sequence_lengths.pdf" , format = "pdf")
	plt.close()

	# Histograms of quality for both runs
	fig, axs = plt.subplots(1, numFiles, sharex=True, sharey=True, tight_layout=False)
	for key in sorted(fileDict.keys()):
		axs.hist(fileDict[key]["qualList"], bins = 100)
	for key in sorted(fileDict.keys()):
		axs.set_ylabel('No. of Reads')
		axs.tick_params(axis='x', labelsize=10)
		axs.set_xlabel('Read Avg. Quality Score', size=10, wrap = True)
		axs.set_title(fileDict[key]["filePrefix"], size = 10, wrap = True)
	plt.savefig("Sequence_qualities.pdf" , format = "pdf")
	plt.close()


	# Histograms of GC% for both runs
	fig, axs = plt.subplots(1, numFiles, sharex=True, sharey=True, tight_layout=False)
	for key in sorted(fileDict.keys()):
		axs.hist(fileDict[key]["gcList"], bins = 100)
	for key in sorted(fileDict.keys()):
		axs.set_ylabel('No. of Reads')
		axs.tick_params(axis='x', labelsize=10)
		axs.set_xlabel('GC%', size=10, wrap = True)
		axs.set_title(fileDict[key]["filePrefix"], size = 10, wrap = True)
	plt.savefig("GC_content.pdf" , format = "pdf")
	plt.close()


	# 2D Histograms of length vs. quality for both runs
	fig, axs = plt.subplots(1, numFiles, sharex=True, sharey=True, tight_layout=False)
	for key in sorted(fileDict.keys()):
		h = axs.hist2d(fileDict[key]["lengthList"], fileDict[key]["qualList"], bins = 100, norm = colors.LogNorm())
	for key in sorted(fileDict.keys()):
		axs.set_ylabel('Read Avg. Quality Score')
		axs.tick_params(axis='x', labelsize=10)
		axs.set_xlabel('Read Length', size=10, wrap = True)
		axs.set_title(fileDict[key]["filePrefix"], size = 10, wrap = True)
	cbar = plt.colorbar(h[3], ax=axs, pad = 0.01, aspect = 50)
	cbar.set_label('# of reads', rotation=90, fontsize = 10)
	cbar.ax.tick_params(labelsize=5) 
	plt.savefig("Length_vs_Quality.pdf" , format = "pdf")
	plt.close()
	
	
	# 2D Histogram of log-length vs. quality for both runs
	fig, axs = plt.subplots(1, numFiles, sharex=True, sharey=True, tight_layout=False)
	for key in sorted(fileDict.keys()):
		h = axs.hist2d(fileDict[key]["logLengthList"], fileDict[key]["qualList"], bins = 100, norm = colors.LogNorm())
	for key in sorted(fileDict.keys()):
		axs.set_ylabel('Read Avg. Quality Score')
		axs.tick_params(axis='x', labelsize=10)
		axs.set_xlabel('Log-transformed Read Length', size=10, wrap = True)
		axs.set_title(fileDict[key]["filePrefix"], size = 10, wrap = True)
	cbar = plt.colorbar(h[3], ax=axs, pad = 0.01, aspect = 50)
	cbar.set_label('# of reads', rotation=90, fontsize = 10)
	cbar.ax.tick_params(labelsize=5) 
	plt.savefig("LogLength_vs_Quality.pdf" , format = "pdf")
	plt.close()
	
	
	# 2D Histograms of length vs. GC% for both runs
	fig, axs = plt.subplots(1, numFiles, sharex=True, sharey=True, tight_layout=False)
	for key in sorted(fileDict.keys()):
		h = axs.hist2d(fileDict[key]["lengthList"], fileDict[key]["gcList"], bins = 100, norm = colors.LogNorm())
	for key in sorted(fileDict.keys()):
		axs.set_ylabel('GC%')
		axs.tick_params(axis='x', labelsize=10)
		axs.set_xlabel('Read Length', size=10, wrap = True)
		axs.set_title(fileDict[key]["filePrefix"], size = 10, wrap = True)
	cbar = plt.colorbar(h[3], ax=axs, pad = 0.01, aspect = 50)
	cbar.set_label('# of reads', rotation=90, fontsize = 10)
	cbar.ax.tick_params(labelsize=5) 
	plt.savefig("Length_vs_GC.pdf" , format = "pdf")
	plt.close()
	
	
	# 2D Histogram of log-length vs. GC% for both runs
	fig, axs = plt.subplots(1, numFiles, sharex=True, sharey=True, tight_layout=False)
	for key in sorted(fileDict.keys()):
		h = axs.hist2d(fileDict[key]["logLengthList"], fileDict[key]["gcList"], bins = 100, norm = colors.LogNorm())
	for key in sorted(fileDict.keys()):
		axs.set_ylabel('GC%')
		axs.tick_params(axis='x', labelsize=10)
		axs.set_xlabel('Log-transformed Read Length', size=10, wrap = True)
		axs.set_title(fileDict[key]["filePrefix"], size = 10, wrap = True)
	cbar = plt.colorbar(h[3], ax=axs, pad = 0.01, aspect = 50)
	cbar.set_label('# of reads', rotation=90, fontsize = 10)
	cbar.ax.tick_params(labelsize=5) 
	plt.savefig("LogLength_vs_GC.pdf" , format = "pdf")
	plt.close()
	
	
	# 2D Histograms of quality vs. GC% for both runs
	fig, axs = plt.subplots(1, numFiles, sharex=True, sharey=True, tight_layout=False)
	for key in sorted(fileDict.keys()):
		h = axs.hist2d(fileDict[key]["qualList"], fileDict[key]["gcList"], bins = 100, norm = colors.LogNorm())
	for key in sorted(fileDict.keys()):
		axs.set_ylabel('GC%')
		axs.tick_params(axis='x', labelsize=10)
		axs.set_xlabel('Read Avg. Quality Score', size=10, wrap = True)
		axs.set_title(fileDict[key]["filePrefix"], size = 10, wrap = True)
	cbar = plt.colorbar(h[3], ax=axs, pad = 0.01, aspect = 50)
	cbar.set_label('# of reads', rotation=90, fontsize = 10)
	cbar.ax.tick_params(labelsize=5) 
	plt.savefig("Quality_vs_GC.pdf" , format = "pdf")
	plt.close()
	
	print(datetime.now())