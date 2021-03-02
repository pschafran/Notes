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

file1 = sys.argv[1]
file1Prefix = file1.split(".fastq")[0]
if "/" in file1Prefix:
	file1Prefix = file1Prefix.split("/")[-1]
file2 = sys.argv[2]
file2Prefix = file2.split(".fastq")[0]
if "/" in file2Prefix:
	file2Prefix = file2Prefix.split("/")[-1]

run1LengthList = []
run2LengthList = []
run1QualList = []
run2QualList = []
run1GCList = []
run2GCList = []

print("Parsing file 1...")
run1Index = SeqIO.index(file1, "fastq")
recordCounter = 0
for key in run1Index.keys():
	record = run1Index[key]
	run1LengthList.append(len(record.seq))
	run1QualList.append(sum(record.letter_annotations["phred_quality"])/len(record.letter_annotations["phred_quality"]))
	run1GCList.append(GC(record.seq))
	
	#Progress bar
	recordCounter += 1
	completionPerc = int(float(100*recordCounter)/float(len(run1Index)))
	sys.stdout.write('\r')
	sys.stdout.write("[%-100s] %d%%" % ('='*completionPerc, completionPerc))
	sys.stdout.flush()
print("\n")

print("Parsing file 2...")
run2Index = SeqIO.index(file2, "fastq") 
recordCounter = 0
for key in run2Index.keys():
	record = run2Index[key]
	run2LengthList.append(len(record.seq))
	run2QualList.append(sum(record.letter_annotations["phred_quality"])/len(record.letter_annotations["phred_quality"]))
	run2GCList.append(GC(record.seq))
	
	#Progress bar
	recordCounter += 1
	completionPerc = int(float(100*recordCounter)/float(len(run2Index)))
	sys.stdout.write('\r')
	sys.stdout.write("[%-100s] %d%%" % ('='*completionPerc, completionPerc))
	sys.stdout.flush()
print("\n")

# Log transform lengths
run1LogLength = [log(record) for record in run1LengthList]
run2LogLength = [log(record) for record in run2LengthList] 
	
### Plotting stuff ###
print("Plotting stuff...")
# Histograms of lengths for both runs
fig, axs = plt.subplots(1, 2, sharex=True, sharey=True, tight_layout=False)
axs[0].hist(run1LengthList, bins = 100)
axs[1].hist(run2LengthList, bins = 100)
axs[0].set(xlabel='Read Length', ylabel='No. of Reads', title=file1)
axs[1].set(xlabel='Read Length', title=file2)
plt.savefig("%s_%s_length.pdf" %(file1Prefix, file2Prefix), format = "pdf")
plt.close()

# Histograms of lengths <50k for both runs
fig, axs = plt.subplots(1, 2, sharex=True, sharey=True, tight_layout=False)
axs[0].hist(run1LengthList, bins = 100, range = [0,50000])
axs[1].hist(run2LengthList, bins = 100, range = [0,50000])
axs[0].set(xlabel='Read Length', ylabel='No. of Reads', title=file1)
axs[1].set(xlabel='Read Length', title=file2)
plt.savefig("%s_%s_length_max50k.pdf" %(file1Prefix, file2Prefix), format = "pdf")
plt.close()

# Histograms of log-lengths for both runs
fig, axs = plt.subplots(1, 2, sharex=True, sharey=True, tight_layout=False)
axs[0].hist(run1LogLength, bins = 100)
axs[1].hist(run2LogLength, bins = 100)
axs[0].set(xlabel='Log-transformed Read Length', ylabel='No. of Reads', title=file1)
axs[1].set(xlabel='Log-transformed Read Length', title=file2)
plt.savefig("%s_%s_loglength.pdf" %(file1Prefix, file2Prefix), format = "pdf")
plt.close()

# Histograms of quality for both runs
fig, axs = plt.subplots(1, 2, sharex=True, sharey=True, tight_layout=False)
axs[0].hist(run1QualList, bins = 100)
axs[1].hist(run2QualList, bins = 100)
axs[0].set(xlabel='Read Avg. Quality Score', ylabel='No. of Reads', title=file1)
axs[1].set(xlabel='Read Avg. Quality Score', title=file2)
plt.savefig("%s_%s_quality.pdf" %(file1Prefix, file2Prefix), format = "pdf")
plt.close()

# Histograms of GC% for both runs
fig, axs = plt.subplots(1, 2, sharex=True, sharey=True, tight_layout=False)
axs[0].hist(run1GCList, bins = 100)
axs[1].hist(run2GCList, bins = 100)
axs[0].set(xlabel='Read GC%', ylabel='No. of Reads', title=file1)
axs[1].set(xlabel='Read GC%', title=file2)
plt.savefig("%s_%s_GCperc.pdf" %(file1Prefix, file2Prefix), format = "pdf")
plt.close()

# 2D Histograms of length vs. quality for both runs

fig, axs = plt.subplots(1, 2, sharey=True, tight_layout=False)
axs[0].hist2d(run1LengthList, run1QualList, bins=100, norm = colors.LogNorm())
h = axs[1].hist2d(run2LengthList, run2QualList, bins=100, norm = colors.LogNorm())
axs[0].set(xlabel='Read Length', ylabel='Read Avg. Quality Score', title=file1)
axs[1].set(xlabel='Read Length', title=file2)
cbar = plt.colorbar(h[3], ax=axs[1], pad = 0.01, aspect = 50)
cbar.set_label('# of reads', rotation=90, fontsize = 10)
cbar.ax.tick_params(labelsize=5) 
plt.savefig("%s_%s_length_vs_qual.pdf" %(file1Prefix, file2Prefix), format = "pdf")
plt.close()

# 2D Histogram of log-length vs. quality for both runs
fig, axs = plt.subplots(1, 2, sharex=True, sharey=True, tight_layout=False)
axs[0].hist2d(run1LogLength, run1QualList, bins=100)
h = axs[1].hist2d(run2LogLength, run2QualList, bins=100)
axs[0].set(xlabel='Log-transformed Read Length', ylabel='Read Avg. Quality Score', title=file1)
axs[1].set(xlabel='Log-transformed Read Length', title=file2)
cbar = plt.colorbar(h[3], ax=axs[1], pad = 0.01, aspect = 50)
cbar.set_label('# of reads', rotation=90, fontsize = 10)
cbar.ax.tick_params(labelsize=5) 
plt.savefig("%s_%s_loglength_vs_qual.pdf" %(file1Prefix, file2Prefix), format = "pdf")
plt.close()

# 2D Histograms of length vs. GC% for both runs
fig, axs = plt.subplots(1, 2, sharey=True, tight_layout=False)
axs[0].hist2d(run1LengthList, run1GCList, bins=100, norm = colors.LogNorm())
h = axs[1].hist2d(run2LengthList, run2GCList, bins=100, norm = colors.LogNorm())
axs[0].set(xlabel='Read Length', ylabel='Read GC%', title=file1)
axs[1].set(xlabel='Read Length', title=file2)
cbar = plt.colorbar(h[3], ax=axs[1], pad = 0.01, aspect = 50)
cbar.set_label('# of reads', rotation=90, fontsize = 10)
cbar.ax.tick_params(labelsize=5) 
plt.savefig("%s_%s_length_vs_GCperc.pdf" %(file1Prefix, file2Prefix), format = "pdf")
plt.close()

# 2D Histogram of log-length vs. GC% for both runs
fig, axs = plt.subplots(1, 2, sharex=True, sharey=True, tight_layout=False)
axs[0].hist2d(run1LogLength, run1GCList, bins=100, norm = colors.LogNorm())
h = axs[1].hist2d(run2LogLength, run2GCList, bins=100, norm = colors.LogNorm())
axs[0].set(xlabel='Log-transformed Read Length', ylabel='Read GC%', title=file1)
axs[1].set(xlabel='Log-transformed Read Length', title=file2)
cbar = plt.colorbar(h[3], ax=axs[1], pad = 0.01, aspect = 50)
cbar.set_label('# of reads', rotation=90, fontsize = 10)
cbar.ax.tick_params(labelsize=5) 
plt.savefig("%s_%s_loglength_vs_GCperc.pdf" %(file1Prefix, file2Prefix), format = "pdf")
plt.close()

# 2D Histograms of quality vs. GC% for both runs
fig, axs = plt.subplots(1, 2, sharey=True, tight_layout=False)
axs[0].hist2d(run1QualList, run1GCList, bins=100, norm = colors.LogNorm())
h = axs[1].hist2d(run2QualList, run2GCList, bins=100, norm = colors.LogNorm())
axs[0].set(xlabel='Read Avg. Quality Score', ylabel='Read GC%', title=file1)
axs[1].set(xlabel='Read Avg. Quality Score', title=file2)
cbar = plt.colorbar(h[3], ax=axs[1], pad = 0.01, aspect = 50)
cbar.set_label('# of reads', rotation=90, fontsize = 10)
cbar.ax.tick_params(labelsize=5) 
plt.savefig("%s_%s_quality_vs_GCperc.pdf" %(file1Prefix, file2Prefix), format = "pdf")
plt.close()

