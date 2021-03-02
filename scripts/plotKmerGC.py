#! /home/ps997/miniconda3/bin/python
# Needs to run with python 3 on our server!!!

# Usage: plotKmerGC.py kmer_counts_table.tsv

# 2020 February 27
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

kmerCounts = open(sys.argv[1], "r")

kmerList = []
kmerFreqList = []
gcList = []

for line in kmerCounts:
	splitline = line.strip("\n").split("\t")
	kmerList.append(splitline[0])
	kmerFreqList.append(int(splitline[1]))
	gcList.append(float(GC(splitline[0])))

fig, axs = plt.subplots(1, 1, sharex=True, sharey=True, tight_layout=False)
h = axs.hist2d(kmerFreqList, gcList, bins = 20, norm = colors.LogNorm())
axs.set_ylabel('GC%')
axs.tick_params(axis='x', labelsize=10)
axs.set_xlabel('Kmer Frequency', size=10, wrap = True)
cbar = plt.colorbar(h[3], ax=axs, pad = 0.01, aspect = 50)
cbar.set_label('Abundance', rotation=90, fontsize = 10)
cbar.ax.tick_params(labelsize=5) 
plt.savefig("%s.pdf" %(sys.argv[1]) , format = "pdf")
plt.close()

