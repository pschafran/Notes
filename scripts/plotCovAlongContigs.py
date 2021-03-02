#! /home/ps997/miniconda3/bin/python
# Needs to run with python 3 on our server!!!

# Usage: plotCovAlongContigs.py depth_output_from_samtools.txt SmoothingFactor (% of data to average across, must be between 1-100. Suggested setting from 1-15)

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import sys
from matplotlib import colors
from numpy import log
import matplotlib.patches as patches
from mpl_toolkits.mplot3d import Axes3D

def running_mean(x, N):
	cumsum = np.cumsum(np.insert(x, 0, 0)) 
	return (cumsum[N:] - cumsum[:-N]) / float(N)

file = open(sys.argv[1],"r")
smoothingFactor = float(sys.argv[2])
if smoothingFactor < 1:
	print("WARNING: Smoothing was set to less than 1. Now set to 1.")
	smoothingFactor = 1
if smoothingFactor > 100:
	print("WARNING: Smooting was set to greater than 100. Now set to 100.")
	smoothingFactor = 100

contigList = []
posList = []
depthList = []
contigDict = {}
assemblySize = 0
contigLengthList = []
indexer = 0
print('Parsing depth file...')
for line in file:
	splitline = line.strip("\n").split("\t")
	contig = splitline[0]
	pos = splitline[1]
	depth = splitline[2]
	try:
		contigDict[contig][0].append(int(pos))
		contigDict[contig][1].append(float(depth))
	except:
		contigDict[contig] = [[],[]]
		contigDict[contig][0].append(int(pos))
		contigDict[contig][1].append(float(depth))
print('Finding high coverage outlier regions...')
for key in contigDict.keys():
	assemblySize += len(contigDict[key][1])
	contigLengthList.append(len(contigDict[key][1]))
	contigDict[key].append(np.median(contigDict[key][1])) # contigDict[key][2]
	contigDict[key].append(np.std(contigDict[key][1])) # contigDict[key][3]
	contigDict[key].append(np.min(contigDict[key][1])) # contigDict[key][4]
	contigDict[key].append(np.max(contigDict[key][1])) # contigDict[key[5]
	if any(x > contigDict[key][2]+(3*contigDict[key][3]) for x in contigDict[key][1]):
		print('	%s' %(key))
		outfile = open("%s_hiCovRegions.txt" %(key), "w")
		outfile.write("position\tcoverage\n")
		for i in range(1, len(contigDict[key][0])):
			if contigDict[key][1][i-1] > contigDict[key][2]+(3*contigDict[key][3]):
				outfile.write("%d\t%f\n" %(contigDict[key][0][i-1], contigDict[key][1][i-1]))
		outfile.close()

n50size = 0
n50list = []
for length in sorted(contigLengthList, reverse=True):
	if n50size < assemblySize*0.5:
		n50list.append(length)
		n50size += length
	else:
		n50list.append(length)
		n50size += length
n50 = np.min(n50list)
n50 = 0
print("Assembly Size: %s" %(assemblySize))

movAvgDict = {}
print('Smoothing data...')
for key in contigDict.keys():
	smoothingBases = int(len(contigDict[key][1]) * (smoothingFactor/100))
	movAvgDict[key] = [np.array(running_mean(contigDict[key][0], smoothingBases)).tolist(), np.array(running_mean(contigDict[key][1], smoothingBases)).tolist()]

print('Plotting charts...')
for key in movAvgDict.keys():
	fig, ax = plt.subplots(1,1)
	ax.plot(movAvgDict[key][0], movAvgDict[key][1])
	ax.set_xlabel("Position")
	ax.set_ylabel("Read Depth")
	ax.set_title("%s" %(key))
	plt.grid(which='major',axis='both',linestyle='dashed')
	#ax.set_ylim(ymin=(contigDict[key][4] - 100), ymax = (contigDict[key][5] + 100))
	ax.set_ylim(ymin=0, ymax=10000)
	ax.set_axisbelow(True)
	ax.text(0.7, 1.1, "Median Depth: %d" %(contigDict[key][2]), horizontalalignment='left', verticalalignment='center', transform=ax.transAxes, bbox=dict(facecolor='white', edgecolor='white', alpha=0.5))
	ax.text(0.7, 1.05, "Standard Deviation: %d" %(contigDict[key][3]), horizontalalignment='left', verticalalignment='center', transform=ax.transAxes, bbox=dict(facecolor='white', edgecolor='white', alpha=0.5))
	plt.savefig("%s_%s.pdf" %(key,smoothingFactor),format = "pdf")
	plt.close()



fig, ax = plt.subplots(len(n50list), 1, sharex=True, sharey=True, tight_layout=True, figsize = (8,len(n50list)))
index = 0
for key in movAvgDict.keys():
	if len(contigDict[key][0]) >= n50:
		ax[index].plot(movAvgDict[key][0], movAvgDict[key][1])
		ax[index].text(1.01, 0.5, key, horizontalalignment='left', verticalalignment='center',transform=ax[index].transAxes, bbox=dict(facecolor='white', edgecolor='white', alpha=0.5))
		ax[index].spines['right'].set_visible(False)
		ax[index].spines['top'].set_visible(False)
		ax[index].set_ylabel("Read Depth")
		index += 1
ax[index-1].set_xlabel("Position")
plt.savefig("%s_n50contigs_%s.pdf" %((sys.argv[1].split(".txt")[0]), smoothingFactor),format = "pdf")
plt.close()

fig, ax = plt.subplots(len(n50list), 1, sharex=True, sharey=True, tight_layout=True, figsize = (8,len(n50list)))
index = 0
for key in movAvgDict.keys():
	halfway = int(len(n50list)/2)
	if len(contigDict[key][0]) >= n50:
		ax[index].plot(movAvgDict[key][0], movAvgDict[key][1])
		ax[index].text(1.01, 0.5, key, horizontalalignment='left', verticalalignment='center',transform=ax[index].transAxes, bbox=dict(facecolor='white', edgecolor='white', alpha=0.5))
		ax[index].spines['right'].set_visible(False)
		ax[index].spines['top'].set_visible(False)
		ax[index].set_yscale('log')
		if index == halfway:
			ax[index].set_ylabel("Read Depth")
		index += 1
ax[index-1].set_xlabel("Position")
plt.savefig("%s_n50contigs_logCov_%s.pdf" %((sys.argv[1].split(".txt")[0]), smoothingFactor),format = "pdf")
plt.close()



fig = plt.figure()
ax = fig.gca(projection='3d')
index = 0
for key in movAvgDict.keys():
	if len(contigDict[key][0]) >= n50:
		ax.plot(movAvgDict[key][0], movAvgDict[key][1], zs=index, zdir='y', label= "%s" %(key))
		index -= 1
plt.savefig("%s_3Dn50contigs_%s.pdf" %(sys.argv[1].split(".txt")[0], smoothingFactor),format = "pdf")
plt.close()







