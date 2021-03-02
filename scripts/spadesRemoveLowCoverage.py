#! /usr/bin/python
#This reports # scaffolds, total length, and average length of scaffolds.fasta produced by SPAdes
#Usage: ./spadesStats.py path(s)/to/files.fasta
#
import sys
import numpy as np


for file in sys.argv[1:]:
	infile = open(file, 'r')
	cov_list=[]
	for line in infile:
		if ">" in line:
			splitline=line.split('_')
			cov=float(splitline[5])
			cov_list.append(cov)
	stdDev = np.std(cov_list, dtype=np.float64)
	maxCov = max(cov_list)
	infile.seek(0)
	outfile = open("%s_hiCov.fasta" %(file), "w")
	writeOut=0
	for line in infile:
		if ">" in line:
			splitline=line.split('_')
			cov=float(splitline[5])
			if cov >= (maxCov - (stdDev*3.0)):
				writeOut=1
				outfile.write(line)
			elif cov < (maxCov - (stdDev*3.0)):
				writeOut=0
		elif ">" not in line and writeOut==1:
			outfile.write(line)
		elif ">" not in line and writeOut==0:
			pass
	outfile.close()
	
