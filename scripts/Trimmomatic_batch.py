#! /usr/bin/python

import sys, os, re, subprocess

'''Usage:
All files must be in directory with trimmomatic-0.3.3.jar. File names must be formatted: UniqueID1_forward.fastq, UniqueID1_reverse.fastq....Call input files in command line.
'''

filedict={}
filelist=sys.argv[1:]
filelistlen= len(filelist)
if filelistlen % 2 != 0:
	print "Uneven number of input files!"
	
#create dictionary keys based on unique ID in file name
for file in filelist:
	#print file
	splitfile = re.split('[_.]', file)
	filedict[splitfile[0]]=[]
#add forward and reverse file names to dictionary
for file in filelist:
	#print file
	splitfile = re.split('[_.]', file)
	filedict[splitfile[0]].append(file)

for key in filedict.keys():
	filename1=filedict[key][0]
	filename2=filedict[key][1]
	subprocess.call(['java', '-jar', 'trimmomatic-0.33.jar', 'PE', '-phred33', '%s' %(filename1), '%s' %(filename2), '%s_TRIM_forward_paired.fq.gz'%(key), '%s_TRIM_forward_unpaired.fq.gz'%(key), '%s_TRIM_reverse_paired.fq.gz'%(key), '%s_TRIM_reverse_unpaired.fq.gz'%(key), 'ILLUMINACLIP:./adapters/MiSeq.fa:2:30:10', 'LEADING:3', 'TRAILING:3', 'SLIDINGWINDOW:4:15', 'MINLEN:36'])