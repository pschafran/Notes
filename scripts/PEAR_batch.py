#! /usr/bin/python

import re, sys, subprocess


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
	if '_paired' in file:
		filedict[splitfile[0]].append(file)
mem=input("Input memory to use in gigabytes:")
thread=input("Input number of threads to use:")	
for key in filedict.keys():
	filename1=filedict[key][0]
	filename2=filedict[key][1]
	subprocess.call(['PEAR', '-f', '%s' %(filename1), '-r', '%s' %(filename2), '-y', '%dG' %(mem),'-j', '%d' %(thread),'-o', '%s_PEAR_TRIM_' %(key)])