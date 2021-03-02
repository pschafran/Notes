#! /usr/bin/python
import sys, subprocess
filelist = sys.argv[1:]
prefix=raw_input("New file prefix:")
filecount=1
for file in filelist:
	fileext = file.split(".")
	subprocess.call(['mv', './%s' %(file), './%s_%s.%s' %(prefix,filecount,fileext[1])])
	filecount+=1