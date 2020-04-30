#! /usr/bin/python

# specify source files, not .out files

import re
import sys
import subprocess

files = sys.argv[1:]
outfile=open("missing_accessions.txt", "w")

for file in files:
	opensourcefile = open(file, "r")
	try:
		openhitfile = open("%s.out" %(file), "r")
	except:
		print("No output file for %s! Running Entrez search..." %(file))
		openhitfile = open("%s.out" %(file), "w")
		entrez_call = subprocess.Popen([ "cat %s | epost -db protein | esummary -db protein | xtract -pattern DocumentSummary -element Caption,TaxId > %s.out" %(file, file)], stdout=openhitfile, shell=True)		
		entrez_call.wait()
		openhitfile.close()	
		openhitfile = open("%s.out" %(file), "r")
	hitList = []
	for line in openhitfile:
		hitList.append(line.strip("\n").split("\t")[0])
#	print(len(hitList))
	for line2 in opensourcefile:
#		print line2.strip("\n").split(".")
		if line2.strip("\n").split(".")[0] not in hitList:
			outfile.write(line2)
		elif line2.strip("\n").split(".")[0] in hitList:
			pass
#			print("Found match for %s" %(line2))
	opensourcefile.close()
	openhitfile.close()
outfile.close()

	
