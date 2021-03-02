#! /usr/bin/python
'''Command Line:getScaffoldsFromFasta.py scaffoldFile.fasta ScaffoldID(.txt)'''
import sys
import os.path

scaffoldFile = sys.argv[1]
scaffoldID = sys.argv[2]

openInFile = open(scaffoldFile, "r")
openOutFile = open("%s.fasta" %(scaffoldID), "w")

if os.path.isfile(scaffoldID):
	infile = open(scaffoldID, "r")
	for scaffoldName in infile:
		writeOut = 0
		for line in openInFile:
			if ">" in line and scaffoldName.strip(">\n").split(" ")[0] == line.strip(">\n").split(" ")[0]:
				print "Found scaffold %s" %(line)
				writeOut = 1
				openOutFile.write(line)
			elif ">" in line and scaffoldName != line.strip(">\n"):
				writeOut = 0
			elif ">" not in line:
				if writeOut == 1:
					openOutFile.write(line)
		openInFile.seek(0)
	infile.close()
else:
	writeOut = 0
	for line in openInFile:
		if ">" in line and scaffoldID == line.strip(">\n"):
			print "Found scaffold %s" %(line)
			writeOut = 1
			openOutFile.write(line)
		elif ">" in line and scaffoldID != line.strip(">\n"):
			writeOut = 0
		elif ">" not in line:
			if writeOut == 1:
				openOutFile.write(line)
		
openInFile.close()
openOutFile.close()
