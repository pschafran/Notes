#! /usr/bin/python

import sys

filelist = sys.argv[1:]

seqNameList = []

for file in filelist:
	infile = open(file, "r")
	filename = ".".join(file.split(".")[:-1])
	fileext = file.split(".")[-1]
	outfile = open("%s.renamed.%s" %(filename, fileext), "w")
	for line in infile:
		if line.startswith(">"):
			counter = 0
			while line.strip(">\n") in seqNameList:
				counter += 1
				line = "%s-%s" %(line.strip(">\n"), counter)
			seqNameList.append(line.strip(">\n"))
			outfile.write(">%s\n" % line.strip(">\n"))
		else:
			outfile.write(line)
	outfile.close()
	infile.close()
