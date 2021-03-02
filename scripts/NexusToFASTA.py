#! /usr/bin/python

import fileinput, sys, os

filelist = sys.argv[1:]

for file in filelist:
	print "Processing: %s" % file
	infile = open(file, 'r')
	outfile = open('%s.fasta' %(file), 'w')
	linecount=1
	for line in infile:
		if linecount >= 9:
			linesplit = line.split(' ')
			#print linesplit[0]
			#print linesplit[1]
			if len(linesplit) < 2:
				#print "BREAK!"
				break
			else:
				if linecount == 9:
					outfile.write('>%s\n%s' %(linesplit[0],linesplit[1]))
				else:
					outfile.write(linesplit[1])
			#print linesplit
		else:
			pass
		linecount += 1
	outfile.close()
	infile.close()

