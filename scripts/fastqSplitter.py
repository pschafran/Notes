#! /usr/bin/python
'''
'''
import sys

file=sys.argv[1]
filename = file.split(".")
splitLength=int(sys.argv[2])

fileCounter = 1
seqCounter = 0
lineCounter = 0
print "Splitting fastq to files with %d sequences" %(splitLength)

openfile = open(file, "r")

outfile = open("%s_%d.fastq" %(filename[0],fileCounter), "w")
print "Writing file %d..." %(fileCounter)
for line in openfile:
	if line.startswith("@") and lineCounter == 0 or line.startswith("@") and lineCounter % 4 == 0:
		seqCounter += 1
		if seqCounter > splitLength:
			fileCounter += 1
			outfile.close()
			outfile = open("%s_%d.fastq" %(filename[0],fileCounter), "w")
			print "Writing file %d..." %(fileCounter)
			seqCounter = 1
	outfile.write(line)
	lineCounter += 1
outfile.close()
