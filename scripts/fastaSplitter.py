#! /usr/bin/python
'''
'''
import sys

file=sys.argv[1]
filename = file.split(".")
splitLength=int(sys.argv[2])

fileCounter = 1
seqCounter = 0

print "Splitting fasta to files with %d sequences" %(splitLength)

openfile = open(file, "r")

outfile = open("%s_%d.fasta" %(filename[0],fileCounter), "w")
print "Writing file %d..." %(fileCounter)
for line in openfile:
	if ">" in line:
		seqCounter += 1
		if seqCounter > splitLength:
			fileCounter += 1
			outfile.close()
			outfile = open("%s_%d.fasta" %(filename[0],fileCounter), "w")
			print "Writing file %d..." %(fileCounter)
			seqCounter = 1
	outfile.write(line)
outfile.close()


filename = filename[0]
trimName = filename.split(".")

for i in range(1,fileCounter):
	trimScript = open("mafft_%s_%d.sh" %(trimName[0],i), "w")
	trimScript.write(
'''#!/bin/bash -l\n
#SBATCH --job-name=mafft    # Job name
#SBATCH --ntasks=1                    # Run on a single CPU
#SBATCH --cpus-per-task=16
#SBATCH --output=mafft_%s_%d.out   # Standard output and error log

pwd; hostname; date

module load mafft/7.309
mafft --thread 16 --adjustdirection --localpair --maxiterate 10 %s_%d.fasta > %s_%d_redirected.fasta

date
''' %(trimName[0], i, trimName[0], i, trimName[0], i))
	trimScript.close()