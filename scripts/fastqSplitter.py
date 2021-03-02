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


filename = filename[0]
trimName = filename.split("_")

for i in range(1,fileCounter):
	trimScript = open("trimmomatic_%s_%d.sh" %(trimName[0],i), "w")
	trimScript.write(
'''#!/bin/bash -l\n
#SBATCH --job-name=trimEng    # Job name
#SBATCH --mail-type=ALL               # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=<pscha005@odu.edu>   # Where to send mail
#SBATCH --ntasks=1                    # Run on a single CPU
#SBATCH --cpus-per-task=16
#SBATCH --output=trimmomatic_%s_%d.out   # Standard output and error log

pwd; hostname; date

enable_lmod
module load java/11
java -jar /scratch-lustre/pscha005/trimmomatic/trimmomatic-0.33.jar PE -threads 16 %s_1_%d.fastq %s_2_%d.fastq %s_1_%d_paired.fq.gz %s_1_%d_unpaired.fq.gz %s_2_%d_paired.fq.gz %s_2_%d_unpaired.fq.gz ILLUMINACLIP:/scratch-lustre/pscha005/trimmomatic/adapters/TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
date
''' %(trimName[0], i, trimName[0], i, trimName[0], i, trimName[0], i, trimName[0], i, trimName[0], i, trimName[0], i))
	trimScript.close()
