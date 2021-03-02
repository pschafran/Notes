#! /usr/bin/python
'''
Separates a multi-fasta file and creates and individual fasta file for each sequence. Also creates a SLURM submission script to run minimap
'''
import sys

file=sys.argv[1]
filename = file.split(".")
openfile = open(file, "r")
for line in openfile:
	if ">" in line:
		try:
			outputfile.close()
		except:
			print(" ")
		finally:
			splitline = line.strip("\n").split(">")
			newFilename = splitline[1]
			outputfile = open("%s.fasta" %(newFilename), "w")
			outputfile.write(">%s\n" %(newFilename))
			mapScript = open("%s.sh" %(newFilename), "w")
			mapScript.write(
'''#!/bin/bash -l\n
#SBATCH --job-name=mapEng    # Job name
#SBATCH --mail-type=FAIL               # Mail events (FAIL)
#SBATCH --mail-user=<pscha005@odu.edu>   # Where to send mail
#SBATCH --ntasks=1                    # Run on a single CPU
#SBATCH --cpus-per-task=16

enable_lmod
module load minimap2
minimap2 -ax asm5n /scratch-lustre/pscha005/Genomes/taiwanensis/%s.fasta /scratch-lustre/pscha005/IridianGenomes/engelmannii/scaffolds/Iengl_SRR8371589_30Mread_scaffolds_unique.fasta > Iengl_to_%s.sam

''' %(newFilename, newFilename))
			mapScript.close()
	if ">" not in line:
		outputfile.write(line)
outputfile.close()

