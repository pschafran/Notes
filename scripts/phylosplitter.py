#! /usr/bin/python
##2016 Nov 1
##Written by: Peter Schafran

##Script to slide along genome alignments (FASTA format) and split into separate files of n length for MrBayes analysis (Nexus formaat).

'''Usage: phylosplitter.py alignment.fasta split_length'''

import sys, re, math

alignment_file=sys.argv[1]
split_length=sys.argv[2]
alignment_dict={}
shells=1

#Parse FASTA file, create dictionary with sequence name key, respective sequence as entry
openfile = open(alignment_file, 'r')
for line in openfile:
	if '>' in line:
		dictkey=line.strip('\n')
		alignment_dict[dictkey]=0
	if not '>' in line:
		dictseq=line.strip('\n')
		alignment_dict[dictkey]=dictseq
		alignment_total_length=len(line)
openfile.close()

#Calculate number of segments alignment will be broken into, rounded down to nearest integer
split_length=int(split_length)
num_divisions = alignment_total_length/split_length
num_divisions = int(num_divisions)

#Set variables to be used in writing output files
ntax = len(alignment_dict.keys())
division_count = 1
startplace = 1
endplace = split_length

#Create individual files for each split of alignment
while division_count <= num_divisions:
	filename = "%s_%s-%s.nex" %(alignment_file, startplace, endplace)
	raxmlfile = "%s_%s-%s.phy" %(alignment_file, startplace, endplace)
#write MrBayes file
	outfile = open(filename, 'w')
	outfile.write("#NEXUS\n")
	outfile.write("begin data;\n")
	outfile.write("\tdimensions ntax=%d nchar=%d;\n" %(ntax, split_length))
	outfile.write("\tformat datatype=dna missing=? gap=-;\n")
	outfile.write("\tmatrix\n")
#write RAxML file
	outfile2 = open(raxmlfile, 'w')
	outfile2.write(" %d %d\n" %(ntax, split_length))
	division_count += 1
	startplace = startplace + split_length
	endplace = endplace + split_length
	outfile.close()
	outfile2.close()

#Split alignment dictionary into separate dictionaries with alignments split by split_length
split_dict={}
for key in alignment_dict.keys():
	division_count = 1
	base_count = 1
	startplace = 1
	endplace = split_length
	split_dict[key]={}
	while division_count <= num_divisions:
		split_dict[key][division_count]=[]
		division_count += 1
#	print key
	for key2 in split_dict[key]:
#		print key2
		base_count = 1
		for base in alignment_dict[key]:
			if base_count >= startplace and base_count <= endplace:
				split_dict[key][key2].append(base)
#				print base_count
#				print endplace
				base_count += 1
			elif base_count < startplace or base_count > endplace:
				base_count += 1
		startplace = startplace + split_length
		endplace = endplace + split_length
	
	
#for key in split_dict.keys():
#	print split_dict[key]
	
#Write homologous segments of alignment to files
startplace = 1
endplace = split_length
division_count = 1

while division_count <= num_divisions:
	filename = "%s_%s-%s.nex" %(alignment_file, startplace, endplace)
	raxmlfile = "%s_%s-%s.phy" %(alignment_file, startplace, endplace)
	outfile = open(filename, 'a')
	outfile2 = open(raxmlfile, 'a')
	for key in split_dict.keys():
		for key2 in split_dict[key]:
			if division_count == key2:
				entry = ''.join(split_dict[key][key2])
				stripkey = key.strip(">") 
				outfile.write("\t%s\t%s\n" %(stripkey, entry))
				outfile2.write("%s %s\n" %(stripkey, entry))
	outfile.write(";\n")
	outfile.write("end;")
	outfile.write("begin mrbayes;\n")
	outfile.write("\tlog start filename=log_%s;\n" %(filename))
	outfile.write("\tset autoclose=yes nowarn=yes;\n")
	outfile.write("\tlset nst=2 rates=inv;\n")
	outfile.write("\tmcmc ngen=50000000 samplefreq=1000 printfreq=10000 nchains=4;\n")
	outfile.write("\tsump burnin=12500;\n")
	outfile.write("\tsumt burnin=12500;\n")
	outfile.write("\tlog stop;\n")
	outfile.write("end;\n")
	outfile.close()
	outfile2.close()
	division_count += 1
	startplace = startplace + split_length
	endplace = endplace + split_length
	

#Write shell scripts for submitting all to computing cluster
if shells == 1:
	startplace = 1
	endplace = split_length
	division_count = 1
	while division_count <= num_divisions:
		shellname = "%s_%s-%s_MrBayes.sh" %(alignment_file, startplace, endplace)
		raxmlshell = "%s_%s-%s_RAxML.sh" %(alignment_file, startplace, endplace)
		outfile = open(shellname, 'w')
		outfile2 = open(raxmlshell, 'w')
		outfile.write("#!/bin/bash -l\n")
		outfile2.write("#!/bin/bash -l\n")
		outfile.write("#$ -o log_%s.txt -j y\n" %(shellname))
		outfile2.write("#SBATCH --job-name=RAxML    # Job name\n")
		outfile.write("#$ -l excl=true\n")
		outfile2.write("#SBATCH --ntasks=1                    # Run on a single CPU\n")
		outfile.write("#$ -q main\n")
		outfile2.write("#SBATCH --output=%s.out   # Standard output and error log\n" %(raxmlshell))
		outfile.write("#$ -m beas\n")
		outfile.write("#$ -M pscha005@odu.edu\n")
		outfile2.write("\n")
		outfile.write("\n")
		outfile2.write("pwd; hostname; date\n")
		outfile.write("module load mrbayes/3.2.6\n")
		outfile.write("mb %s_%s-%s.nex\n" %(alignment_file, startplace, endplace))
		outfile2.write("module load /cm/shared/modulefiles/openmpi/gcc/64/1.6.5\n")
		outfile2.write("/cm/shared/apps/raxml/7.3.0/raxmlHPC -s %s_%s-%s.phy -n output_%s_%s-%s -m GTRGAMMA -f a -x 1 -N 1000 -d -p 3556\n" %(alignment_file, startplace, endplace, alignment_file, startplace, endplace))
		outfile.close()
		outfile2.close()
		division_count += 1
		startplace = startplace + split_length
		endplace = endplace + split_length
	


