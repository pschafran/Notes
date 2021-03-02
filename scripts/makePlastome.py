#! /usr/bin/python

import re
import sys
import os
import subprocess

sampleDict = {} # Structure is {'Sample1':['Input_File_R1.fastq.gz','Input_File_R2.fastq.gz'], 'Sample2'...}
hpc = 0
fileList= []
if "-hpc" in sys.argv:
	hpc = 1
for item in sys.argv:
	if ".fastq" in item: 
		fileList.append(item)
	if ".fasta" in item:
		referenceFile = item
		referenceName = ".".join(referenceFile.split(".")[:-1])

try:
	if len(referenceFile) == 1:
		print "Using reference: %s" %(referenceFile)
except:
	print "*" *10
	print "WARNING: No reference genome provided. SPAdes will be run on all data."
	print "*" *10

if len(fileList) == 0:
	print "!" *10
	print "ERROR: No data (fastq files) supplied"
	print "!" *10
	sys.exit()

for fileName in fileList:
	splitFileName = fileName.strip("\n").split("_")
	sampleName = "_".join(splitFileName[:-2])
	sampleDict[sampleName] = []
	
for fileName in fileList:
	splitFileName = fileName.strip("\n").split("_")
	sampleName = "_".join(splitFileName[:-2])
	sampleDict[sampleName].append(fileName)

if hpc == 1:
	print "*"*10
	print "Running in HPC mode..."

'''Trimmomatic'''
def trimmomatic():
	print "*"*10
	print "Running trimmomatic..."
	print "*"*10
	global trimOutputDict
	trimOutputDict = {}
	for sample in sorted(sampleDict.keys()):
		subprocess.call(['mkdir', sample])
		if len(sampleDict[sample]) == 1 and hpc == 0:
			print "*" * 10
			print "WARNING: %s does not have 2 input files. Trimmomatic will be run in SE mode" %(sample)
			print "*" * 10
			try:
				subprocess.call(['java', '-jar', 'trimmomatic-0.33.jar', 'SE', '-phred33', '%s' %(sampleDict[sample][0]), './%s/%s_TRIM.fq.gz'%(sample,sample), 'ILLUMINACLIP:./adapters/TruSeq3-SE.fa:2:30:10', 'LEADING:3', 'TRAILING:3', 'SLIDINGWINDOW:4:15', 'MINLEN:36'])
			except:
				sys.exit("ERROR: Trimmomatic.jar not found in current directory")
			trimOutputDict[sample] = "%s_TRIM.fq.gz" %(sample)
		elif len(sampleDict[sample]) == 1 and hpc == 1:
			trimOutputDict[sample] = "%s_TRIM.fq.gz" %(sample)
			trimScript = open("trimmomatic_%s.sh" %(sample), "w")
			trimScript.write(
'''#!/bin/bash -l\n
#SBATCH --job-name=trim    # Job name
#SBATCH --ntasks=1                    # Run on a single CPU
#SBATCH --cpus-per-task=16
#SBATCH --output=trimmomatic_%s.out   # Standard output and error log

pwd; hostname; date

enable_lmod
module load java/11.0
java -jar /scratch-lustre/pscha005/trimmomatic/trimmomatic-0.33.jar SE -threads 16 %s ./%s/%s_TRIM.fq.gz ILLUMINACLIP:/scratch-lustre/pscha005/trimmomatic/adapters/TruSeq3-SE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
date
'''%(sample, sampleDict[sample][0], sample, sample))
			trimScript.close()
		elif len(sampleDict[sample]) == 2 and hpc == 0:
			try:
				subprocess.call(['java', '-jar', 'trimmomatic-0.33.jar', 'PE', '-phred33', '%s' %(sampleDict[sample][0]), '%s' %(sampleDict[sample][1]), './%s/%s_TRIM_forward_paired.fq.gz'%(sample,sample), './%s/%s_TRIM_forward_unpaired.fq.gz'%(sample,sample), './%s/%s_TRIM_reverse_paired.fq.gz'%(sample,sample), './%s/%s_TRIM_reverse_unpaired.fq.gz'%(sample,sample), 'ILLUMINACLIP:./adapters/TruSeq3-PE.fa:2:30:10', 'LEADING:3', 'TRAILING:3', 'SLIDINGWINDOW:4:15', 'MINLEN:36'])
			except:
				sys.exit("ERROR: Trimmomatic.jar not found in current directory")
			trimOutputDict[sample] = ['%s_TRIM_forward_paired.fq.gz'%(sample), '%s_TRIM_reverse_paired.fq.gz'%(sample)]
		elif len(sampleDict[sample]) == 2 and hpc == 1:
			trimOutputDict[sample] = ['%s_TRIM_forward_paired.fq.gz'%(sample), '%s_TRIM_reverse_paired.fq.gz'%(sample)]
			trimScript = open("trimmomatic_%s.sh" %(sample), "w")
			trimScript.write(
'''#!/bin/bash -l\n
#SBATCH --job-name=trim    # Job name
#SBATCH --ntasks=1                    # Run on a single CPU
#SBATCH --cpus-per-task=16
#SBATCH --output=trimmomatic_%s.out   # Standard output and error log

pwd; hostname; date

enable_lmod
module load java/11.0
java -jar /scratch-lustre/pscha005/trimmomatic/trimmomatic-0.33.jar PE -threads 16 %s %s ./%s/%s_TRIM_forward_paired.fq.gz ./%s/%s_TRIM_forward_unpaired.fq.gz ./%s/%s_TRIM_reverse_paired.fq.gz ./%s/%s_TRIM_reverse_unpaired.fq.gz ILLUMINACLIP:/scratch-lustre/pscha005/trimmomatic/adapters/TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
date
'''%(sample, sampleDict[sample][0], sampleDict[sample][1], sample, sample, sample, sample, sample, sample, sample, sample))
			trimScript.close()
			novoplastyScript = open("novoplasty_%s.sh" %(sample), "w")
			novoplastyScript.write(
'''#!/bin/bash -l\n
#SBATCH --job-name=novoplast    # Job name
#SBATCH --ntasks=1                    # Run on a single CPU
#SBATCH --cpus-per-task=1
#SBATCH --output=novoplasty_%s.out   # Standard output and error log

pwd; hostname; date

enable_lmod
module load perl/5.8
perl /scratch-lustre/pscha005/scripts/NOVOPlasty/NOVOPlasty2.7.2.pl -c novoplasty_config_%s_plastome.txt
date
'''%(sample, sample))
			novoplastyScript.close()
			novoplastConfig = open("novoplasty_config_%s_plastome.txt" %(sample), "w")
			novoplastConfig.write(
'''Project:
-----------------------
Project name          = %s_plastome
Type                  = chloro
Genome Range          = 140000-150000
K-mer                 = 39
Max memory            = 
Extended log          = 0
Save assembled reads  = no
Seed Input            = /scratch-lustre/pscha005/scripts/NOVOPlasty/Isoetes_echinospora_rbcL.fasta
Reference sequence    = /scratch-lustre/pscha005/Illumina/Isoetes_flaccida_GU191333.fasta
Variance detection    = no
Heteroplasmy          =
HP exclude list       =
Chloroplast sequence  =

Dataset 1:
-----------------------
Read Length           = 151
Insert size           = 300
Platform              = illumina
Single/Paired         = PE
Combined reads        =
Forward reads         = %s
Reverse reads         = %s

Optional:
-----------------------
Insert size auto      = yes
Insert Range          = 1.8
Insert Range strict   = 1.3
Use Quality Scores    = no


Project:
-----------------------
Project name         = Choose a name for your project, it will be used for the output files.
Type                 = (chloro/mito/mito_plant) "chloro" for chloroplast assembly, "mito" for mitochondrial assembly and 
                       "mito_plant" for mitochondrial assembly in plants.
Genome Range         = (minimum genome size-maximum genome size) The expected genome size range of the genome.
                       Default value for mito: 12000-20000 / Default value for chloro: 120000-200000
                       If the expected size is know, you can lower the range, this can be useful when there is a repetitive
                       region, what could lead to a premature circularization of the genome.
K-mer                = (integer) This is the length of the overlap between matching reads (Default: 39). 
                       If reads are shorter then 90 bp or you have low coverage data, this value should be decreased down to 23. 
                       For reads longer then 101 bp, this value can be increased, but this is not necessary.
Max memory           = You can choose a max memory usage, suitable to automatically subsample the data or when you have limited                      
                       memory capacity. If you have sufficient memory, leave it blank, else write your available memory in GB
                       (if you have for example a 8 GB RAM laptop, put down 7 or 7.5 (don't add the unit in the config file))
Extended log         = Prints out a very extensive log, could be useful to send me when there is a problem  (0/1).
Save assembled reads = All the reads used for the assembly will be stored in seperate files (yes/no)
Seed Input           = The path to the file that contains the seed sequence.
Reference (optional) = If a reference is available, you can give here the path to the fasta file.
                       The assembly will still be de novo, but references of the same genus can be used as a guide to resolve 
                       duplicated regions in the plant mitochondria or the inverted repeat in the chloroplast. 
                       References from different genus haven't beeen tested yet.
Variance detection   = If you select yes, you should also have a reference sequence (previous line). It will create a vcf file                
                       with all the variances compared to the give reference (yes/no)
Heteroplasmy         = If yo uwant to detect heteroplasmy,first assemble the genome without this option. Then give the resulting                         
                       sequence as a reference and as a seed input. And give the minimum minor allele frequency for this option 
                       (0.01 will detect heteroplasmy of >1%%)
HP exclude list      = Option not yet available  
Chloroplast sequence = The path to the file that contains the chloroplast sequence (Only for mito_plant mode).
                       You have to assemble the chloroplast before you assemble the mitochondria of plants!

Dataset 1:
-----------------------
Read Length          = The read length of your reads.
Insert size          = Total insert size of your paired end reads, it doesn't have to be accurate but should be close enough.
Platform             = illumina is for now the only option
Single/Paired        = For the moment only paired end reads are supported.
Combined reads       = The path to the file that contains the combined reads (forward and reverse in 1 file)
Forward reads        = The path to the file that contains the forward reads (not necessary when there is a merged file)
Reverse reads        = The path to the file that contains the reverse reads (not necessary when there is a merged file)

Optional:
-----------------------
Insert size auto     = (yes/no) This will finetune your insert size automatically (Default: yes)
Insert Range         = This variation on the insert size, could lower it when the coverage is very high or raise it when the
                       coverage is too low (Default: 1.6). 
Insert Range strict  = Strict variation to resolve repetitive regions (Default: 1.2).                                
Use Quality Scores   = It will take in account the quality scores, only use this when reads have low quality, like with the    
                       300 bp reads of Illumina (yes/no)
''' %(sample, sampleDict[sample][0],sampleDict[sample][1]))
			novoplastConfig.close()

		elif len(sampleDict[sample]) > 2:
			print "!" * 10
			print "ERROR: Too many files specified for %s" %(sample)
			print "!" * 10
			sys.exit()

'''Bowtie2'''
def bowtie():
	if 'referenceName' in globals():
		print "*"*10
		print "Running bowtie..."
		print "*"*10

		reference = "reference"
		subprocess.call(['mkdir', reference])
		referenceOutput = "./reference/%s" %(referenceName)
		try:
			subprocess.call(['bowtie2-build', '-q', referenceFile, referenceOutput])
		except:
			sys.exit("ERROR: bowtie2 not found in PATH")
		global bowtieOutputDict
		bowtieOutputDict = {}
		for sample in sorted(trimOutputDict.keys()):
			if len(trimOutputDict[sample]) == 1 and hpc == 0:
				subprocess.call(['bowtie2','--quiet', '--very-sensitive-local', '-x', referenceOutput, '-U', './%s/%s' %(sample,trimOutputDict[sample][0]), '--al-gz', './%s/%s_align.fq.gz' %(sample,sample), '--un-gz', './%s/%s_unalign.fq.gz' %(sample,sample)])
				bowtieOutputDict[sample] = ["%s_align.fq.gz" %(sample), "%s_unalign.fq.gz" %(sample)]
			elif len(trimOutputDict[sample]) == 1 and hpc == 1:
				bowtieOutputDict[sample] = ["%s_align.fq.gz" %(sample), "%s_unalign.fq.gz" %(sample)]
				bowtieScript = open("bowtie_%s.sh" %(sample), "w")
				bowtieScript.write(
'''#!/bin/bash\n
#SBATCH --job-name=bwtie    # Job name
#SBATCH --ntasks=1                    # Run on a single CPU
#SBATCH --output=bowtie_%s.out

pwd; hostname; date

module load bowtie2/2.2.4
bowtie2 --quiet --very-sensitive-local -x %s -U ./%s/%s --al-gz ./%s/%s_align.fq.gz --un-gz ./%s/%s_unalign.fq.gz

date
'''%(sample, referenceOutput, sample, trimOutputDict[sample][0], sample, sample, sample, sample))
				bowtieScript.close()
			elif len(trimOutputDict[sample]) == 2 and hpc == 0:
				subprocess.call(['bowtie2','--quiet', '--very-sensitive-local', '-x', referenceOutput, '-1', './%s/%s' %(sample,trimOutputDict[sample][0]), '-2', './%s/%s' %(sample,trimOutputDict[sample][1]), '--al-conc-gz', './%s/%s_R%%_align.fq.gz' %(sample,sample), '--un-conc-gz', './%s/%s_R%%_unalign.fq.gz' %(sample,sample)])
				bowtieOutputDict[sample] = ["%s_R1_align.fq.gz" %(sample), "%s_R2_align.fq.gz" %(sample), "%s_R1_unalign.fq.gz" %(sample), "%s_R2_unalign.fq.gz" %(sample), "%s_R1_align_rDNA.fq.gz" %(sample), "%s_R2_align_rDNA.fq.gz" %(sample)]
			elif len(trimOutputDict[sample]) == 2 and hpc == 1:
				bowtieOutputDict[sample] = ["%s_R1_align.fq.gz" %(sample), "%s_R2_align.fq.gz" %(sample), "%s_R1_unalign.fq.gz" %(sample), "%s_R2_unalign.fq.gz" %(sample), "%s_R1_align_rDNA.fq.gz" %(sample), "%s_R2_align_rDNA.fq.gz" %(sample)]
				bowtieScript = open("bowtie_%s.sh" %(sample), "w")
				bowtieScriptRDNA = open("bowtie_%s_rDNA.sh" %sample, "w")
				bowtieScript.write(
'''#!/bin/bash\n
#SBATCH --job-name=bwtie    # Job name
#SBATCH --ntasks=1                    # Run on a single CPU
#SBATCH --output=bowtie_%s.out

pwd; hostname; date

module load bowtie2/2.2.4
bowtie2 --quiet --very-sensitive-local -x %s -1 ./%s/%s -2 ./%s/%s --al-conc-gz ./%s/%s_R%%_align.fq.gz --un-conc-gz ./%s/%s_R%%_unalign.fq.gz

date
'''%(sample, referenceOutput, sample, trimOutputDict[sample][0], sample, trimOutputDict[sample][1], sample, sample, sample, sample))
				bowtieScript.close()
				bowtieScriptRDNA.write(
'''#!/bin/bash\n
#SBATCH --job-name=bwtie    # Job name
#SBATCH --ntasks=1                    # Run on a single CPU
#SBATCH --output=bowtie_%s_rDNA.out

pwd; hostname; date

module load bowtie2/2.2.4
bowtie2 --quiet --very-sensitive-local -x ./reference/Isoetes_rDNA -1 ./%s/%s -2 ./%s/%s --al-conc-gz ./%s/%s_R%%_align_rDNA.fq.gz

date
'''%(sample, sample, trimOutputDict[sample][0], sample, trimOutputDict[sample][1], sample, sample))
				bowtieScriptRDNA.close()
	else:
		print "*"*10
		print "WARNING: No reference file found, skipping bowtie filtering..."
		print "*"*10


'''SPAdes'''
def spades():
	for sample in sorted(bowtieOutputDict.keys()):
		if hpc == 0 and len(bowtieOutputDict[sample]) == 2:
			subprocess.call(['spades.py', '-k', '21,33,55,77', '--careful', '-s', './%s/%s' %(sample, bowtieOutputDict[sample][0]),'-o', './%s/spades_%s' %(sample, sample)])
		elif hpc == 1 and len(bowtieOutputDict[sample]) == 2:
			spadesScript = open("spades_%s_plastome.sh" %(sample), "w")
			nonplastScript = open("spades_%s_nonPlastome.sh" %(sample), "w")
			spadesScript.write(
'''#!/bin/bash -l
#SBATCH --job-name=spds    # Job name
#SBATCH --mail-type=ALL               # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=<pscha005@odu.edu>   # Where to send mail   
#SBATCH --ntasks=1                    # Run on a single CPU
#SBATCH --cpus-per-task=16
#SBATCH --output=spades_%s_plastome.out   # Standard output and error log

pwd; hostname; date

enable_lmod
module load spades/3.13
spades.py -t 16 -k 21,33,55,77 --careful -s ./%s/%s -o ./%s/spades_%s_plastome

date
'''%(sample, sample, bowtieOutputDict[sample][0], sample, sample))
			spadesScript.close()
			nonplastScript.write(
'''#!/bin/bash -l
#SBATCH --job-name=spds    # Job name
#SBATCH --mail-type=ALL               # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=<pscha005@odu.edu>   # Where to send mail   
#SBATCH --ntasks=1                    # Run on a single CPU
#SBATCH --cpus-per-task=16
#SBATCH --output=spades_%s_nonPlastome.out   # Standard output and error log

pwd; hostname; date

enable_lmod
module load spades/3.13
spades.py -t 16 -k 21,33,55,77 -s ./%s/%s -o ./%s/spades_%s_nonPlastome

date
'''%(sample, sample, bowtieOutputDict[sample][1], sample, sample))
			nonplastScript.close()
			allReadsScript = open("spades_%s_allReads.sh" %(sample), "w")
			allReadsScript.write(
'''#!/bin/bash -l
#SBATCH --job-name=spds    # Job name
#SBATCH --mail-type=ALL               # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=<pscha005@odu.edu>   # Where to send mail   
#SBATCH --ntasks=1                    # Run on a single CPU
#SBATCH --cpus-per-task=16
#SBATCH --output=spades_%s_allReads.out   # Standard output and error log

pwd; hostname; date

enable_lmod
module load spades/3.13
spades.py -t 16 -k 21,33,55,77 -s ./%s/%s -o ./%s/spades_%s_allReads

date
'''%(sample, sample, trimOutputDict[sample][0], sample, sample))
			allReadsScript.close()
		elif hpc == 0 and len(bowtieOutputDict[sample]) == 6:
			subprocess.call(['spades.py', '-k', '21,33,55,77', '--careful', '-1', './%s/%s' %(sample, bowtieOutputDict[sample][0]), '-2', './%s/%s' %(sample, bowtieOutputDict[sample][1]),'-o', './%s/spades_%s' %(sample, sample)])
		elif hpc == 1 and len(bowtieOutputDict[sample]) == 6:
			spadesScript = open("spades_%s_plastome.sh" %(sample), "w")
			spadesScript.write(
'''#!/bin/bash -l
#SBATCH --job-name=spds    # Job name
#SBATCH --mail-type=ALL               # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=<pscha005@odu.edu>   # Where to send mail   
#SBATCH --ntasks=1                    # Run on a single CPU
#SBATCH --cpus-per-task=16
#SBATCH --output=spades_%s_plastome.out   # Standard output and error log

pwd; hostname; date

enable_lmod
module load spades/3.13
spades.py -t 16 -k 21,33,55,77 --careful -1 ./%s/%s -2 ./%s/%s -o ./%s/spades_%s_plastome

date
'''%(sample, sample, bowtieOutputDict[sample][0], sample, bowtieOutputDict[sample][1], sample, sample))
			spadesScript.close()
			nonplastScript = open("spades_%s_nonPlastome.sh" %(sample), "w")
			nonplastScript.write(
'''#!/bin/bash -l
#SBATCH --job-name=spds    # Job name
#SBATCH --mail-type=ALL               # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=<pscha005@odu.edu>   # Where to send mail   
#SBATCH --ntasks=1                    # Run on a single CPU
#SBATCH --cpus-per-task=16
#SBATCH --output=spades_%s_nonPlastome.out   # Standard output and error log

pwd; hostname; date

enable_lmod
module load spades/3.13
spades.py -t 16 -k 21,33,55,77 -1 ./%s/%s -2 ./%s/%s -o ./%s/spades_%s_nonPlastome

date
'''%(sample, sample, bowtieOutputDict[sample][2], sample, bowtieOutputDict[sample][3], sample, sample))
			nonplastScript.close()
			allReadsScript = open("spades_%s_allReads.sh" %(sample), "w")
			allReadsScript.write(
'''#!/bin/bash -l
#SBATCH --job-name=spds    # Job name
#SBATCH --mail-type=ALL               # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=<pscha005@odu.edu>   # Where to send mail   
#SBATCH --ntasks=1                    # Run on a single CPU
#SBATCH --cpus-per-task=16
#SBATCH --output=spades_%s_allReads.out   # Standard output and error log

pwd; hostname; date

enable_lmod
module load spades/3.13
spades.py -t 16 -k 21,33,55,77 -1 ./%s/%s -2 ./%s/%s -o ./%s/spades_%s_allReads

date
'''%(sample, sample, trimOutputDict[sample][0], sample, trimOutputDict[sample][1], sample, sample))
			allReadsScript.close()
			rDNAScript = open("spades_%s_rDNA.sh" %(sample), "w")
			rDNAScript.write(
'''#!/bin/bash -l
#SBATCH --job-name=spds    # Job name
#SBATCH --mail-type=ALL               # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=<pscha005@odu.edu>   # Where to send mail   
#SBATCH --ntasks=1                    # Run on a single CPU
#SBATCH --cpus-per-task=16
#SBATCH --output=spades_%s_rDNA.out   # Standard output and error log

pwd; hostname; date

enable_lmod
module load spades/3.13
spades.py -t 16 -k 21,33,55,77 -1 ./%s/%s -2 ./%s/%s -o ./%s/spades_%s_rDNA

date
'''%(sample, sample, bowtieOutputDict[sample][4], sample, bowtieOutputDict[sample][5], sample, sample))
			rDNAScript.close()

trimmomatic()
bowtie()
spades()
