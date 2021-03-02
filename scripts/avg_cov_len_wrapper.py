#! /usr/bin/python
#Wrapper for avg_cov_len_fasta_DJB.py
#
#Usage: avg_cov_len_fasta_DJB.py files.fasta

import sys,subprocess,glob,os


for file in sys.argv[1:]:
	split = file.split('/')
	filename = split[1]
	outfile=open('OUTPUT_%s_avg_cov_len_fasta_DJB.txt' %(filename),'w')
	subprocess.call(['/Users/Peter/python/scripts/avg_cov_len_fasta_DJB.py','%s' %(file)], stdout=outfile)
masteroutputfile=open("FINAL_OUTPUT_avg_cov_len_fasta_DJB.txt", 'w')
for file2 in glob.glob("OUTPUT*txt"):
	openfile=open(file2, 'r')
	print file2
	for line in openfile:
		print line
