#! /usr/bin/python
'''Sorts .fasta files into sequences <500bp, >500 and <1000bp, and >1000bp
Usage: ./thisscript.py files.fasta
'''
import sys

filelist=sys.argv[1:]

for file in filelist:
	splitfile=file.split(".")
	infile = open(file, 'r')
	rightsize_outfile = open("%s_rightsize.fasta" %(splitfile[0]), 'w')
	toosmall_outfile = open("%s_toosmall.fasta" %(splitfile[0]), 'w')
	toobig_outfile = open("%s_toobig.fasta" %(splitfile[0]), 'w')
	toosmallseqs=0
	rightsizeseqs=0
	toobigseqs=0
	for line in infile:
		if ">" in line:
			seqname=line
		else:
			seq=line
			if len(seq) < 500:
				toosmall_outfile.write("%s\n" %(seqname))
				toosmall_outfile.write("%s\n" %(seq))
				toosmallseqs+=1
			elif len(seq) >=500 and len(seq) <1000:
				rightsize_outfile.write("%s\n" %(seqname))
				rightsize_outfile.write("%s\n" %(seq))
				rightsizeseqs+=1
			elif len(seq) >=1000:
				toobig_outfile.write("%s\n" %(seqname))
				toobig_outfile.write("%s\n" %(seq))
				toobigseqs+=1
	totalseqs=toosmallseqs+rightsizeseqs+toobigseqs
	print "Output Stats"
	print "============"
	print "Total Seqs Processed: %d" %(totalseqs)
	print "Seqs <500bp: %d" %(toosmallseqs)
	print "Seqs >500bp and <1000bp: %d" %(rightsizeseqs)
	print "Seqs >1000bp: %d" %(toobigseqs)
