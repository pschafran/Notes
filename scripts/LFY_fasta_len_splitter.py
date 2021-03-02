#! /usr/bin/python
'''Sorts .fasta files into sequences <900bp, >900 and <1400bp, and >1400bp
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
	
	fastaDict = {}
	for line in infile:
		if ">" in line:
			key = line
			fastaDict[key] = []
			joinLine = "".join(fastaDict[key])
			fastaDict[key] = joinLine
		if ">" not in line:
				stripLine = line.strip("\n")
				try:
					fastaDict[key].append(stripLine)
				except:
					fastaDict[key] = stripLine
	joinLine = "".join(fastaDict[key])
	fastaDict[key] = joinLine
	infile.close()
	for key in fastaDict.keys():
			if len(fastaDict[key]) < 1000:
				toosmall_outfile.write("%s\n" %(key))
				toosmall_outfile.write("%s\n" %(fastaDict[key]))
				toosmallseqs+=1
			elif len(fastaDict[key]) >=1000 and len(fastaDict[key]) <1200:
				rightsize_outfile.write("%s\n" %(key))
				rightsize_outfile.write("%s\n" %(fastaDict[key]))
				rightsizeseqs+=1
			elif len(fastaDict[key]) >=1200:
				toobig_outfile.write("%s\n" %(key))
				toobig_outfile.write("%s\n" %(fastaDict[key]))
				toobigseqs+=1
	totalseqs=toosmallseqs+rightsizeseqs+toobigseqs
	percenttoosmall = (float(toosmallseqs)/float(totalseqs))*100
	percentrightsize = (float(rightsizeseqs)/float(totalseqs))*100
	percenttoobig = (float(toobigseqs)/float(totalseqs))*100
	print "============"
	print "Output Stats: %s" %(file)
	print "Total Seqs Processed: %d" %(totalseqs)
	print "Seqs <900bp: %d (%.1f%%)" %(toosmallseqs,percenttoosmall)
	print "Seqs >900bp and <1400bp: %d (%.1f%%)" %(rightsizeseqs,percentrightsize)
	print "Seqs >1400bp: %d (%.1f%%)" %(toobigseqs,percenttoobig)
	toosmall_outfile.close()
	rightsize_outfile.close()
	toobig_outfile.close()
