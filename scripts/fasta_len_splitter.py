#! /usr/bin/python
'''Sorts .fasta files into sequences ....
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
	totalseqs=0
	key = 0
	bin0 = 0
	bin200 = 0
	bin400 = 0
	bin600 = 0
	bin800 = 0
	bin1000 = 0
	bin1200 = 0
	bin1400 = 0
	bin1600 = 0
	bin1800 = 0
	bin2000 = 0
	fastaDict = {}
	for line in infile:
		if ">" in line:
			while key != 0:
				joinLine = "".join(fastaDict[key])
				fastaDict[key] = joinLine
				key = 0
			key = line.strip("\n")
			totalseqs += 1
			fastaDict[key] = []
			if totalseqs == 1000:
				print "1k seqs processed..."
			if totalseqs == 10000:
				print "10k seqs processed..."
			if totalseqs == 50000:
				print "50k seqs processed..."
			if totalseqs == 100000:
				print "100k seqs processed..."
			if totalseqs == 250000:
				print "250k seqs processed..."
			if totalseqs == 500000:
				print "500k seqs processed..."
			if totalseqs == 750000:
				print "750k seqs processed..."
			if totalseqs == 1000000:
				print "1M seqs processed..."
		if ">" not in line:
			stripLine = line.strip("\n")
			fastaDict[key].append(stripLine)
	joinLine = "".join(fastaDict[key])
	fastaDict[key] = joinLine
	for key in fastaDict.keys():
		seq = fastaDict[key]
		if len(seq) < 700:
			toosmall_outfile.write("%s\n" %(key))
			toosmall_outfile.write("%s\n" %(seq))
			toosmallseqs+=1
		elif len(seq) >=700 and len(seq) <1500:
			rightsize_outfile.write("%s\n" %(key))
			rightsize_outfile.write("%s\n" %(seq))
			rightsizeseqs+=1
		elif len(seq) >=1500:
			toobig_outfile.write("%s\n" %(key))
			toobig_outfile.write("%s\n" %(seq))
			toobigseqs+=1
		if len(seq) < 200:
			bin0 += 1
		if len(seq) >= 200 and len(seq) < 400:
			bin200 += 1
		if len(seq) >= 400 and len(seq) < 600:
			bin400 += 1
		if len(seq) >= 600 and len(seq) < 800:
			bin600 += 1
		if len(seq) >= 800 and len(seq) < 1000:
			bin800 += 1
		if len(seq) >= 1000 and len(seq) < 1200:
			bin1000 += 1
		if len(seq) >= 1200 and len(seq) < 1400:
			bin1200 += 1
		if len(seq) >= 1400 and len(seq) < 1600:
			bin1400 += 1
		if len(seq) >= 1600 and len(seq) < 1800:
			bin1600 += 1
		if len(seq) >= 1800 and len(seq) < 2000:
			bin1800 += 1
		if len(seq) >= 2000:
			bin2000 += 1
	totalseqs=toosmallseqs+rightsizeseqs+toobigseqs
	percenttoosmall = (float(toosmallseqs)/float(totalseqs))*100
	percentrightsize = (float(rightsizeseqs)/float(totalseqs))*100
	percenttoobig = (float(toobigseqs)/float(totalseqs))*100
	print "============"
	print "Output Stats"
	print "============"
	print "Total Seqs Processed: %d" %(totalseqs)
	print "Seqs <800bp: %d (%.1f%%)" %(toosmallseqs,percenttoosmall)
	print "Seqs >800bp and <1400bp: %d (%.1f%%)" %(rightsizeseqs,percentrightsize)
	print "Seqs >1400bp: %d (%.1f%%)" %(toobigseqs,percenttoobig)
	print "============"
	print "bin0\t%s" %(bin0)
	print "bin200\t%s" %(bin200)
	print "bin400\t%s" %(bin400)
	print "bin600\t%s" %(bin600)
	print "bin800\t%s" %(bin800)
	print "bin1000\t%s" %(bin1000)
	print "bin1200\t%s" %(bin1200)
	print "bin1400\t%s" %(bin1400)
	print "bin1600\t%s" %(bin1600)
	print "bin1800\t%s" %(bin1800)
	print "bin2000\t%s" %(bin2000)
	

