#! /usr/bin/python
import sys
#if sys.argv[1] == 'usage':
print '''==========================================================
Usage: ConcatenateNexus.py OutfileName.nex InputFileName1.nex InputFileName2.nex ...
-Input Nexus files must have same taxa names
-Only combines 2 files at a time
=========================================================='''
	


filelist= sys.argv[2:]
#print sys.argv[1:2]
outputfile = sys.argv[1:2]
#print outputfile
outfile = open(outputfile[0], 'w')
filecounter = 1
infiledict = {}
for file in filelist:
	print "Processing: %s" % file
	infile = open(file, 'r')
	linecounter=1
	taxacounter=1
	taxalist=[]
	seqlist={}
	ntax=[]
	nchar=[]
	infiledict[file] = {}
	for line in infile:
		if linecounter==3:
			splitline = line.strip(';\n').split(' ')
			#print splitline
			splittax = splitline[1].split('=')
			#print splittax
			ntax = int(splittax[1])
			infiledict[file]['ntax'] = ntax
			print "Number of taxa: %s" % ntax
		if linecounter>=5 and taxacounter<=ntax:
			splitlabel = line.strip('\t').split('[')
			taxalist.append(splitlabel[0])
			
			taxacounter+=1
		if linecounter >=5 and linecounter == (ntax + 9):
			#print line
			splitline2 = line.strip(';\n').split(' ')
			splitchar = splitline2[1].split('=')
			nchar = int(splitchar[1])
			infiledict[file]['nchar'] = nchar
			print "Number of characters: %s" % nchar
		if linecounter >=5 and linecounter >= (ntax + 12):
			splitseq = line.strip('\n').split('\t')
			if len(splitseq) < 3:
				break
			seqlist[splitseq[1]] = splitseq[2]
			
		linecounter+=1
		infiledict[file]['seqlist'] = seqlist
		infiledict[file]['taxalist'] = taxalist
		#print infiledict
	#print taxalist
	#print seqlist
	#print infiledict
	filecounter += 1
	
outfile.write("#NEXUS\n")
outfile.write("begin data;\n")
outfile.write("\tdimensions ntax=%s nchar=%s;\n" % (infiledict[filelist[0]]['ntax'], infiledict[filelist[0]]['nchar']+infiledict[filelist[1]]['nchar']))
outfile.write("\tformat datatype=dna missing=? gap=-;\n")
outfile.write("matrix\n")
for key in infiledict[filelist[0]]['seqlist']:
	if key in infiledict[filelist[1]]['seqlist'].keys():
		outfile.write("%s %s%s\n" %(key, infiledict[filelist[0]]['seqlist'][key[0:]], infiledict[filelist[1]]['seqlist'][key[0:]]))
	else:
		print '''==========================================================
ERROR: NAMES IN INPUT FILES DO NOT MATCH!
=========================================================='''
		break
		#print key
	#print infiledict[filelist[0]]['seqlist'][key[0:]]
	#print infiledict[filelist[1]]['seqlist'][key[0:]]
	#print ' '
#print infiledict
outfile.write(";\n")
outfile.write("end;\n")
outfile.write("\n")
outfile.write("begin mrbayes;\n")
outfile.write("charset %s = 1-%s;\n" %(filelist[0], infiledict[filelist[0]]['nchar']))
outfile.write("charset %s = %s-%s;\n" %(filelist[1], infiledict[filelist[0]]['nchar']+1, infiledict[filelist[0]]['nchar']+infiledict[filelist[1]]['nchar']))
outfile.write("partition combined = 2: %s,%s;\n" %(filelist[0], filelist[1]))
outfile.write("end;\n")