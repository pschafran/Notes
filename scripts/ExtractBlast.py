#! /usr/bin/python
'''Usage: blastresultfile1.out...blastresultfile2.out...databasesource.fasta'''
import sys, re, gc
filelist=sys.argv[1:]
blastfiles=[]
dbfiles=[]
for file in filelist:
	splitfile=file.split('.')
	#print splitfile[1:]
	fileextension1='out'
	fileextension2='fasta'
	if fileextension1 in splitfile[1:]:
		blastfiles.append(file)
	if fileextension2 in splitfile[1:]:
		dbfiles.append(file)
#print filelist
#print blastfiles
#print dbfiles
hitdict={}
query="Empty"
for file in blastfiles:
	count_on=0
	openfile = open(file, 'r')
	for line in openfile:
		#print query
		#print line
		if 'Query=' in line:
			#print line
			splitline=line.strip('\n').split('=')
			query=splitline[1:]
			#print query
			hitdict[query[0]]=[]
		if 'significant' in line:
			count_on=1
		if count_on==1:
			hitdict[query[0]].append(line.strip('\n'))
			#print line
		if '>' in line:
			count_on=0
			
	openfile.close()
hitdict2={}
for k in hitdict.keys():
	if k not in hitdict2.keys():
		hitdict2[k]=[]
#print hitdict2
for key in hitdict.keys():
	#print key
	del hitdict[key][0:2]
	del hitdict[key][-3:]
	for item in hitdict[key]:
		itemsplit=item.split(' ')
		wantthis=itemsplit[2:3]
		#print wantthis
		hitdict2[key].append(wantthis)
#print hitdict2
resultsdict={}
for k in hitdict.keys():
	if k not in resultsdict.keys():
		resultsdict[k]=[]
print "Extracting hits..."
hitnum=0
gc.enable()
for file in dbfiles:
	print file
	with open(file, 'r') as searchfile:
		count_on=0
		for line in searchfile:
			#print line
			if count_on==1:
				#print "if count_on==1"
				count_on=0
				for key2 in resultsdict.keys():
					for item3 in resultsdict[key2]:
						#print "Item3: %s" %(item3)
						#print "Item2: %s" %(item2)
						
						if item3==item2:
							#sprint "Item3 Match!"
							resultsdict[key2].append(line)
							print "Hit number: %d   Key: %s" %(hitnum, key2)
							hitnum += 1
							break
			if '>' in line:
				#print "if > in line"
				for key in hitdict2.keys():
					for item in hitdict2[key]:
						for item2 in item:
							#print "Item2: %s" %(item2)
							if item2 in line:
								count_on=1
								#print "Item2 Match!", key
								resultsdict[key].append(item2)
								
								break
						else:
							continue
						break
outputfile = open('%s_BLASTEXTRACT.fasta'%(blastfiles[0]), 'w')
print '''
=================================================
Writing to output file: %s_BLASTEXTRACT.fasta
=================================================''' %(blastfiles[0])
counter=0
for key in resultsdict.keys():
	for item in resultsdict[key]:
		readstrip = resultsdict[key][counter:counter+1]
		readstrip1 = str(readstrip).strip("[]")
		readstrip2 = readstrip1.strip("\'")
		print readstrip2
		#readstrip1 = readstrip[0].translate(None, "[']")
		seqstrip = resultsdict[key][counter+1:counter+2]
		seqstrip1 = str(seqstrip).strip("[]")
		seqstrip2 = seqstrip1.strip("\'")
		#seqstrip1 = seqstrip[0].translate(None, "[']")
		outputfile.write('>%s\n' %(readstrip2))
		outputfile.write('%s\n' %(seqstrip2))
		counter += 2
searchfile.close()
outputfile.close()
		