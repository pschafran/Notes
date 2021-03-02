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
outputfile = open('%s_BLASTEXTRACT.fasta'%(blastfiles[0]), 'w')
hitlist=[]
query="Empty"
for file in blastfiles:
	openfile = open(file, 'r')
	for line in openfile:
		#print query
		splitline = line.split(',')
		query=splitline[1]
		hitlist.append(query)
		#print line
	#print hitdict
	openfile.close()
#hitdict2={}
#for k in hitdict.keys():
#	if k not in hitdict2.keys():
#		hitdict2[k]=[]
#print hitdict2
#for key in hitdict.keys():
	#print key
#	del hitdict[key][0:2]
#	del hitdict[key][-3:]
#	for item in hitdict[key]:
#		itemsplit=item.split(' ')
#		wantthis=itemsplit[2:3]
		#print wantthis
#		hitdict2[key].append(wantthis)
#print hitdict2
resultsdict={}
for k in hitlist:
	if k not in resultsdict.keys():
		resultsdict[k]=[]
print "Extracting hits..."
#hitnum=1
gc.enable()
for file in dbfiles:
	print file
	for key in set(hitlist):
		counter=0
		with open(file, 'r') as searchfile:
			for line in searchfile:
				if ">" in line:
					counter=0
					if key in line:
						outputfile.write(line)
						counter=1
				elif counter>0:
					outputfile.write(line)
					counter += 1
		
			#print line
#			if count_on==1:
#				#print "if count_on==1"
#				count_on=0
#				for key2 in resultsdict.keys():
#					for item3 in resultsdict[key2]:
#						#print "Item3: %s" %(item3)
#						#print "Item2: %s" %(item2)
#						
#						if item3==item2:
#							#print "Item3 Match!"
#							resultsdict[key2].append(line)
#							print "Hit number: %d   Key: %s" %(hitnum, key2)
#							hitnum += 1
#							break
#			if '>' in line:
#				#print "if > in line"
#				for item2 in hitdict.keys():
#					#print key
#					for item in hitdict[key]:
#						print item
#						for item2 in item:
#							print item2
#							print "Item2: %s" %(item2)
#					if item2 in line:
#						count_on=1
#						#print "Item2 Match!", key
#						resultsdict[item2].append(item2)
#						print "Loop 1"
searchfile.close()
outputfile.close()
		