#! /usr/bin/python

file = "BarcodedLFYPrimers_30July2015.csv"

primerlist=[]
openfile = open(file, 'r')
count1=0
count2=0
for line in openfile:
	splitline = line.strip('\n')
	splitline =splitline.split(',')
	if splitline[2] in primerlist:
		print "Duplicate primer! %s" %(line)
		primerlist.append(splitline[2])
		count1-=1
		continue
	else:
		primerlist.append(splitline[2])
		count1+=1
		count2+=1
		continue
	continue

if count1 == count2:
	print "All barcodes unique!"