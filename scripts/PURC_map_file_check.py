#! /usr/bin/python
##2017 Nov 16
##Peter W. Schafran

##Checks for duplicate barcode combinations and sample names in tab-delimited locus-barcode-taxon-map file for PURC (Pipeline for Unraveling Reticulate Complexes [Rothfels et al. 2017])
##Usage: PURC_map_file_check.py mapFile.txt

import sys


file = sys.argv[1]
openfile = open(file, 'r')

barcodedict ={}
linenumber = 1
ncol = 0
for line in openfile:
	splitline = line.strip('\n').split('\t')
	barcodedict[linenumber] = splitline
	linenumber += 1
	ncol = len(splitline)
barcodefail = 1
samplefail = 1
if ncol == 3:
	for key in barcodedict.keys():
		BCF = barcodedict[key][0]
		BCR = barcodedict[key][1]
		sample = barcodedict[key][2]
		for key2 in barcodedict.keys():
			if key2 == key:
				pass
			else:
				if BCF == barcodedict[key2][0] and BCR == barcodedict[key2][1]:
					print "Duplicate barcode combination found: %s %s at lines %s and %s" %(BCF,BCR,key,key2)
					passfail = 0
				if sample == barcodedict[key2][2]:
					print "Duplicate sample name found: %s at lines %s and %s" %(sample,key,key2)
					samplefail = 0
elif ncol == 2:
	for key in barcodedict.keys():
		BCF = barcodedict[key][0]
		sample = barcodedict[key][1]
		for key2 in barcodedict.keys():
			if key2 == key:
				pass
			else:
				if BCF == barcodedict[key2][0]:
					print "Duplicate barcode combination found: %s at lines %s and %s" %(BCF,key,key2)
					passfail = 0
				if sample == barcodedict[key2][1]:
					print "Duplicate sample name found: %s at lines %s and %s" %(sample,key,key2)
					samplefail = 0
if barcodefail == 1:
	print "No duplicate barcodes found!"
if samplefail == 1:
	print "No duplicate samples found!"
