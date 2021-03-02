#! /usr/bin/python

## 2019 September 2 
## Last updated 2019 September 23
## Peter W. Schafran
## This script takes barcode sequences and searches for and counts them in a sequence file. Reports number of each barcode found at forward, reverse, combined (sum of forward + reverse)


## Files: (1) PURC_barcode_file.fasta
##		   (2) PURC
### Usage: PacBio_Seq_Barcodes.py PURC_barcode_file.fasta PURC_map_file.txt CCS_Sequence_File_1.fasta[.fastq] CCS_Sequence_File_2.fasta[.fastq]...


try:
	from Bio.Seq import Seq
	from Bio import SeqIO
	from Bio.Alphabet import generic_dna
except:
	print "ERROR: BioPython not installed"
import sys

### read in PURC Barcode File
purc_barcode_file = open(sys.argv[1], "r")
purcBarcodeDict = {}
purcBarcodeLength = []
for line in purc_barcode_file:
	if ">" in line:
		BCname = line.strip("\n")
	elif ">" not in line:
		purcBarcodeDict[line.strip("\n")] = BCname
		purcBarcodeLength.append(len(line.strip("\n")))
purc_barcode_file.close()

if len(set(purcBarcodeLength)) > 1:
	print "!" * 10
	print "ERROR: Barcodes not all same length"
	print "!" * 10
	sys.exit()
else:
	purcBarcodeLength = purcBarcodeLength[0]

### check for dupe barcodes (reverse complements)
for key in purcBarcodeDict:
	revComp = Seq(key, generic_dna)
	revComp = revComp.reverse_complement()
	if revComp in purcBarcodeDict.keys():
		print "WARNING: %s [%s] and %s [%s] are reverse complements!" %(purcBarcodeDict[key], key, purcBarcodeDict[revComp], revComp)
### read in PURC Map File
purc_map_file = open(sys.argv[2], "r")
purcMapDict = {}
for line in purc_map_file:
	splitline = line.strip("\n").split("\t")
	purcMapDict[splitline[2]] = [splitline[0],splitline[1]]
purc_map_file.close()

### Loop through input sequence files
for file in sys.argv[3:]:
	FbarcodeDict = {}
	RbarcodeDict = {}
	CombinedBarcodeDict = {}
	PairedBarcodeDict = {}
	fastaDict = {}
	totalseqs = 0
	key = 0
	
### parse fastq
	if ".fastq" in file or ".fq" in file:
		SeqIO.convert(file, "fastq", "PacBioBarcodes.temp.fasta", "fasta")
		infile = open("PacBioBarcodes.temp.fasta", "r")
### parse fasta
	elif ".fasta" in file or ".fa" in file:
		infile = open(file, "r")
### continue processing -- read fasta into dict to deal with interleaved format
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
			stripLine = line.strip("\n").upper()
			fastaDict[key].append(stripLine)
	joinLine = "".join(fastaDict[key])
	fastaDict[key] = joinLine
	print "%s seqs processed...DONE" %(totalseqs)
	
### parse fastaDict to read barcodes
	for key in fastaDict.keys():
		Fbarcode = fastaDict[key][0:purcBarcodeLength]
		endBC = fastaDict[key][-purcBarcodeLength:]
		endBC = Seq(endBC, generic_dna)
		Rbarcode = endBC.reverse_complement()
		Pairedbarcode = Fbarcode + Rbarcode
		try:
			FbarcodeDict[Fbarcode] += 1
		except:
			FbarcodeDict[Fbarcode] = 1
		try:
			RbarcodeDict[Rbarcode] += 1
		except:
			RbarcodeDict[Rbarcode] = 1
		try:
			CombinedBarcodeDict[Fbarcode] += 1
		except:
			CombinedBarcodeDict[Fbarcode] = 1
		try:
			CombinedBarcodeDict[Rbarcode] += 1
		except:
			CombinedBarcodeDict[Rbarcode] = 1
		try:
			PairedBarcodeDict[Pairedbarcode] += 1
		except:
			PairedBarcodeDict[Pairedbarcode] = 1
	
### output results
	Foutfile = open("Forward_PacBio_Barcodes.csv", "w")
	for key in FbarcodeDict.keys():
		Foutfile.write("%s,%s\n" %(key, FbarcodeDict[key]))
	Foutfile.close()
	
	Routfile = open("Reverse_PacBio_Barcodes.csv", "w")
	for key in RbarcodeDict.keys():
		Routfile.write("%s,%s\n" %(key, RbarcodeDict[key]))
	Routfile.close()

	Combinedoutfile = open("Combined_PacBio_Barcodes.csv", "w")
	for key in CombinedBarcodeDict.keys():
		Combinedoutfile.write("%s,%s\n" %(key, CombinedBarcodeDict[key]))
	Combinedoutfile.close()
	
	Pairedoutfile = open("Paired_PacBio_Barcodes.csv", "w")
	
### compare PURC barcodes to barcodes from CCS sequence file
	for key in PairedBarcodeDict.keys():
		if key[0:purcBarcodeLength] in purcBarcodeDict:
			FBCname = purcBarcodeDict[key[0:purcBarcodeLength]]
		else:
			try:
				FBCreversed = purcBarcodeDict[key[0:purcBarcodeLength]]
				FBCname = FBCreversed.reverse_complement() + "_reversed"
			except:
				FBCname = "N/A"
		if key[purcBarcodeLength:] in purcBarcodeDict.keys():
			RBCname = purcBarcodeDict[key[purcBarcodeLength:]]
		else:
			try:
				RBCreversed = purcBarcodeDict[key[purcBarcodeLength:]]
				RBCname = RBCreversed.reverse_complement() + "_reversed"
			except:
				RBCname = "N/A"
		sample = "No sample match!"
		for mapKey in purcMapDict.keys():
			if FBCname.strip(">") in purcMapDict[mapKey] and RBCname.strip(">") in purcMapDict[mapKey]:
				sample = mapKey
		Pairedoutfile.write("%s,%s,%s,%s,%s,%s\n" %(sample, key[0:purcBarcodeLength], key[purcBarcodeLength:], FBCname.strip(">"), RBCname.strip(">"), PairedBarcodeDict[key]))
	Pairedoutfile.close()
	
	
	
