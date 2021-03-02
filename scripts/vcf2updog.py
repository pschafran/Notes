#! /usr/bin/env python

# Usage: vcf2updog.py ipyrad.vcf (optional) min_num_samples_with_data

import sys

vcfFile = sys.argv[1] 
minSamples = 0
try:
	minSamples = int(sys.argv[2])
except:
	pass
vcfDict = {}

with open(vcfFile, "r") as infile:
	snpNum = 1
	for line in infile:
		if line.startswith("#"):
			if line.startswith("#CHROM"):
				splitline = line.strip("\n").split("\t")
				if splitline[7] != "INFO":
					print("Formatting error -- is this an ipyrad file?")
					exit(1)
				sampleList = splitline[9:]
		else:
			splitline = line.strip("\n").split("\t")
			if int(splitline[7].split(";")[0].split("=")[1]) >= minSamples:
				snpName = "SNP%s" % snpNum
				vcfDict[snpName] = {"refAllele" : "N"}
				vcfDict[snpName].update({"altAllele" : "N"})
				#vcfDict[snpName].update({"sampleDict" : {} })
				vcfDict[snpName].update({"refCountList" : []})
				vcfDict[snpName].update({"altCountList" : []})
				vcfDict[snpName].update({"totalCountList" : []})
				splitline = line.strip("\n").split("\t")
				refAllele = splitline[3]
				altAllele = splitline[4]
				samples = splitline[9:]
				vcfDict[snpName]["refAllele"] = refAllele
				vcfDict[snpName]["altAllele"] = altAllele
				if refAllele == "C":
					refAlleleField = 0
				elif refAllele == "A":
					refAlleleField = 1
				elif refAllele == "T":
					refAlleleField = 2
				elif refAllele == "G":
					refAlleleField = 3
				if altAllele == "C":
					altAlleleField = 0
				elif altAllele == "A":
					altAlleleField = 1
				elif altAllele == "T":
					altAlleleField = 2
				elif altAllele == "G":
					altAlleleField = 3
				for sample in samples:
					sampleField = sample.split(":")
					sampleAlleleCounts = sampleField[2].split(",")
					vcfDict[snpName]["refCountList"].append(sampleAlleleCounts[refAlleleField])
					vcfDict[snpName]["altCountList"].append(sampleAlleleCounts[altAlleleField])
					vcfDict[snpName]["totalCountList"].append(sampleField[1])
		snpNum += 1

with open("%s_updog_reference_matrix.tsv" % vcfFile, "w") as outfile:
	outfile.write("\t%s\n" % "\t".join(sampleList))
	dictKeys = vcfDict.keys()
	for snpName in dictKeys:
		outfile.write("%s\t%s\n" %(snpName, "\t".join(vcfDict[snpName]["refCountList"])))
		
with open("%s_updog_alternate_matrix.tsv" % vcfFile, "w") as outfile:
	outfile.write("\t%s\n" % "\t".join(sampleList))
	dictKeys = vcfDict.keys()
	for snpName in dictKeys:
		outfile.write("%s\t%s\n" %(snpName, "\t".join(vcfDict[snpName]["altCountList"])))

with open("%s_updog_readDepth_matrix.tsv" % vcfFile, "w") as outfile:
	outfile.write("\t%s\n" % "\t".join(sampleList))
	dictKeys = vcfDict.keys()
	for snpName in dictKeys:
		outfile.write("%s\t%s\n" %(snpName, "\t".join(vcfDict[snpName]["totalCountList"])))