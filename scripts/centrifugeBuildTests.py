#! /usr/bin/python

# centrifugeBuildTests.py seq.fna taxid.map

import sys

fasta = open(sys.argv[1],"r")
taxIDmap = open(sys.argv[2],"r")
outfile = open("centrifuge_build_tests.out","w")

taxIDlist = []
seqNameList = []
testsFailed = 0
for line in taxIDmap:
	taxIDlist.append(line.split("\t")[0])
	if len(line.split("\t")) > 2:
		print ("%s contains more than 2 fields" % line)
		outfile.write("%s contains more than 2 fields" % line)
		testsFailed += 1
for seq in fasta:
	if seq.startswith(">"):
		if seq.strip(">\n").split(" ")[0] not in taxIDlist:
			print("%s found in sequence file but missing from tax ID map" % seq.strip(">\n"))
			outfile.write("%s found in sequence file but missing from tax ID map" % seq.strip(">\n"))
			testsFailed += 1
		seqNameList.append(seq.strip(">\n").split(" ")[0])
for item in taxIDlist:
	if item not in seqNameList:
		print("%s found in tax ID map but missing from sequence file" % item )
		outfile.write("%s found in tax ID map but missing from sequence file" % item )
		testsFailed += 1
if len(taxIDlist) != len(seqNameList):
	print("Tax ID map and sequence file have difference number of entries")
	print("Tax ID map: %d" % len(taxIDlist))
	print("Sequence file: %d" % len(seqNameList))
	outfile.write("Tax ID map and sequence file have difference number of entries")
	outfile.write("Tax ID map: %d" % len(taxIDlist))
	outfile.write("Sequence file: %d" % len(seqNameList))
	testsFailed += 1
if testsFailed == 0:
	print("All tests passed!")
	outfile.write("All tests passed!")
else:
	print("Failed %d tests" % testsFailed)
	outfile.write("Failed %d tests" % testsFailed)

fasta.close()
taxIDmap.close()
outfile.close()
