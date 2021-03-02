#! /usr/bin/python
'''Remove duplicate sequences by name and sequence. If duplicates found, new file output as sequences.dupSeqsCombined.fasta. If no duplicates found, no new file created.
Usage: removeDuplicateSeqFasta.py sequences.fasta
'''

import sys

nameDict = {}
seqDict = {}

dupnamecount = 0
with open(sys.argv[1], "r") as infile:
	for line in infile:
		if line.startswith(">"):
			duplicate = 0
			lineName = line.strip(">\n").split(" ")[0]
			if lineName in nameDict:
				print("Duplicate name: %s" % lineName)
				dupnamecount += 1
				duplicate = 1
			else:
				nameDict.update({lineName : []})
			try:
				seq = "".join(nameDict[previousLineName])
				nameDict[previousLineName] = seq.strip("*")
			except:
				pass
		else:
			if duplicate == 0:
				previousLineName = lineName
				nameDict[lineName].append(line.strip("\n"))
			else:
				previousLineName = lineName
	seq = "".join(nameDict[previousLineName])
	nameDict[previousLineName] = seq

for name in nameDict:
	seq = nameDict[name]
	if seq in seqDict:
		seqDict[seq].append(name)
	else:
		seqDict.update({seq : [name]})
fileNoExt = ".".join(sys.argv[1].split(".")[:-1])

multiseqcount = 0
for seq in seqDict:
	if len(seqDict[seq]) > 1:
		multiseqcount += 1
#if multiseqcount == 0 and dupnamecount == 0:
#	exit(0)
#else:
outfileName =  fileNoExt + ".dupSeqsCombined.fasta"
with open(outfileName, "w") as outfile:
	for seq in seqDict:
		if len(seqDict[seq]) > 1:
			csvNames = ",".join(seqDict[seq])
			joinNames = "---".join(seqDict[seq])
			print("Duplicate sequences: %s\t%s" %(len(seqDict[seq]), csvNames))
			outfile.write(">%s\n" % joinNames)
			outfile.write("%s\n" % seq)
		else:
			outfile.write(">%s\n" % seqDict[seq][0])
			outfile.write("%s\n" % seq)
