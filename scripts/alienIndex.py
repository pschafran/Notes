#! /usr/bin/env python

import sys
import numpy as np

usageMsg='''alienIndex.py

Calculate alien index (AI) for a Diamond output file that includes taxonomy info. MUST create using output format 6 command:
--outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue evalue staxids sskingdoms skingdoms sphylums sscinames

Uses equation AI = (ln(bbhG + 1 * 10e-200)-ln(bbhO + 1 * 10e-200)) from Fan et al. (2020) Science Advances 6: eaba0111

AI ranges from approximately +/- 466, with AI >0 if evalue higher in ingroup, <0 if evalue higher in outgroup. Also reports percentage of hits that fell into ingroup/outgroup.'''

helpMsg='''
Required parameters
	--ingroup, -i	Name of taxonomic ingroup. Can be at any taxonomic level listed below. Set --taxon to search only a particular taxonomic level (use if ingroup name is shared among multiple taxonomic levels)
	--file, -f	Diamond output file with format specified above
Optional parameters
	--output, -o	Name of output file. If not specified, output is printed to terminal
	--taxon, -t	Taxonomic level for ingroup. Available options: 'superkingdom' , 'kingdom', 'phylum', 'genus'
	--missing, -m	How to treat N/A taxonomic annotations. Available options: 'outgroup' (default) or 'ingroup'.
	--help, -h	Display full usage
'''

# Parse command line and set vars

if "-h" in sys.argv or "--help" in sys.argv:
	print(usageMsg)
	print(helpMsg)
	exit(1)
if "-i" not in sys.argv and "--ingroup" not in sys.argv:
	print("ERROR: Ingroup not specified")
	print(helpMsg)
	exit(1)
if "-f" not in sys.argv and "--file" not in sys.argv:
	print("ERROR: BLAST results file not specified")
	print(helpMsg)
	exit(1)
if "-m" not in sys.argv and "--missing" not in sys.argv:
	missingData = "outgroup"
for item in sys.argv:
	if "-i" == item or "--ingroup" == item:
		ingroup = sys.argv[sys.argv.index(item)+1]
	if "-f" == item or "--file" == item:
		infile = sys.argv[sys.argv.index(item)+1]
	if item in ["-o", "--output", "--out", "-out"]:
		outfile = sys.argv[sys.argv.index(item)+1]
	if "-t" == item or "--taxon" == item:
		taxonRank = sys.argv[sys.argv.index(item)+1]
		if taxonRank not in ['superkingdom' , 'kingdom', 'phylum', 'genus']:
			print("ERROR: Not an accepted taxonomic level")
			print(helpMsg)
			exit(1)
	if "-m" == item or "--missing" == item:
		missingData = sys.argv[sys.argv.index(item)+1]
		if missingData != "outgroup" and missingData != "ingroup":
			print("ERROR: Not an accepted missing data option")
			print(helpMsg)
			exit(1)
# Check that ingroup is actually in file - warn against spelling errors
with open(infile) as openInfile:
	test_ingroup = []
	for line in openInfile:
		if 'taxonRank' in locals():
			if taxonRank == "superkingdom":
				for item in line.strip("\n").split("\t")[13].split(";"):
					test_ingroup.append(item)
			elif taxonRank == "kingdom":
				for item in line.strip("\n").split("\t")[14].split(";"):
					test_ingroup.append(item)
			elif taxonRank == "phylum":
				for item in line.strip("\n").split("\t")[15].split(";"):
					test_ingroup.append(item)
			elif taxonRank == "genus":
				test_ingroup.append(line.strip("\n").split("\t")[16].split(" ")[0].split(";")[0])
		else:
			for x in line.strip("\n").split("\t")[13:16]:
				test_ingroup.append(x)
			test_ingroup.append(line.strip("\n").split("\t")[16].split(" ")[0].split(";")[0])
	if ingroup.lower() not in [x.lower() for x in set(test_ingroup)]:
		print("*"*20)
		print("WARNING: Ingroup not found in file. Did you mispell or select wrong taxonomic level?")
		if 'taxonRank' in locals():
			print("Options for selected taxonomic level are: %s" % set(test_ingroup))
		print("*"*20)


# Run

mainDict = {} # Structure: { qseqid1 : {"bestBlastHitIngroup" = line}, {"bestBlastHitOutgroup" = line}, {"AI" : NUM} }, {"numIngroup" : int}, {"numOutgroup" = int} } ; qseqid2...}

with open(infile, 'r') as openInfile:
	for line in openInfile:
		qseqid = line.split("\t")[0]
		evalue = line.split("\t")[10]
		try: # needed to skip reassigning 0 to counts after first instance of qseqid
			mainDict[qseqid]["numIngroup"]
		except KeyError:
			mainDict[qseqid] = {"numIngroup" : float(0)}
			mainDict[qseqid].update({"numOutgroup" : float(0)})
		# set 'staxon' as the element to compare to ingroup name if --taxon is set
		if 'taxonRank' in locals():
			if taxonRank == "superkingdom":
				staxon = [x for x in line.strip("\n").split("\t")[13].split(";") if x != 0]
			elif taxonRank == "kingdom":
				staxon = [x for x in line.strip("\n").split("\t")[14].split(";") if x != 0]
			elif taxonRank == "phylum":
				staxon = [x for x in line.strip("\n").split("\t")[15].split(";") if x != 0]
			elif taxonRank == "genus":
				staxon = [line.strip("\n").split("\t")[16].split(" ")[0].split(";")[0]]
			try:
				staxon.remove("0")
			except:
				pass
		# set 'staxon' to all taxonomy fields if not specified in command line
		else:
			staxon = []
			for x in line.strip("\n").split("\t")[13:16]:
				try:
					templist = [y for y in x.split(";") if y != 0]
					for z in templist:
						staxon.append(z)
				except:
					staxon.append(x)
			staxon.append(line.strip("\n").split("\t")[16].split(" ")[0].split(";")[0])
		uniqstaxon = list(filter(None, [y for y in set(staxon)]))
		try:
			uniqstaxon.remove("0")
		except:
			pass
		# Compare evalue from current line to previous dictionary entry for ingroup/outgroup, replace with line if previous evalue lower than current
		if ingroup.lower() in [j.lower() for j in staxon] or uniqstaxon[0] == "N/A" and missingData == "ingroup":
			mainDict[qseqid]["numIngroup"] += 1
			try:
				previousBBHG = mainDict[qseqid]["bestBlastHitIngroup"].split("\t")[10] # retrieve evalue of best ingroup hit stored in dictionary for qseqid
				if evalue > previousBBHG:
					mainDict[qseqid]["bestBlastHitIngroup"] = line
			except KeyError: # if no previous entry in dictionary
				try:
					mainDict[qseqid].update({"bestBlastHitIngroup" : line})
				except:
					mainDict[qseqid] = {"bestBlastHitIngroup" : line}
		elif ingroup not in staxon or uniqstaxon[0] == "N/A" and missingData == "outgroup":
			mainDict[qseqid]["numOutgroup"] += 1
			try:
				previousBBHO = mainDict[qseqid]["bestBlastHitOutgroup"].split("\t")[10] # retrieve evalue of best outgroup hit stored in dictionary for qseqid
				if evalue > previousBBHO:
					mainDict[qseqid]["bestBlastHitOutgroup"] = line
			except KeyError: # if no previous entry in dictionary
				try:
					mainDict[qseqid].update({"bestBlastHitOutgroup" : line})
				except KeyError:
					mainDict[qseqid] = {"bestBlastHitOutgroup" : line}
		else:
			continue

		# Debugging precaution
		#else:
		#	print("This shouldn't happen - error comparing ingroup to staxon")
		#	print("Ingroup: %s , staxon: %s" %(ingroup, staxon))
		#	print(line)
		#	exit(1)

# Calculate AI and write to file or screen
if 'outfile' in locals():
	with open(outfile, 'w') as openOutfile:
		openOutfile.write("#QueryID\tNumber-Ingroup\tNumber-Outgroup\tPercent-Ingroup\tPercent-Outgroup\tBestBlastHit-Ingroup\tBestBlastHit-Outgroup\tAlienIndex\n")
		for qseqid in mainDict:
			percIngroup = 100 * mainDict[qseqid]["numIngroup"] / (mainDict[qseqid]["numIngroup"] + mainDict[qseqid]["numOutgroup"])
			percOutgroup = 100 * mainDict[qseqid]["numOutgroup"] / (mainDict[qseqid]["numIngroup"] + mainDict[qseqid]["numOutgroup"])
			try:
				bbhG = float(mainDict[qseqid]["bestBlastHitIngroup"].split("\t")[10])
			except KeyError:
				bbhG = 1
			try:
				bbhO = float(mainDict[qseqid]["bestBlastHitOutgroup"].split("\t")[10])
			except KeyError:
				bbhO = 1

			AI = (np.log(bbhG + (1*10e-200)))-(np.log(bbhO + (1*10e-200)))
			mainDict[qseqid].update({"AI" : AI}) # add AI to main dictionary in case I want to use later
			openOutfile.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" %(qseqid, mainDict[qseqid]["numIngroup"], mainDict[qseqid]["numOutgroup"], percIngroup, percOutgroup, bbhG, bbhO, AI))
else:
	print("#QueryID\tNumber-Ingroup\tNumber-Outgroup\tPercent-Ingroup\tPercent-Outgroup\tBestBlastHit-Ingroup\tBestBlastHit-Outgroup\tAlienIndex")
	for qseqid in mainDict:
		percIngroup = 100 * mainDict[qseqid]["numIngroup"] / (mainDict[qseqid]["numIngroup"] + mainDict[qseqid]["numOutgroup"])
		percOutgroup = 100 * mainDict[qseqid]["numOutgroup"] / (mainDict[qseqid]["numIngroup"] + mainDict[qseqid]["numOutgroup"])
		try:
			bbhG = float(mainDict[qseqid]["bestBlastHitIngroup"].split("\t")[10])
		except KeyError:
			bbhG = 1
		try:
			bbhO = float(mainDict[qseqid]["bestBlastHitOutgroup"].split("\t")[10])
		except KeyError:
			bbhO = 1

		AI = (np.log(bbhG + (1*10e-200)))-(np.log(bbhO + (1*10e-200)))
		mainDict[qseqid].update({"AI" : AI}) # add AI to main dictionary in case I want to use later
		try:
			print("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s" %(qseqid, mainDict[qseqid]["numIngroup"], mainDict[qseqid]["numOutgroup"], percIngroup, percOutgroup, bbhG, bbhO, AI))
		except:
			pass
