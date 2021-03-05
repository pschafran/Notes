#! /usr/bin/env python

import sys
import numpy as np
import time

usageMsg='''alienIndex.py

Calculate alien index (AI) for a Diamond output file that includes taxonomy info. MUST create using output format 6 command:
--outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue evalue staxids sskingdoms skingdoms sphylums sscinames

Uses equation AI = (ln(bbhG + 1 * 10e-200)-ln(bbhO + 1 * 10e-200)) from Fan et al. (2020) Science Advances 6: eaba0111

AI ranges from approximately +/- 466, with AI > 0 if evalue higher in ingroup, < 0 if evalue higher in outgroup (So AI > 0 means query is more similar to outgroup than ingroup). Also reports raw numbers and percentage of hits that fell into ingroup/outgroup.'''

helpMsg='''
Required parameters
	--ingroup, -i	Name of taxonomic ingroup. Can be at any taxonomic level listed below. Set --taxon to search only a particular taxonomic level (use if ingroup name is shared among multiple taxonomic levels)
	--file, -f	Diamond output file with format specified above
Optional parameters
	--outgroup, -g	Specify outgroup. If not set, all sequences not in the ingroup are treated as outgroup. Set this if you want to leave a 'gap' between ingroup and outgroup.
	--ignore, -n	Specify clade within the ingroup to ignore. E.g. Ingroup is Eukaryota, but do not consider hits to Polypodiopsida.
	--output, -o	Name of output file. If not specified, output is printed to stdout
	--missing, -m	How to treat N/A taxonomic annotations. Available options: 'ignore' (default), 'outgroup', 'ingroup.
	--help, -h	Display full usage
	--log, -l	NOT ACTIVE. File to write log containing warnings generated and offending lines in BLAST file for debugging. Does not affect what is printed to STDOUT.
	--taxon, -t	Usually not necessary to set. Taxonomic level for ingroup. Available options: 'superkingdom' , 'kingdom', 'phylum', 'genus'

'''

# function to parse dictionary of NCBI taxonomy searching upward until the queryName (or its higher taxonomic level) is found in the outgroupList or ingroupList (or root of all life is reached)
def parseTaxonomy(queryName, outgroupList, ingroupList, ignoregroupList):
	currentName = queryName
	keyFound = False
	ingroupCount = 0
	outgroupCount = 0
	missingCount = 0
	ignoregroupCount = 0
	try:
		currentNode = firstNodesDict[queryName][0]
		currentRank = firstNodesDict[queryName][1]
		parentNode = firstNodesDict[queryName][2]
		keyFound = True
	except KeyError:
		keyfound = False
	#for key in nodesDict:
	#	if nodesDict[key][1] == currentName:
	#		currentNode = key
	#		currentRank = nodesDict[key][0]
	#		currentName = nodesDict[key][1]
	#		parentNode = nodesDict[key][2]
	#		keyFound = True
	#		break
	if keyFound == False:
		currentNode = 1
		currentRank = "root"
		currentName = "root"
		parentNode = 1
		warning = "WARNING: %s not found in NCBI taxonomy" % (queryName)
		print(warning, file=sys.stderr)
		if 'logfile' in globals():
			global logfile
			logfile.write("WARNING: %s not found in NCBI taxonomy\n" % queryName)
	while currentName not in ingroupList and currentName not in outgroupList and currentName not in ignoregroupList:
		try:
			#print(currentName, currentRank)
			currentNode = parentNode
			currentRank = nodesDict[parentNode][0]
			currentName = nodesDict[parentNode][1]
			parentNode = nodesDict[parentNode][2]
			if currentNode == 1:
				#print("Reached root without finding corrent taxonomic rank" %(currentNode, parentNode))
				#currentRank = taxonGrouping
				break
			elif currentNode != 1 and currentNode == parentNode:
				#print("Error: current node and parent node are the same: %s\t%s" %(currentNode, parentNode))
				#currentRank = taxonGrouping
				break
		except KeyError:
			break
	if currentName in ingroupList:
		ingroupCount += 1
	elif currentName in outgroupList or len(outgroupList) == 0 and currentNode == 1:
		outgroupCount += 1
	elif currentName in ignoregroupList:
		ignoregroupCount += 1
	else:
		missingCount += 1
	return([queryName, ingroupCount, outgroupCount, ignoregroupCount, missingCount]) # focalTipsPresent probably unecessary and can be removed b/c it can be inferred from 0/1 in other results

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
		ingroupList = [ingroup]
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
		if missingData != "outgroup" and missingData != "ingroup" and missingData != "ignore":
			print("ERROR: Not an accepted missing data option")
			print(helpMsg)
			exit(1)
	if "-g" == item or "--outgroup" == item:
		outgroupList = sys.argv[sys.argv.index(item)+1].split(",")
	if "-n" == item or "--ignore" == item:
		ignoregroupList = sys.argv[sys.argv.index(item)+1].split(",")
	if "-l" == item or "--log" == item:
		logfileName = sys.argv[sys.argv.index(item)+1]
		logfile = open(logfileName,"w")

if 'outgroupList' not in locals():
	outgroupList = []
if 'ignoregroupList' not in locals():
	ignoregroupList = []

# Check that ingroup is actually in file - warn against spelling errors. Initiate extensive taxon search if any of the specified ingroup/outgroup/ignoregroup are missing from BLAST file
extensiveTaxonSearch = False
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
	for ingroup in ingroupList:
		if ingroup.lower() not in [x.lower() for x in set(test_ingroup)]:
			extensiveTaxonSearch = True
	if 'outgroupList' in locals():
		for outgroup in outgroupList:
			if outgroup.lower() not in [x.lower() for x in set(test_ingroup)]:
				extensiveTaxonSearch = True
	if 'ignoregroupList' in locals():
		for ignoregroup in ignoregroupList:
			if ignoregroup.lower() not in [x.lower() for x in set(test_ingroup)]:
				extensiveTaxonSearch = True
	if extensiveTaxonSearch == True:
		print("WARNING: Ingroup, outgroup, or ignoregroup not found in file. Will perform extensive (slow) taxon search.", file=sys.stderr)
		if 'taxonRank' in locals():
			print("Options for selected taxonomic level are: %s" % set(test_ingroup), file=sys.stderr)

# Create dictionary of taxonomic heirarchy (only if identifed as needed above)
if extensiveTaxonSearch == True:
	nodesDB = "/home/ps997/bin/blobtools/data/nodesDB.txt"
	lineCount = 0
	nodesDictNamesList = []
	nodesDict = {}
	with open(nodesDB, "r") as openNodesDB:
		for line in openNodesDB:
			lineCount += 1
			splitline = line.strip("\n").split("\t")
			if lineCount > 1:
				try:
					nodesDict[splitline[0]] = [splitline[1],splitline[2],splitline[3]]
					nodesDictNamesList.append(splitline[2])
				except:
					print("ERROR on line %s of %s: Incorrectly formatted, must have 4 tab-separated columns" %(lineCount, nodesDB))
					exit(1)
	firstNodesDict = {}
	lineCount = 0
	with open(nodesDB, "r") as openNodesDB:
		for line in openNodesDB:
			lineCount += 1
			splitline = line.strip("\n").split("\t")
			if lineCount > 1:
				try:
					firstNodesDict[splitline[2]] = [splitline[0],splitline[1],splitline[3]]
				except:
					print("ERROR on line %s of %s: Incorrectly formatted, must have 4 tab-separated columns" %(lineCount, nodesDB))
					exit(1)


####### Run #######
# Parse BLAST result file. In each line, use taxonomy info in columns 13,14,15,16 to determine if line is member of ingroup, outgroup, or neither.
# Creates dictionary with lowest evalue for ingroup and outgroup, and numbers of hits in ingroup and outgroup for each query sequence.
mainDict = {} # Structure: { qseqid1 : {"bestBlastHitIngroup" = line}, {"bestBlastHitOutgroup" = line}, {"AI" : NUM} }, {"numIngroup" : int}, {"numOutgroup" = int} } ; qseqid2...}
totalLines = 0
with open(infile, 'r') as openInfile:
	for line in openInfile:
		totalLines += 1
lineCount = 0
startTime = time.time()
with open(infile, 'r') as openInfile:
	for line in openInfile:
		lineCount += 1
		percDone = float(lineCount/totalLines*100)
		currentTime = time.time()
		runTime = currentTime - startTime
		linesPerSec = lineCount/runTime
		sys.stderr.write("\r")
		sys.stderr.write("Percent Complete: %0.2f\t\tAvg. speed: %0.2f lines/sec" % (percDone, linesPerSec))
		sys.stderr.flush()
		qseqid = line.split("\t")[0]
		evalue = line.split("\t")[10]
		try: # needed to skip reassigning 0 to counts after first instance of qseqid
			mainDict[qseqid]["numIngroup"]
		except KeyError:
			mainDict[qseqid] = {"numIngroup" : float(0)}
			mainDict[qseqid].update({"numOutgroup" : float(0)})

		# subsection if extensive taxon search needed
		if extensiveTaxonSearch == True:
			# Start search with genus from BLAST line
			staxon = [line.strip("\n").split("\t")[16].split(" ")[0].split(";")[0]]
			# Exception to handle Candidatus names because they have a space
			if "Candidatus" in staxon:
				staxon = [" ".join(line.strip("\n").split("\t")[16].split(" ")[0:2])]
			try:
				staxon.remove("0")
			except:
				pass
			uniqstaxon = list(filter(None, [y for y in set(staxon)]))
			try:
				uniqstaxon.remove("0")
			except:
				pass
			ingroupFound = False
			outgroupFound = False
			ignoregroupFound = False
			missingGroup = False
			# If an invalid genus name found, try the phylum instead
			#if uniqstaxon[0] not in nodesDictNamesList:
			#if uniqstaxon[0] == "N/A" or uniqstaxon[0] == "synthetic" or uniqstaxon[0] == "unclassified":
			#	staxon = [x for x in line.strip("\n").split("\t")[15].split(";") if x != 0]
			#	try:
			#		staxon.remove("0")
			#	except:
			#		pass
			#	uniqstaxon = list(filter(None, [y for y in set(staxon)]))
			#	try:
			#		uniqstaxon.remove("0")
			#	except:
			#		pass
			# As long as taxon name is not N/A, try parse taxonomy to find the name to determine if ingroup, outgroup, or ignored
			if uniqstaxon[0] != "N/A":
				for taxon in uniqstaxon:
					# parseTaxonomy expects lists for ingroup, outgroup so have to create them. Should probably be removed and just use strings for this script.
					if taxon in nodesDictNamesList:
						results = parseTaxonomy(taxon, outgroupList, ingroupList, ignoregroupList)
						# results == [queryName, ingroupCount, outgroupCount, missingCount, focalTipsPresent]
						sumResults = results[1] + results[2] + results[3] + results[4]
						if results[1] == 1:
							ingroupFound = True
						elif results[2] == 1:
							outgroupFound = True
						elif results[3] == 1:
							ignoregroupFound = True
						elif results[4] == 1:
							missingGroup = True
					else:
						staxon2 = [x for x in line.strip("\n").split("\t")[15].split(";") if x != 0]
						try:
							staxon2.remove("0")
						except:
							pass
						uniqstaxon2 = list(filter(None, [y for y in set(staxon2)]))
						try:
							uniqstaxon2.remove("0")
						except:
							pass
						for taxon2 in uniqstaxon2:
							if taxon2 in nodesDictNamesList:
								results2 = parseTaxonomy(taxon2, outgroupList, ingroupList, ignoregroupList)
								sumResults = results2[1] + results2[2] + results2[3] + results2[4]
								if results2[1] == 1:
									ingroupFound = True
								elif results2[2] == 1:
									outgroupFound = True
								elif results2[3] == 1:
									ignoregroupFound = True
								elif results[4] == 1:
									missingGroup = True
			if ingroupFound == False and outgroupFound == False: # don't do anything, assuming query falls in gap between ingroup and outgroup or in ignoregroup
				pass
			elif ingroupFound == True and outgroupFound == False  and ignoregroupFound == False or uniqstaxon[0] == "N/A" and missingData == "ingroup":
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
			elif ingroupFound == False and outgroupFound == True and ignoregroupFound == False or uniqstaxon[0] == "N/A" and missingData == "outgroup":
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
			elif sumResults > 1: # catch errors where both multiple groups are found in the same BLAST line
				print("WARNING: Multiple group names found. Check your groups are mutually exclusive in NCBI taxonomy. Offending taxon list: %s" % uniqstaxon , file=sys.stderr)

		# subsection if only standard taxon searching needed
		else:
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
			if ingroup.lower() in [j.lower() for j in staxon] and len(set(staxon)&set(ignoregroupList)) == 0 or uniqstaxon[0] == "N/A" and missingData == "ingroup":
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
			elif len(outgroupList) >= 1:
				if len(set(staxon)&set(outgroupList)) >= 1 or uniqstaxon[0] == "N/A" and missingData == "outgroup":
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
			elif len(outgroupList) == 0:
				if ingroup.lower() not in [j.lower() for j in staxon] or uniqstaxon[0] == "N/A" and missingData == "outgroup":
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
			numIngroup = mainDict[qseqid]["numIngroup"]
			numOutgroup = mainDict[qseqid]["numOutgroup"]
			if numIngroup == 0 and numOutgroup == 0:
				percIngroup = 0
				percOutgroup = 0
			else:
				percIngroup = 100 * numIngroup / (numIngroup + numOutgroup)
				percOutgroup = 100 * numOutgroup / (numIngroup + numOutgroup)
			try:
				bbhG = float(mainDict[qseqid]["bestBlastHitIngroup"].split("\t")[10])
				#taxG = mainDict[qseqid]["bestBlastHitIngroup"].split("\t")[13:16]
			except KeyError:
				bbhG = 1
				#taxG = "No hit"
			try:
				bbhO = float(mainDict[qseqid]["bestBlastHitOutgroup"].split("\t")[10])
				#taxO = mainDict[qseqid]["bestBlastHitOutgroup"].split("\t")[13:16]
			except KeyError:
				bbhO = 1
				#taxO = "No hit"

			AI = (np.log(bbhG + (1*10e-200)))-(np.log(bbhO + (1*10e-200)))
			mainDict[qseqid].update({"AI" : AI}) # add AI to main dictionary in case I want to use later
			openOutfile.write("%s\t%d\t%d\t%0.2f\t%0.2f\t%s\t%s\t%0.2f\n" %(qseqid, int(mainDict[qseqid]["numIngroup"]), int(mainDict[qseqid]["numOutgroup"]), percIngroup, percOutgroup, bbhG, bbhO, AI))
else:
	print("#QueryID\tNumber-Ingroup\tNumber-Outgroup\tPercent-Ingroup\tPercent-Outgroup\tBestBlastHit-Ingroup\tBestBlastHit-Outgroup\tAlienIndex")
	for qseqid in mainDict:
		numIngroup = mainDict[qseqid]["numIngroup"]
		numOutgroup = mainDict[qseqid]["numOutgroup"]
		if numIngroup == 0 and numOutgroup == 0:
			percIngroup = 0
			percOutgroup = 0
		else:
			percIngroup = 100 * numIngroup / (numIngroup + numOutgroup)
			percOutgroup = 100 * numOutgroup / (numIngroup + numOutgroup)
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
			print("%s\t%d\t%d\t%0.2f\t%0.2f\t%s\t%s\t%0.2f" %(qseqid, int(mainDict[qseqid]["numIngroup"]), int(mainDict[qseqid]["numOutgroup"]), percIngroup, percOutgroup, bbhG, bbhO, AI))
		except:
			pass
if 'logfileName' in locals():
	logfile.close()
