#!/usr/bin/env python

# Tree tip labels expected to be genus_species_...etc. Labels will be split on first underscore and only genus used for searching taxonomic hierarchy

from ete3 import Tree
import glob
import sys
import re

try:
	tree_list = sys.argv[1:]
except:
	tree_list = glob.glob('*contree')

focal_taxon = 'Ceratopteris_richardii' # tip starts with this string (no splitting into genus-species)
exclude_list = ['Eukaryota'] # Taxonomic group(s) to be considered focal groups
donor_list = ['Bacteria', 'Archaea', 'Viruses'] # Taxonomic group(s) expected to be outgroup (i.e. HGT donors)
excluded_sample_list = []
donor_sample_list = []
missing_sample_list = []
error_tree_list = []
too_few_tips_list = []
no_focal_tips_list = []
all_bact_list = []
nest_in_bact_list = []
euk_sister_list = []
mix_sister_list = []

# function to parse NCBI taxonomy searching upward until the queryName (or its higher taxonomic level) is found in the excludeList or donorList (or root of all life is reached)
def parseTaxonomy(queryName, excludeList, donorList):
	currentName = queryName
	keyFound = False
	donorCount = 0
	excludeCount = 0
	missingCount = 0
	focalTipsPresent = False
	for key in nodesDict:
		if nodesDict[key][1] == currentName:
			currentNode = key
			currentRank = nodesDict[key][0]
			currentName = nodesDict[key][1]
			parentNode = nodesDict[key][2]
			keyFound = True
			break
	if keyFound == False:
		currentNode = 1
		currentRank = "root"
		currentName = "root"
		parentNode = 1
	while currentName not in donorList and currentName not in excludeList:
		try:
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
	if currentName in donorList:
		donorCount += 1
		global donor_sample_list
		donor_sample_list.append(queryName)
	elif currentName in excludeList:
		excludeCount += 1
		focalTipsPresent = True
		global excluded_sample_list
		excluded_sample_list.append(queryName)
	else:
		missingCount += 1
		global missing_sample_list
		missing_sample_list.append(queryName)
	return([queryName, donorCount, excludeCount, missingCount, focalTipsPresent])

# Create dictionary of taxonmic heirarchy
nodesDB = "/home/ps997/bin/blobtools/data/nodesDB.txt"
lineCount = 0
nodesDict = {}
with open(nodesDB, "r") as openNodesDB:
	for line in openNodesDB:
		lineCount += 1
		splitline = line.strip("\n").split("\t")
		if lineCount > 1:
			try:
				nodesDict[splitline[0]] = [splitline[1],splitline[2],splitline[3]]
			except:
				print("ERROR on line %s of %s: Incorrectly formatted,ust have 4 tab-separated columns" %(lineCount, nodesDB))
				exit(1)


for tree in tree_list:
	try:
		t = Tree(tree, format=0)
	except:
		print('error reading ', tree)
		error_tree_list.append(tree)
		continue

	### Get basic info from the tree ###
	root = t.get_tree_root()
	all_leaf_name = root.get_leaf_names()
	total_leaf_number = len(all_leaf_name)

	### Criterion 1 ###
	### have at least 5 tips ###
	if total_leaf_number < 5:
		too_few_tips_list.append(tree)
		continue

	### Criterion 2 ###
	### tips are all bacteria ###
	bacteria_count = 0
	total_count = 0
	missing_count = 0
	donor_count = 0
	focal_tips_present = False
	for k in all_leaf_name:
		k2 = k.strip(">\n").split("_")[0]
		queryName = k2.split("-")[0]
		results = parseTaxonomy(queryName, exclude_list, donor_list)
		bacteria_count += results[1]
		donor_count += results[2]
		missing_count += results[3]
		if results[4] == True:
			focal_tips_present = True
		total_count += 1
	#	if k.startswith(tuple(donor_list)):
	#		bacteria_count = bacteria_count + 1
	#		total_count = total_count + 1
	#	elif k.startswith(tuple(exclude_list)):
	#		focal_tips_present = True
	#	else:
	#		total_count = total_count + 1
	if focal_tips_present == True:
		if bacteria_count > 0 and total_count > 0:
			if bacteria_count/float(total_count) == 1.0:
				all_bact_list.append(tree)
				focal_leaves = [i for i in all_leaf_name if i.startswith(focal_taxon)]
				for leaf in focal_leaves:
					print(leaf, '\tall_bact')
				continue
	elif focal_tips_present == False:
		no_focal_tips_list.append(tree)
		continue

	### Criterion 3 ###
	### focal tip nested within bacteria ###
	focal_leaves = [i for i in all_leaf_name if i.startswith(focal_taxon)]
	for focal_leaf in focal_leaves:
		node = t.search_nodes(name=focal_leaf)[0]

		# reroot the tree with the furthest node from the focal tip
		node_furthest = node.get_farthest_node()
		t.set_outgroup(node_furthest[0])

		# loop through the list of tips, traverse toward the root
		pass_number = 0 # control the number of nodes to traverse
		with_high_support_node = False # at least one node needs to be highly supported
		while node and pass_number <= 2:
			des_leaves = node.get_leaf_names() # all the descendent leaves from this node

			# check if the node is supported, one supported node will satisfy
			if float(node.support) > 90:
				with_high_support_node = True

			# get sister tips
			sister_node = node.get_sisters() # get the sister node (as a list oddly)
			sister_node_leaves = []
			for i in sister_node: # loop through the sister node list and get tips
				for j in i.get_leaf_names():
					sister_node_leaves.append(j)

			# check if the sister tips contain any of the eustig seqs, defined by exclude_list
			pass_increase = True
			for des_leaf in sister_node_leaves:
				des_leaf2 = des_leaf.strip(">\n").split("_")[0]
				des_leaf_query = des_leaf2.split("-")[0]
				results = parseTaxonomy(des_leaf_query, exclude_list, donor_list)
				if results[4] == True:
					pass_increase = False
			if pass_increase and with_high_support_node:
				pass_number = pass_number + 1

			node_previous = node
			node = node.up

		# check if all bac sisters
		bacteria_count = 0
		total_count = 0
		for k in des_leaves:
			k2 = k.strip(">\n").split("_")[0]
			queryName = k2.split("-")[0]
			results = parseTaxonomy(queryName, exclude_list, donor_list)
			if results[4] == False:
				bacteria_count = bacteria_count + 1
				total_count = total_count + 1
			elif results[4] == True:
				pass
			else:
				total_count = total_count + 1
		if bacteria_count > 0 and total_count > 0:
			if bacteria_count/float(total_count) == 1.0:
				nest_in_bact_list.append(tree)
				print(focal_leaf, '\tnested_in_bact')
				#print('#', tree, total_leaf_number, focal_leaf)
				#print(node_previous)
			else:
				mix_sister_list.append(tree)
		else:
			euk_sister_list.append(tree)
			#print('##', tree, total_leaf_number, focal_leaf)
			#print(node_previous)

print('trees analyzed: ', len(tree_list))
print('error_tree_list', len(error_tree_list))
print('too_few_tips_list', len(too_few_tips_list))
print('no_focal_tips_list', len(no_focal_tips_list))
print('all_bact_list', len(all_bact_list))
print('nest_in_bact_list', len(nest_in_bact_list))
print('euk_sister_list', len(euk_sister_list))
print('mix_sister_list', len(mix_sister_list))

with open("traverse_tree.out", "w") as outfile:
	outfile.write("Trees analyzed: %s\n" % len(tree_list))
	outfile.write("Euk_sister trees: %s\n" % len(euk_sister_list))
	for ES_tree in euk_sister_list:
		outfile.write("\t%s\n" % ES_tree)
	outfile.write("\n")
	outfile.write("Mix_sister trees: %s\n" % len(mix_sister_list))
	for MS_tree in mix_sister_list:
		outfile.write("\t%s\n" % MS_tree)
	outfile.write("\n")
	outfile.write("Nest_in_bact trees: %s\n" % len(nest_in_bact_list))
	for NIB_tree in nest_in_bact_list:
		outfile.write("\t%s\n" % NIB_tree)
	outfile.write("\n")
	outfile.write("All_bact trees: %s\n" % len(all_bact_list))
	for AB_tree in all_bact_list:
		outfile.write("\t%s\n" % AB_tree)
	outfile.write("\n")
	outfile.write("No_focal_tips trees: %s\n" % len(no_focal_tips_list))
	for NFT_tree in no_focal_tips_list:
		outfile.write("\t%s\n" % NFT_tree)
	outfile.write("\n")
	outfile.write("Too_few_tips trees: %s\n" % len(too_few_tips_list))
	for TFT_tree in too_few_tips_list:
		outfile.write("\t%s\n" % TFT_tree)
	outfile.write("\n")
	outfile.write("Error trees: %s\n" % len(error_tree_list))
	for E_tree in error_tree_list:
		outfile.write("\t%s\n" % E_tree)
