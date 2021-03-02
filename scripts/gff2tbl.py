#! /usr/bin/python
##2017 Nov 7 (C) Peter Schafran

##Converts GFF3 annotation file from Geneious to GenBank 5-column annotation file
##WARNINGS: Assumes all annotations uniquely named -- will join annotations with the same name (see note below)
##          Requires user to manually edit output file to fix input name, truncated sequences
## 8 Nov 2017 -- modified to write individual filenames and table headers as GenBank accessions; creates single merged file at end
## 27 Nov 2017 -- modified to deal with inverted repeats of CPgenome (i.e. regions with same name but in opposite directions should not be joined)
## 28 Nov 2017 -- added automatic protein_id numbering, unique for each run of gff2tbl.py (i.e. output from multiple input files will get unique IDs)
## 30 Nov 2017 -- added function to create .pep file with protein_id, gene, and protein names when an "RNA editing" exception is encountered

import sys, re

filelist=sys.argv[1:]
newfilelist = []
protein_id_num = int(input("Protein ID starting number: "))
for gff_file in filelist:
	outputfilename=gff_file.split('.')
	infile = open(gff_file, 'r')
	#Create output tbl file and write first header line
	outfile = open("%s.tbl" %(outputfilename[0]), 'w')
	pepoutfile = open("%s.pep" %(outputfilename[0]), 'w')
	newfilelist.append("%s.tbl" %(outputfilename[0]))
	outfile.write(">Feature\t%s\tTable\n" %(outputfilename[0]))
	line_count=1
	gff_dict={}
	feature_dict={}
	#Parse input GFF3 file. Splits on tabs and enters this into gff_dict
	for line in infile:
		splitline=line.strip('\n').split('\t')
		gff_dict[line_count]=splitline
		#This if line skips the header in GFF3 files
		if line.startswith("#"):
			pass
		else:
			try:
				#print "Processing line %d - %s" %(line_count, line)
				#Breaks up items in descriptor column of GFF3 file
				splitline_2 = gff_dict[line_count][8].split(";")
				gff_dict[line_count][8]=splitline_2
				feature_dict[line_count]={}
				for item in gff_dict[line_count][8]:
					#print item
					splititem=item.strip('\n').split('=')
					feature_dict[line_count][splititem[0]]=splititem[1]
				
				#write elements from dictionaries to output tbl file
				#First two if lines determine if annotation is for forward or reverse strand and switches start/end location for reverse. Writes line with start/end points and NCBI feature class
				#print gff_dict[line_count]
				#print feature_dict[line_count]
	#			if gff_dict[line_count][6] == "+":
	#				outfile.write("%s\t%s\t%s\n" %(gff_dict[line_count][3],gff_dict[line_count][4],gff_dict[line_count][2]))
	#			if gff_dict[line_count][6] == "-":
	#				outfile.write("%s\t%s\t%s\n" %(gff_dict[line_count][4],gff_dict[line_count][3],gff_dict[line_count][2]))
	#			if "gene" in feature_dict[line_count].keys():
	#				outfile.write("\t\t\tgene\t%s\n" %(feature_dict[line_count]["gene"]))
	#			if "product" in feature_dict[line_count].keys():
	#				outfile.write("\t\t\tproduct\t%s\n" %(feature_dict[line_count]["product"]))
	#			if "note" in feature_dict[line_count].keys():
	#				outfile.write("\t\t\tnote\t%s\n" %(feature_dict[line_count]["note"]))
	#			if "protein_id" in feature_dict[line_count].keys():
	#				outfile.write("\t\t\tprotein_id\t%s\n" %(feature_dict[line_count]["protein_id"]))
	#			if "codon_start" in feature_dict[line_count].keys():
	#				outfile.write("\t\t\tcodon_start\t%s\n" %(feature_dict[line_count]["codon_start"]))
	#			if "exception" in feature_dict[line_count].keys():
	#				outfile.write("\t\t\texception\t%s\n" %(feature_dict[line_count]["exception"]))
	#			if "Name" in feature_dict[line_count].keys():
	#				outfile.write("\t\t\tname\t%s\n" %(feature_dict[line_count]["Name"]))
			except:
				#print "Cannot process line #%d - %s" %(line_count, line)
				pass
		line_count+=1
	infile.close()
	#print feature_dict
	#print "%s\n" %(gff_dict)
	features = []
	for key in feature_dict.keys():
		
		try:
			features.append(feature_dict[key]["Name"])
		except:
			pass
	unique_features = list(set(features))
##Converts feature types to case that tbl2asn likes
	for gffline in gff_dict:
		try:
			if str(gff_dict[gffline][2]).upper() == "CDS":
				gff_dict[gffline][2] = "CDS"
			if str(gff_dict[gffline][2]).upper() == "gene":
				gff_dict[gffline][2] = "gene"
			if str(gff_dict[gffline][2]).upper() == "exon":
				gff_dict[gffline][2] = "exon"
			if str(gff_dict[gffline][2]).upper() == "intron":
				gff_dict[gffline][2] = "intron"
			if str(gff_dict[gffline][2]).upper() == "TRNA":
				gff_dict[gffline][2] = "tRNA"
			if str(gff_dict[gffline][2]).upper() == "MRNA":
				gff_dict[gffline][2] = "mRNA"
			if str(gff_dict[gffline][2]).upper() == "RRNA":
				gff_dict[gffline][2] = "rRNA"
		except:
			pass
##Checks for number of times a name appears in annotation list, checks strand direction, creates list (unjoin_features) of those names that code in both directions (and shouldn't be joined)
	join_dict = {}
	unjoin_features = []
	for key1 in feature_dict.keys():
		try:
			name = "Name=%s" %(feature_dict[key1]["Name"])
			splitname = feature_dict[key1]["Name"]
			join_dict[splitname] = []
		except:
			name = "null"
			join_dict[name] = []
		for gffline in gff_dict:
			try:
				if name in gff_dict[gffline][8]:
					join_dict[splitname].append(gff_dict[gffline][6])
			except:
				pass
	for key2 in join_dict.keys():
		if len(set(join_dict[key2])) == 2:
			unjoin_features.append(key2)

##Create .tbl file from entries in feature_dict and gff_dict
	for feature in unique_features:
		#print feature
		firstline = 1
##First, use unjoin_features list to loop through dicts and join only those regions that are oriented in same direction
		if feature in unjoin_features:
			firstline = 1
			for line in feature_dict:
				try:
					if feature == feature_dict[line]["Name"] and gff_dict[line][6] == "+":
						line_pass = line
						if firstline == 1:
							outfile.write("%s\t%s\t%s\n" %(gff_dict[line][3], gff_dict[line][4], gff_dict[line][2]))
							firstline = 0
						elif firstline == 0:
							outfile.write("%s\t%s\n" %(gff_dict[line][3], gff_dict[line][4]))
				except:
					pass
			try:
				outfile.write("\t\t\tgene\t%s\n" %(feature_dict[line_pass]["gene"]))
			except:
				pass
			try:
				outfile.write("\t\t\tproduct\t%s\n" %(feature_dict[line_pass]["product"]))
			except:
				pass
			try:
				outfile.write("\t\t\tnote\t%s\n" %(feature_dict[line_pass]["note"]))
			except:
				pass
			try:
				outfile.write("\t\t\texception\t%s\n" %(feature_dict[line_pass]["exception"]))
				if feature_dict[line_pass]["exception"] == "RNA editing":
					pepoutfile.write(">gnl|ZimmerNMNH|PWS_%05d [gene=%s][protein=%s]\n\n" %(protein_id_num, feature_dict[line_pass]["gene"], feature_dict[line_pass]["product"]))
			except:
				pass
			try:
				outfile.write("\t\t\tcodon_start\t%s\n" %(feature_dict[line_pass]["codon_start"]))
			except:
				pass
			try:
				if "pseudo" in feature_dict[line_pass]:
					outfile.write("\t\t\tpseudo\n")
				else:
					pass
			except:
				pass
				if "CDS" in feature:
					outfile.write("\t\t\tprotein_id\tgnl|ZimmerNMNH|PWS_%05d\n" %(protein_id_num))
					protein_id_num += 1
#			try:
#				outfile.write("\t\t\tprotein_id\t%s\n" %(feature_dict[line_pass]["protein_id"]))
#			except:
#				pass
#			try:
#				outfile.write("\t\t\tdb_xref\t%s\n" %(feature_dict[line_pass]["db_xref"]))
#			except:
#				pass
			firstline = 1
			for line in feature_dict:
				try:
					if feature == feature_dict[line]["Name"] and gff_dict[line][6] == "-":
						line_pass = line
						#print line
						if firstline == 1:
							outfile.write("%s\t%s\t%s\n" %(gff_dict[line][4], gff_dict[line][3], gff_dict[line][2]))
							firstline = 0
						elif firstline == 0:
							outfile.write("%s\t%s\n" %(gff_dict[line][4], gff_dict[line][3]))
				except:
					pass
			try:
				outfile.write("\t\t\tgene\t%s\n" %(feature_dict[line_pass]["gene"]))
			except:
				pass
			try:
				outfile.write("\t\t\tproduct\t%s\n" %(feature_dict[line_pass]["product"]))
			except:
				pass
			try:
				outfile.write("\t\t\tnote\t%s\n" %(feature_dict[line_pass]["note"]))
			except:
				pass
			try:
				outfile.write("\t\t\texception\t%s\n" %(feature_dict[line_pass]["exception"]))
				if feature_dict[line_pass]["exception"] == "RNA editing":
					pepoutfile.write(">gnl|ZimmerNMNH|PWS_%05d [gene=%s][protein=%s]\n\n" %(protein_id_num, feature_dict[line_pass]["gene"], feature_dict[line_pass]["product"]))
			except:
				pass
			try:
				outfile.write("\t\t\tcodon_start\t%s\n" %(feature_dict[line_pass]["codon_start"]))
			except:
				pass
			try:
				if "pseudo" in feature_dict[line_pass]:
					outfile.write("\t\t\tpseudo\n")
				else:
					pass
			except:
				pass
				if "CDS" in feature:
					outfile.write("\t\t\tprotein_id\tgnl|ZimmerNMNH|PWS_%05d\n" %(protein_id_num))
					protein_id_num += 1
#			try:
#				outfile.write("\t\t\tprotein_id\t%s\n" %(feature_dict[line_pass]["protein_id"]))
#			except:
#				pass
#			try:
#				outfile.write("\t\t\tdb_xref\t%s\n" %(feature_dict[line_pass]["db_xref"]))
#			except:
#				pass
##Second, for all "simple" regions with unique names, write to outfile
		else:
			for line in feature_dict:
				#print line
				#print gff_dict[line][2]
				try:
					if feature == feature_dict[line]["Name"]:
						line_pass = line
						if firstline == 1:
							#print gff_dict[line][6], gff_dict[line][3], gff_dict[line][4], gff_dict[line][2]
							if gff_dict[line][6] == "+":
								outfile.write("%s\t%s\t%s\n" %(gff_dict[line][3], gff_dict[line][4], gff_dict[line][2]))
							if gff_dict[line][6] == "-":
								outfile.write("%s\t%s\t%s\n" %(gff_dict[line][4], gff_dict[line][3], gff_dict[line][2]))
							firstline = 0
						elif firstline == 0:
							#print gff_dict[line][3], gff_dict[line][4]
							if gff_dict[line][6] == "+":
								outfile.write("%s\t%s\n" %(gff_dict[line][3], gff_dict[line][4]))
							if gff_dict[line][6] == "-":
								outfile.write("%s\t%s\n" %(gff_dict[line][4], gff_dict[line][3]))
				except:
					pass
				
			try:
				outfile.write("\t\t\tgene\t%s\n" %(feature_dict[line_pass]["gene"]))
			except:
				pass
			try:
				outfile.write("\t\t\tproduct\t%s\n" %(feature_dict[line_pass]["product"]))
			except:
				pass
			try:
				outfile.write("\t\t\tnote\t%s\n" %(feature_dict[line_pass]["note"]))
			except:
				pass
			try:
				outfile.write("\t\t\texception\t%s\n" %(feature_dict[line_pass]["exception"]))
				if feature_dict[line_pass]["exception"] == "RNA editing":
					pepoutfile.write(">gnl|ZimmerNMNH|PWS_%05d [gene=%s][protein=%s]\n\n" %(protein_id_num, feature_dict[line_pass]["gene"], feature_dict[line_pass]["product"]))
			except:
				pass
			try:
				outfile.write("\t\t\tcodon_start\t%s\n" %(feature_dict[line_pass]["codon_start"]))
			except:
				pass
		#	try:
			if "pseudo" in feature_dict[line_pass]:
				outfile.write("\t\t\tpseudo\n")
			else:
				pass
		#	except:
		#		pass
#			try:
#				outfile.write("\t\t\tprotein_id\t%s\n" %(feature_dict[line_pass]["protein_id"]))
#			except:
#				pass
#			try:
#				outfile.write("\t\t\tdb_xref\t%s\n" %(feature_dict[line_pass]["db_xref"]))
#			except:
#				pass
				if "CDS" in feature:
					outfile.write("\t\t\tprotein_id\tgnl|ZimmerNMNH|PWS_%05d\n" %(protein_id_num))
					protein_id_num += 1
	outfile.close()
	pepoutfile.close()

##Merge all single output files into one
if len(newfilelist) > 1:
	mergedoutput = open("MergedOutput.tbl", "w")
	for newfile in newfilelist:
		newopenfile = open(newfile, "r")
		for newline in newopenfile:
			mergedoutput.write(newline)
		newopenfile.close()
	mergedoutput.close()

