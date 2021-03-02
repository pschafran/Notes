#! /usr/bin/python
'''
USAGE: makeTargetFile.py blastDatabase list_of_reference_files.txt transcriptomes_used_for_BLASTdb.fasta

To create target file for HybPiper:
	BLASTn bait files from GoFlag against Isoetes transcriptomes
	Extract BLAST hits from transcriptomes
	Rename scaffolds with locus ID
	Remove any duplicates
'''
import sys, subprocess

blast_database = sys.argv[1]
reference_file_list = sys.argv[2]
transcriptome_file = sys.argv[3]


def getReferences():
	openReference_file_list = open(reference_file_list, "r")
	IsoetesRef = open("Isoetes_GoFlag_references.fasta", "w")
	LycophytesRef = open("Lycophytes_GoFlag_references.fasta", "w")
	AllRef = open("All_GoFlag_references.fasta", "w")
	for file in openReference_file_list:
		openRef = open(file.strip("\n"), "r")
		for line in openRef:
			AllRef.write(line)
			if "Isoetes" in line:
				writeOutIsoetes = 1
				IsoetesRef.write(line)
			if "Lycophytes" in line:
				writeOutLycos = 1
				LycophytesRef.write(line)
			if ">" in line and "Isoetes" not in line:
				writeOutIsoetes = 0
			if ">" in line and "Lycophytes" not in line:
				writeOutLycos = 0
			if ">" not in line:
				if writeOutLycos == 1:
					LycophytesRef.write(line)
				if writeOutIsoetes == 1:
					IsoetesRef.write(line)
	IsoetesRef.close()
	LycophytesRef.close()
	AllRef.close()
	openReference_file_list.close()
	

def getIsoetesScaffolds():
	output_file = "Isoetes_BLAST_to_transcriptomes.out"
	subprocess.call(['blastn', '-query','Isoetes_GoFlag_references.fasta', '-db', blast_database, '-out', output_file, '-outfmt', '6'])
	scaffoldsOut = open("Isoetes_GoFlag_target_scaffolds.fasta", 'w')
	transcriptomesIn = open(transcriptome_file, "r")
	scaffoldNamesIn = open(output_file, "r")
	for scaffoldLine in scaffoldNamesIn:
		splitScaffold = scaffoldLine.split("\t")
		scaffoldName = splitScaffold[1]
		refName = splitScaffold[0].split("_")
		locusName = refName[0]
		writeOut = 0
		for line in transcriptomesIn:
			if ">" in line and scaffoldName.strip(">\n") == line.strip(">\n"):
				print "Found scaffold %s" %(line)
				writeOut = 1
				scaffoldsOut.write("%s-%s\n" %(line.strip("\n"), locusName))
			if ">" in line and scaffoldName not in line:
				writeOut = 0
			if ">" not in line:
				if writeOut == 1:
					scaffoldsOut.write(line)
		transcriptomesIn.seek(0)
	scaffoldNamesIn.close()
	scaffoldsOut.close()
	transcriptomesIn.close()

def getLycophyteScaffolds():
	output_file = "Lycophytes_BLAST_to_transcriptomes.out"
	subprocess.call(['blastn', '-query','Lycophytes_GoFlag_references.fasta', '-db', blast_database, '-out', output_file, '-outfmt', '6'])
	scaffoldsOut = open("Lycophyte_GoFlag_target_scaffolds.fasta", 'w')
	transcriptomesIn = open(transcriptome_file, "r")
	scaffoldNamesIn = open(output_file, "r")
	for scaffoldLine in scaffoldNamesIn:
		splitScaffold = scaffoldLine.split("\t")
		scaffoldName = splitScaffold[1]
		refName = splitScaffold[0].split("_")
		locusName = refName[0]
		writeOut = 0
		for line in transcriptomesIn:
			if ">" in line and scaffoldName.strip(">\n") == line.strip(">\n"):
				print "Found scaffold %s" %(line)
				writeOut = 1
				scaffoldsOut.write("%s-%s\n" %(line.strip("\n"), locusName))
			if ">" in line and scaffoldName not in line:
				writeOut = 0
			if ">" not in line:
				if writeOut == 1:
					scaffoldsOut.write(line)
		transcriptomesIn.seek(0)
	scaffoldNamesIn.close()
	scaffoldsOut.close()
	transcriptomesIn.close()
	
def getAllScaffolds():
	output_file = "AllRef_BLAST_to_transcriptomes.out"
	subprocess.call(['blastn', '-query','All_GoFlag_references.fasta', '-db', blast_database, '-out', output_file, '-outfmt', '6'])
	scaffoldsOut = open("All_GoFlag_target_scaffolds.fasta", 'w')
	transcriptomesIn = open(transcriptome_file, "r")
	scaffoldNamesIn = open(output_file, "r")
	for scaffoldLine in scaffoldNamesIn:
		splitScaffold = scaffoldLine.split("\t")
		scaffoldName = splitScaffold[1]
		refName = splitScaffold[0].split("_")
		locusName = refName[0]
		writeOut = 0
		for line in transcriptomesIn:
			if ">" in line and scaffoldName.strip(">\n") == line.strip(">\n"):
				print "Found scaffold %s" %(line)
				writeOut = 1
				scaffoldsOut.write("%s-%s\n" %(line.strip("\n"), locusName))
			if ">" in line and scaffoldName not in line:
				writeOut = 0
			if ">" not in line:
				if writeOut == 1:
					scaffoldsOut.write(line)
		transcriptomesIn.seek(0)
	scaffoldNamesIn.close()
	scaffoldsOut.close()
	transcriptomesIn.close()

def getUniqueIsoetesScaffolds():
	inScaffolds = open("Isoetes_GoFlag_target_scaffolds.fasta", "r")
	outScaffolds = open("Isoetes_Unique_GoFlag_target_scaffolds.fasta", "w")
	inScaffoldDict = {}
	for line in inScaffolds:
		if ">" in line:
			try:
				joinLine = "".join(inScaffoldDict[scaffoldName]) # This must be before next line to merge previously parsed seq
				inScaffoldDict[scaffoldName] = joinLine
			except:
				pass
			scaffoldName = line
			inScaffoldDict[scaffoldName] = []
		if ">" not in line:
			stripLine = line.strip("\n")
			inScaffoldDict[scaffoldName].append(stripLine)
	joinLine = "".join(inScaffoldDict[scaffoldName]) # This line must be after for loop to join final sequence
	inScaffoldDict[scaffoldName] = joinLine
	for key in sorted(inScaffoldDict.keys()):
		outScaffolds.write(key)
		outScaffolds.write("%s\n" %(inScaffoldDict[key]))
	inScaffolds.close()
	outScaffolds.close()

def getUniqueLycophyteScaffolds():
	inScaffolds = open("Lycophyte_GoFlag_target_scaffolds.fasta", "r")
	outScaffolds = open("Lycophyte_Unique_GoFlag_target_scaffolds.fasta", "w")
	inScaffoldDict = {}
	for line in inScaffolds:
		if ">" in line:
			try:
				joinLine = "".join(inScaffoldDict[scaffoldName]) # This must be before next line to merge previously parsed seq
				inScaffoldDict[scaffoldName] = joinLine
			except:
				pass
			scaffoldName = line
			inScaffoldDict[scaffoldName] = []
		if ">" not in line:
			stripLine = line.strip("\n")
			inScaffoldDict[scaffoldName].append(stripLine)
	joinLine = "".join(inScaffoldDict[scaffoldName]) # This line must be after for loop to join final sequence
	inScaffoldDict[scaffoldName] = joinLine
	for key in sorted(inScaffoldDict.keys()):
		outScaffolds.write(key)
		outScaffolds.write("%s\n" %(inScaffoldDict[key]))
	inScaffolds.close()
	outScaffolds.close()

def getUniqueAllScaffolds():
	inScaffolds = open("All_GoFlag_target_scaffolds.fasta", "r")
	outScaffolds = open("All_Unique_GoFlag_target_scaffolds.fasta", "w")
	inScaffoldDict = {}
	for line in inScaffolds:
		if ">" in line:
			try:
				joinLine = "".join(inScaffoldDict[scaffoldName]) # This must be before next line to merge previously parsed seq
				inScaffoldDict[scaffoldName] = joinLine
			except:
				pass
			scaffoldName = line
			inScaffoldDict[scaffoldName] = []
		if ">" not in line:
			stripLine = line.strip("\n")
			inScaffoldDict[scaffoldName].append(stripLine)
	joinLine = "".join(inScaffoldDict[scaffoldName]) # This line must be after for loop to join final sequence
	inScaffoldDict[scaffoldName] = joinLine
	for key in sorted(inScaffoldDict.keys()):
		outScaffolds.write(key)
		outScaffolds.write("%s\n" %(inScaffoldDict[key]))
	inScaffolds.close()
	outScaffolds.close()

getReferences()
#getIsoetesScaffolds()
#getLycophyteScaffolds()
getAllScaffolds()
#getUniqueIsoetesScaffolds()
#getUniqueLycophyteScaffolds()
getUniqueAllScaffolds()







































