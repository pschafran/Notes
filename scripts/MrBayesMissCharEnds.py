#! /usr/bin/python
###Usage: MrBayesMissCharEnds.py Alignment.nex
###Replace leading and tailing dashes in a nexus alignment with question marks

import re, sys

file = sys.argv[1]
infile = open(file, "r")
outfile = open("%s_missingchars.nex" %(file), 'w')
search = re.compile("(\t.*\t)(-*)([ACTGactg].*[ACTGactg])(-*)$")

for line in infile:
	if re.search(search,line) == None:
		outfile.write(line)
	else:
		name = re.search(search,line).group(1)
		head = re.search(search,line).group(2)
		seq = re.search(search,line).group(3)
		tail = re.search(search,line).group(4)
		
		head_questions = ""
		tail_questions = ""
		for i in range(len(head)):
			head_questions = head_questions + "?"
		for i in range(len(tail)):
			tail_questions = tail_questions + "?"
		outfile.write("%s%s%s%s\n" %(name , head_questions , seq , tail_questions))
infile.close()
outfile.close()