#! /usr/bin/python
#This reports # scaffolds, total length, and average length of scaffolds.fasta produced by SPAdes
#Usage: ./spadesStats.py path(s)/to/files.fasta
#
import sys

def calc_N50(seq_length_list):
	seq_length_list.sort()
	total_num_base = sum(seq_length_list)
	base_sums = 0
	current_length = 0
	while base_sums < total_num_base/2.0:
		current_length = seq_length_list.pop()
		base_sums += current_length
	return current_length

print "Scaffolds\tFileName\tAvgLen\tTotalLen\tN50"
for file in sys.argv[1:]:
	infile=open(file, 'r')
	scaffoldcount=0
	totallength=0
	avglength=0
	lngths_lst=[]
	for line in infile:
		if ">" in line:
			splitline=line.split('_')
			length=int(splitline[3])
			lngths_lst.append(length)
			totallength=totallength+length
			#print line
			scaffoldcount+=1
	avglength=totallength/scaffoldcount
	print "%d\t%s\t%d\t%d\t%d" %(scaffoldcount,file,avglength,totallength, calc_N50(lngths_lst))
	
	
	