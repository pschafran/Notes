#!/usr/bin/env python

import sys
from read_write_files import read_fasta_lists

def calc_N50(seq_length_list):
	seq_length_list.sort()
	total_num_base = sum(seq_length_list)
	base_sums = 0
	current_length = 0
	while base_sums < total_num_base/2.0:
		current_length = seq_length_list.pop()
		base_sums += current_length
	return current_length

names, seqs = read_fasta_lists(sys.argv[1])

lngths_lst = [len(seq) for seq in seqs]

print "The total number of sequences is %d" % (len(seqs))
print "The average sequence length is %d" % (sum(lngths_lst)/float(len(seqs)))
print "The total number of bases is %d" % (sum(lngths_lst))
print "The minimum sequence length is %d" % (min(lngths_lst))
print "The maximum sequence length is %d" % (max(lngths_lst))
print "The N50 is %d" % (calc_N50(lngths_lst))
lngths_lst = [len(seq) for seq in seqs]
under150=0
over500=0
over1000=0
over2000=0

print "Median Length = %d" % (lngths_lst[len(lngths_lst)/2])

for item in lngths_lst:
	if item < 150:
		under150+=1
	if item >= 500:
		over500+=1
		if item >= 1000:
			over1000+=1
			if item >= 2000:
				over2000+=1

print "contigs < 150bp = " + str(under150)
print "contigs >= 500bp = " + str(over500)
print "contigs >= 1000bp = " + str(over1000)
print "contigs >= 2000bp = " + str(over2000)


# total_cov=0
# for name in names:
# 	cov = float(name.split()[-1].replace(",", ""))
# 	total_cov+=cov
# print "The average coverage is %.2f" % (total_cov/len(names))	
