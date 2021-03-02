#!/usr/bin/env python

import sys
from Bio import SeqIO

prefix = sys.argv[1]
evalue_setting = float('1e-' + sys.argv[2])
foreign_ratio_setting = float(sys.argv[3])
#prefix = 'Vischeria_C74_genome_v1_annotation'
#diamond_file = open(sys.argv[1], 'rU')
#diamond_file = open('Monodopsis_C73_genome_v1_annotation_prot_sub_diamond.out', 'rU')
diamond_file = open(prefix + '.blastp', 'r')
eggnog_file = open(prefix + '.emapper.annotations', 'r')

## Parse the diamond output, into seq2skingdom_dict
seq2skingdom_dict = {}
for line in diamond_file:
	line = line.strip('\n')
	qseqid,sseqid,pident,length,mismatch,gapopen,qstart,qend,sstart,send,evalue,bitscore,staxids,sskingdoms,skingdoms,sphylums,sscinames,full_sseq = line.split('\t')

	if sskingdoms != 'N/A' and float(evalue) < evalue_setting:
		try:
			seq2skingdom_dict[qseqid].append(sskingdoms)
		except:
			seq2skingdom_dict[qseqid] = [sskingdoms]

## Parse the eggnog output, into seq2annotation_dict
seq2annotation_dict = {}
for line in eggnog_file:
	line = line.strip('\n')
	seq = line.split('\t')[0]
	seq2annotation_dict[seq] = line

## Find seq with predominately bacterial hits
seq2seq_dict = SeqIO.index(prefix + '.fasta', 'fasta')
out_fasta = open(prefix + '_aa_filtered_t1_HGT_candidate_e' +  str(evalue_setting).replace('1e-','') + 'p' + str(foreign_ratio_setting).replace('0.','') + '.fas', 'w')
out_table = open(prefix + '_aa_filtered_t1_HGT_candidate_e' +  str(evalue_setting).replace('1e-','') + 'p' + str(foreign_ratio_setting).replace('0.','') + '.txt', 'w')
HGTcandidate_list = []
for rec in seq2skingdom_dict:
	bacteria_hit_count = seq2skingdom_dict[rec].count('Bacteria') + seq2skingdom_dict[rec].count('Archaea')
	total_hit_count = len(seq2skingdom_dict[rec])
	foreign_ratio = bacteria_hit_count / float(total_hit_count)
	if foreign_ratio > foreign_ratio_setting:
		try:
			#print(seq2annotation_dict[rec])
			out_table.write(seq2annotation_dict[rec] + '\n')
		except:
			#print(rec)
			out_table.write(rec + '\n')
		HGTcandidate_list.append(rec)
		out_fasta.write('>' + rec + '\n' + str(seq2seq_dict[rec].seq) + '\n')

print("No. of HGT candidates found:", len(HGTcandidate_list))
