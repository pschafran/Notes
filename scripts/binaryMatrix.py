#! /usr/bin/python

import sys, random

cutoff = int(sys.argv[1])
listlen = int(sys.argv[2])
matrix = []
binmatrix = []
count = 1
outfile = open("binmatrix.csv", 'w')

while int(count) <= int(listlen):
	randnum = random.randint(1,listlen)
	matrix.append(randnum)
	count += 1
bincount = 1
for item in matrix:
	if int(item) < int(cutoff):
		if bincount < listlen:
			binmatrix.append(1)
			outfile.write("1,")
			bincount += 1
		elif bincount == listlen:
			binmatrix.append(1)
			outfile.write("1")
	else:
		if bincount < listlen:
			binmatrix.append(0)
			outfile.write("0,")
			bincount += 1
		elif bincount == listlen:
			binmatrix.append(0)
			outfile.write("0")
outfile.write("\n")
#print matrix
#print binmatrix



