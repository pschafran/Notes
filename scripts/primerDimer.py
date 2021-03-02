#! /usr/bin/python

import sys


file = sys.argv[1]
openfile = open(file, 'r')
outfile = "PrimerDimers.csv"
openoutfile = open(outfile, 'w')
forward =[]
reverse = []
for line in openfile:
	if "Forward" in line and "Reverse" in line:
		splitline = line.split("_")
		openoutfile.write("%s,%s\n" %(splitline[2],splitline[5]))
		forward.append(splitline[2])
		reverse.append(splitline[5])
openoutfile.close()
openfile.close()
uniquesForward = set(forward)
uniquesReverse = set(reverse)
print "Forward Uniques"
print sorted(uniquesForward)
print "Reverse Uniques"
print sorted(uniquesReverse)
