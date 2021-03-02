#! /home/ps997/miniconda3/bin/python
### example usage: python plotKs.py <tsv_file> <output.png>
import sys
import numpy as np
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from scipy import stats

inFile = open(sys.argv[1], 'r')
KsValues = []
for line in inFile:
	if line.startswith("\t"):
		pass
	else:
		spLine = line.strip().split("\t")
		if spLine[1] and spLine[3]:
			if float(spLine[1]) <= 0.25 or float(spLine[2]) <= 0.25:  # you can change these to select for percent aligned/percent identity
				pass
			else:	
				if spLine[8]:
#					print (spLine[8])  # debug
					if float(spLine[8]) < 10:
						KsValues.append(float(spLine[8]))
		else:
			pass
#print(KsValues)
bins = np.arange(-10, 10, 0.0075) # fixed bin size at 0.075
plt.xlim([0, 10])  # might have to adjust x lim and y lim based on your data or you can just comment these out
plt.ylim([0, 2])
#plt.hist(KsValues, bins=bins, alpha=1, color = 'cornflowerblue')
#plt.title("I. taiwanensis whole paranome Ks") # change title etc...
plt.xlabel('Ks')
plt.ylabel('Frequency')
sns.distplot(KsValues, bins = bins)
plt.savefig(sys.argv[2], dpi=300)
# plt.show()  # uncomment this if you want to see the plot- doesn't work on my python installation for some reason
