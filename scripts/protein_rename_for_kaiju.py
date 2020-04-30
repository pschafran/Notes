#! /usr/bin/python

import re

taxFile = open("~/bin/blobtools/data/nodesDB.txt","r")
proteinFile = open("complete.protein.all.faa","r")
proteinOut = open("complete.protein.all.renamed.faa","w")
taxDict = {}

for taxLine in taxFile:
	splitline = taxLine.strip("\n").split("\t")
	taxDict[splitline[2]] = splitline[0]
counter=1
errorCounter = 0
genusCounter = 0
for protLine in proteinFile:
	if ">" in protLine:
		try:
			if "[[" in protLine:
				splitline=re.findall("\[\[(.*?)\](.*?)\]", protLine)
				sciName = "".join(splitline[0])
				try:
					taxID=taxDict[sciName]
					proteinOut.write(">%s_%s\n" %(counter,taxID))
				except KeyError:
					print("WARNING: %s not found in taxID list. Trying genus only..." %(sciName))
					try:
						genus=sciName.split("\s")[0]
						genusTaxID=taxDict[genus]
						proteinOut.write(">%s_%s\n" %(counter,genusTaxID))
					except:
						print("ERROR: %s not found in taxID list. TaxID will be assigned to root. New seq id is >%s_1" %(sciName, protLine.strip("\n")))
						proteinOut.write(">%s_1\n" %(counter))
						genusCounter += 1
			else:
				splitline=re.split("\[|\]", protLine)
				sciName2=splitline[-2]
				try:
					taxID2=taxDict[sciName2]
					proteinOut.write(">%s_%s\n" %(counter,taxID2))
				except KeyError:
					print("WARNING: %s not found in taxID list. Trying genus only..." (sciName))
					try:
						genus2=sciName2.split("\s")[0]
						genusTaxID2=taxDict[genus2]
						proteinOut.write(">%s_%s\n" %(counter,genusTaxID2))
					except:
						print("ERROR: %s not found in taxID list. TaxID will be assigned to root. New seq id is >%s_1" %(sciName2, protLine.strip("\n")))
						proteinOut.write(">%s_1\n" %(counter))
						genusCounter += 1
		except IndexError:
			print("ERROR: %s not formatted correctly. TaxID will be assigned to root. New seq id is >%s_1" %(protLine.strip("\n"), counter))
			proteinOut.write(">%s_1\n" %(counter))
			errorCounter += 1
		counter+=1
	else:
		proteinOut.write(protLine)
taxFile.close()
proteinFile.close()
proteinOut.close()

print("Renaming finished!")
print("%s sequences processed" %(counter-1))
print("%s species could not be identified, were assigned to genus taxID" %(genusCounter))
print("%s erroneous sequences assigned to root" %(errorCounter))
errorPct = (errorCounter/counter)*100
print("%f2.2 \% error rate" %(errorPct))