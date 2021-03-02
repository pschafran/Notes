#! /usr/bin/env python

# Supply PfamScan output files as command line input

import sys
import subprocess

filelist = sys.argv[1:]
goodOutFile = open("pfam_single_clan.txt", "w")
badOutFile = open("pfam_multi_clan.txt", "w")
for file in filelist:
    baseFileName = ".".join(file.split(".")[:-1])
    with open(file, "r") as infile:
        clanList = []
        clanDict = {}
        seqList = []
        for line in infile:
            if line.startswith("#"):
                pass
            else:
                splitline = line.strip("\n").split()
                try:
                    clanList.append(splitline[14])
                    seqList.append(splitline[0])
                except:
                    pass
                try:
                    clanDict[splitline[14]].append(splitline[0])
                except:
                    try:
                        clanDict.update({splitline[14] : [splitline[0]]})
                    except IndexError:
                        pass
        if len(set(clanList)) == 1:
            goodOutFile.write("%s\n" % baseFileName)
        if len(set(clanList)) > 1:
            badOutFile.write("%s\n" % baseFileName)
            seqNum = float(len(set(seqList)))
            for clan in clanDict:
                perc = 100*(float(len(set(clanDict[clan])))/seqNum)
                outfileName = "%s.%s%0.2fperc" %(baseFileName, clan, perc)
                with open(outfileName, "w") as outfile:
                    for seq in set(clanDict[clan]):
                        outfile.write("%s\n" % seq)
                cmd = ["/home/ps997/scripts/getFromFasta.py", "%s" % baseFileName, "%s" % outfileName]
                result = subprocess.run(cmd, capture_output=True, text = True)
                outFastaName = "%s.%s.fasta" %(baseFileName, clan)
                with open(outFastaName, "w") as outfasta:
                    outfasta.write(result.stdout)
goodOutFile.close()
badOutFile.close()
