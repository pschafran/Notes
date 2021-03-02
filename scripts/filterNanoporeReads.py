count = 0
outfile = open("Phaeomegaceros_split_reads.fastq","w")
for record in dict:
	print(record)
	mean = np.mean(dict[record].letter_annotations["phred_quality"])
	std = np.std(dict[record].letter_annotations["phred_quality"])
	limit = mean - std
	baseBegin = 1
	window = 100
	baseEnd = window
	outlierWindows = []
	if mean > 20:
		for base in dict[record].letter_annotations["phred_quality"]:
			windowMean = np.mean(dict[record].letter_annotations["phred_quality"][baseBegin:baseEnd])
			if windowMean < 10 and baseEnd < len(dict[record].seq)-window:
				outlierWindows.extend(range(baseBegin,baseEnd))
			baseBegin += 1
			baseEnd += 1
	if len(outlierWindows) > 1:
		segmentDict = {}
		segmentNum = 1
		count += 1
		outlierList = sorted(set(outlierWindows))
		startBase = 1
		stopBase = 1
		for base in range(1,len(dict[record].seq)):
			if base not in outlierList:
				if int(base-1) in outlierList:
					startBase = base
				if int(base+1) in outlierList:
					stopBase = base
					segment = "segment%s" % segmentNum
					segmentDict[segment] = [startBase,stopBase]
					SeqIO.write(dict[record][startBase:stopBase],outfile,"fastq")
					segmentNum += 1
		print(outlierList)
		print(segmentDict)
	else:
		SeqIO.write(dict[record],outfile,"fastq")
outfile.close()
#		print(baseCount)