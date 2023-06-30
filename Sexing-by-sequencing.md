# Sexing-by-sequencing

This is my approach for determine if an individual represents a different sex than a reference genome from another individal using shotgun sequence DNA. Because shotgun sequencing theoretically samples randomly from the whole genome, the data should have even coverage across the whole genome. In practice, PCR bias, sequencing bias, and repeats cause some unevenness in this distribution.

## Materials
* Trimmed/filtered DNA reads from individual 1
* Reference genome from individual 2

## Steps
0. Remove unscaffolded contigs from genome (optional)
Make a list of sequences you want to use in a text file, one per line.
```
/home/ps997/scripts/getFromFasta.py genome.fasta sequence_list.txt > genome_subset.fasta
```

1. Read mapping - find best alignment of each read to the genome
```
bwa index genome_subset.fasta
bwa mem -t 12 genome_subset.fasta samplename_R1.fq.gz samplename_R2.fq.gz | samtools sort -o samplename_mapped_to_genome.bam
```

2. Calculate read depth at each position in the genome
```
samtools depth -a -o samplename_mapped_to_genome.depth.txt
```

3. Plot the depth. Look for "n50contigs" PDF file.
```
/home/ps997/scripts/plotCovAlongContigs.py samplename_mapped_to_genome.depth.txt
```

## Read Normalization
It might be necessary to adjust the number of input reads to make easier comparisons between samples. This can be done by randomly subsampling the samples with larger amounts of data.
1. Figure out what amount of data to use:

| Sample | # of reads (M) | amount of data (Gbp) |
|--------|----------------|----------------------|
|Leiosporoceros H23 |  133.3 |       18.62       |
|Leiosporoceros JC1 |  52.6  |       7.79        |
|Leiosporoceros JC2 |  64.7  |       9.40        |

So we will downsample H23 and JC2 to 52.6 M reads.

2. Subsampling
This can be done with `seqtk`. See example down the page <a href="https://github.com/Li-Lab-BTI/Li-Lab-BTI/wiki/Protocols_wet#flow-cytometry-lb01-buffer">here</a>.
```
seqtk sample -s100 Leiosporoceros_H23.fq.gz 52600000 > Leiosporoceros_H23.fq.gz
seqtk sample -s100 Leiosporoceros_H23.fq.gz 52600000 > Leiosporoceros_H23.fq.gz
```

3. Use these new read files in `bwa` step above.

4. Generate depth .txt file as above with new BAM file.

5. With a depth .txt file for each sample, calculated the difference between the two at each nucleotide position in the genome. For example, if the two sample files have these values:

Sample 1:
```
Ledus.S1	1	25
Ledus.S1	2	21
Ledus.S1	3	23
Ledus.S1	4	19
Ledus.S1	5	18
```
Sample 2:
```
Ledus.S1	1	29
Ledus.S1	2	22
Ledus.S1	3	23
Ledus.S1	4	25
Ledus.S1	5	21
```
Then the new file with difference between the two would look like (which gets subtracted from which depends on the context of your experiment):
```
Ledus.S1	1	4
Ledus.S1	2	1
Ledus.S1	3	0
Ledus.S1	4	6
Ledus.S1	5	3
```
I have done this by loading each input file into its own python dictionary and then looping through and subtracting each line-by-line, but it is VERY slow (hours).

```python
# Read in first file and store to dictionary
phaeom147655 = {}
with open("Phaeomegaceros_fimbriatus_genome.Phaeom_14765-5.20M_mapped.bam.depth.txt.scaffolds.txt","r") as infile:
	for line in infile:
		splitline = line.strip("\n").split("\t")
		try:
			phaeom147655[splitline[0]].update({ splitline[1] : splitline[2]})
		except:
			try:
				phaeom147655[splitline[0]] = { splitline[1] : splitline[2] }
			except:
				phaeom147655.update({splitline[0] : {} })

# Read in second file and store to dictionary
phaeom147658 = {}
with open("Phaeomegaceros_fimbriatus_genome.Phaeom_14765-8.20M_mapped.bam.depth.txt.scaffolds.txt","r") as infile:
	for line in infile:
		splitline = line.strip("\n").split("\t")
		try:
			phaeom147658[splitline[0]].update({ splitline[1] : splitline[2]})
		except:
			try:
				phaeom147658[splitline[0]] = { splitline[1] : splitline[2] }
			except:
				phaeom147658.update({splitline[0] : {} })

# Loop through all sites in genome, subtract one from the other, write to new file
with open("Phaeomegaceros_fimbriatus_genome.Phaeom_14765-8.20M_mapped.bam.depth.txt.scaffolds.delta.txt","w") as outfile:
	for scaffold in phaeom147658:
		for base in phaeom147658[scaffold]:
			diff = int(phaeom147655[scaffold][base]) - int(phaeom147658[scaffold][base])
			outfile.write("%s\t%s\t%s\t%s\n" %(scaffold, base, phaeom147658[scaffold][base], diff))
```

6. Plot this new depth file with `plotCovAlongContigs.py` as above.
