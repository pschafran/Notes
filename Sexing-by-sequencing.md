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
