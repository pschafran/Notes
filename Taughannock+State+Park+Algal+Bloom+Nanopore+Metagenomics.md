
# Taughannock State Park Algal Bloom Nanopore Metagenomics

## Draft paper: 

## Outline:
1. Trim reads with PoreChop
2. Assemble trimmed reads with Flye v. 2.7 metagenomic mode
3. Polish assembly with several iterations of racon followed by Medaka
4. Assign taxonomy to contigs using Centrifuge, Kraken2, Kaiju, BLASTn, and BLASTx
5. Assign taxonomy to reads with programs that can handle it (not BLAST)
6. Cross-compare taxonomic annotations between programs

------------------

### 1. Trim reads


```python
cd /home/ps997/cyanobloom
porechop -i taughannock_HAB_all.fastq -o taughannock_HAB_all_porechop.fastq -v 3 > porechop.out
```

### Trimming Results:
```
grep "reads" porechop.out

Loading reads
3,742,983 reads loaded
  425,084 / 3,742,983 reads had adapters trimmed from their start (8,560,371 bp removed)
2,680,062 / 3,742,983 reads had adapters trimmed from their end (39,091,431 bp removed)
Splitting reads containing middle adapters
150 / 3,742,983 reads were split based on middle adapters
```

| filename | total_length | number | mean_length | longest | shortest | N_count | Gaps | N50 | N50n | N70 | N70n | N90 | N90n |
|----------|--------------|--------|-------------|---------|----------|---------|------|-----|------|-----|------|----|------|
| taughannock_HAB_all.fastq | 9395949756 | 3742983 | 2510.28 | 110215 | 1 | 0 | 0 | 3645 | 630162 | 2140 | 1315530 | 1212 | 2468748 |
| taughannock_HAB_all_porechop.fastq | 9348234939 | 3743030 | 2497.50 | 110215 | 1 | 0 | 0 | 3651 | 625673 | 2138 | 1307369 | 1205 | 2458719 |

### 2. Assemble
Note: Genome size estimate comes from previous assembly


```python
conda activate flye_2.7
flye --nano-raw taughannock_HAB_all_porechop.fastq -t 12 -g 50m --meta -o taughannock_HAB_all_porechop_flye_meta_2.7
```

### 3. Polish


```python
# racon
cd /home/ps997/cyanobloom/taughannock_HAB_all_porechop_flye_meta_2.7/
ln -s assembly.fasta assembly.racon-iter0.fasta
for i in {1..5} ; do j=$((i-1)) ; \
minimap2 -t 12 -ax map-ont assembly.racon-iter"$j".fasta ../taughannock_HAB_all_porechop.fastq > assembly.racon-iter"$i".sam ; \
racon -t 12 ../taughannock_HAB_all_porechop.fastq assembly.racon-iter"$i".sam assembly.racon-iter"$j".fasta > assembly.racon-iter"$i".fasta 2> racon-iter"$i".out ; \
done
```


```python
# medaka
conda activate medaka
medaka_consensus -i ../taughannock_HAB_all_porechop.fastq -d assembly.racon-iter5.fasta -t 8
```

### Polishing Results

| filename | total_length | number | mean_length | longest | shortest | N_count | Gaps | N50 | N50n | N70 | N70n | N90 | N90n |
|----------|--------------|--------|-------------|---------|----------|---------|------|-----|------|-----|------|----|------|
| assembly.racon-iter5.medaka-consensus.fasta | 48916416 | 2183 | 22407.89 | 2298548 | 80 | 0 | 0 | 90015 | 91 | 42028 | 252 | 12372 | 646 |


### 4.  Taxonomic annotation


```python
mkdir annotations ; mkdir annotation/blastn ; mkdir annotation/centrifuge ; mkdir annotation/kaiju ; mkdir annotation/kraken2

### Commands
blastn -query assembly.racon-iter5.medaka-consensus.fasta -db ~/blastdbs/nt/nt -out blastn.results -outfmt "6 qseqid staxids bitscore" -num_threads 12

centrifuge -x ~/centrifugeDBs/nt/nt -k -p 12 -f assembly.racon-iter5.medaka-consensus.fasta -S centrifuge.results --report-file centrifuge.report

kaiju 

kraken2 

```


```python
RANK="species"
while read contig ; do grep $contig */split_"$RANK"/*fasta | awk -v contig=$contig 'BEGIN { FS=":" }  { print $2"\t"$1 }' ; done < contig.names > contigs.split_"$RANK".tmp
awk -F" |\t|\\\/" '{print $1"\t"$3"\t"$NF}' contigs.split_"$RANK".tmp > contigs.split_"$RANK".tmp2
# for species level use this awk + sed combo instead of line above
# awk -F"|\t|\\\/" '{print $1"\t"$2"\t"$NF}' contigs.split_"$RANK".tmp > contigs.split_"$RANK".tmp2
# sed -i 's/ \(contig\|scaffold\)_.\+\(    .\+     .\+\)$/\1/' contigs.split_species.tmp2 # will have to manually change tabs to tabs in shell

sed -i 's/>//' contigs.split_"$RANK".tmp2
sed -i "s/assembly_species_//" contigs.split_"$RANK".tmp2 
sed -i 's/.fasta//' contigs.split_"$RANK".tmp2 
mv contigs.split_"$RANK".tmp2 contigs.split_"$RANK"
rm contigs.split_"$RANK".tmp

### python
openfile = open("contigs.split_species","r")
contigDict = {}

for line in openfile:
	contig = line.strip("\n").split("\t")[0]
	program = line.strip("\n").split("\t")[1]
	kingdom = line.strip("\n").split("\t")[2]
	try:
		contigDict[contig].update({program:kingdom})
	except:
		contigDict[contig] = {program:kingdom}

outfile = open("classifier_species.table", "w")
outfile.write("Contig\tBlast_species\tCentrifuge_species\tKaiju_species\tKraken_species\n")

for contig in contigDict:
	try:
		centrifuge = contigDict[contig]["centrifuge"]
	except KeyError:
		centrifuge = "NA"
	try:
		blast = contigDict[contig]["blastn"]
	except KeyError:
		blast = "NA"
	try:
		kaiju = contigDict[contig]["kaiju"]
	except KeyError:
		kaiju = "NA"
	try:
		kraken = contigDict[contig]["kraken2"]
	except KeyError:
		kraken ="NA"
	outfile.write("%s\t%s\t%s\t%s\t%s\n" %(contig, blast, centrifuge, kaiju, kraken))


```

### Annotation Results

| Taxonomic Rank | No. of contigs with same annotation by all 4 programs | % of assembly represented by contigs w/ perfect agreement |
|----------------|:-----------------------------------------------------:|:------------------------:
| classifier_superkingdom.table | 1255 | 80.0 |
| classifier_phylum.table | 1119 | 71.2 |
| classifier_family.table | 466 | 29.3 |
| classifier_genus.table | 342 | 20.2 |
| classifier_species.table | 77 | 9.5 | 

### 5. rRNA extraction


```python
barrnap -outseq rRNA.fa --kingdom bac ../assembly.racon-iter5.medaka-consensus.fasta > rRNA.gff
```

barrnap identified 80 rRNA segments
