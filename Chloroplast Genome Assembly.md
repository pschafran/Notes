# Chloroplast Genome Assembly

## What is the chloroplast genome aka plastome?
* Genome separate from the plant's nuclear genome,
* Mostly contains genes related to photosynthesis and protein synthesis (but some crucial genes have relocated to the nuclear genome)
* Derived from ancient bacteria, circular structure
* Mostly have a 4 part structure in this order: Large single copy subunit (LSC), inverted repeat A (IRA), small single copy subunit (SSC), inverted repeat B (IRB)
* Inverted repeats have identical sequence, but reversed directions
* Most plants have roughly 50:50 ratio of two isoforms of the plastome, with SSC in both directions in relation to other units (LSC-IRA-SSC-IRB and LSC-IRA-rSSC-IRB)
* Structure and gene content mostly conserved among land plants, occasional inversions, expansions of inverted repeat regions

## How to sequence a chloroplast genome
* Genome skimming takes advantage of higher abundance of organelles compared to nuclear genome
* Randomly (shotgun) sequencing 5-10 Gbp of genomic DNA usually yields enough chloroplast data for a good assembly

## Assembly
Requirements:
* Unix-like operating system (macOS or Linux). On Windows, can run Linux in a virtual machine (see [here](https://www.lifewire.com/install-ubuntu-linux-windows-10-steps-2202108) and [here](https://itsfoss.com/install-linux-in-virtualbox/) for instructions)
* Anaconda/Miniconda package manager [link](https://docs.conda.io/en/latest/miniconda.html)
* Geneious [link](https://www.geneious.com/)

### Step 1.
First always check your data to make sure no major issues with sequencing. [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) is a common tool for this. See link for download, instructions, and examples of interpretation. Some warnings are ok, but should not get any errors.

### Step 2.
Use [fastp](https://github.com/OpenGene/fastp) to trim low quality bases from the ends of reads and poly-G tails if sequenced with Illumina 2-color chemistry. Remove any reads that are too short.
```
# To install fastp
conda install -c bioconda fastp

# To run
fastp --in1 R1.fastq --in2 R2.fastq --out1 fastp_output_R1.fastq --out2 fastp_output_R2.fastq --cut_front --cut_tail --trim_poly_g --detect_adapter_for_pe --length_required 50 --thread 8
```
After fastp, look at report in fastp.html. Can also evaluate output fastq file with FastQC.

### Two alternative assembly methods
#### NOVOPlasty
NOVOPlasty assemblies organelle genomes directly from total genomic data. Usually it's pretty good if you have sufficient coverage.
1. Install
```
git clone https://github.com/ndierckx/NOVOPlasty.git
```
2. Select a seed sequence. rbcL gene from closest available species is a good choice. Save sequence in FASTA file.
3. Set up NOVOPlasty config file. Change project name, read files for each sample. See  example file for explanation of settings
```
Project:
#-----------------------
Project name          = Project_Name
Type                  = chloro
Genome Range          = 120000-200000
K-mer                 = 33
Max memory            =
Extended log          = 0
Save assembled reads  = no
Seed Input            = Seed.fasta
Extend seed directly  = no
Reference sequence    =
Variance detection    =
Chloroplast sequence  =
#
Dataset 1:
#-----------------------
Read Length           = 151
Insert size           = 300
Platform              = illumina
Single/Paired         = PE
Combined reads        =
Forward reads         = fastp_output_R1.fastq
Reverse reads         = fastp_output_R2.fastq
Store Hash            =
#
Heteroplasmy:
#-----------------------
MAF                   =
HP exclude list       =
PCR-free              =
#
Optional:
#-----------------------
Insert size auto      = yes
Use Quality Scores    = no
Output path           =
```
4. Run NOVOPlasty
```
./NOVOPlasty4.3.1.pl config.txt
```
5. Check [output files](https://github.com/ndierckx/NOVOPlasty/wiki/Output-files).
* `Circularized_assembly_...fasta` If lucky, NOVOPlasty will identify one circular structure and output in this file.
* `Option_number_...fasta` More likely you will have multiple of these files. One of these files will hopefully contain a correct assembly.
* `Merged_contigs_ProjectName.txt` contains information about how contigs were put together to create the `Option...fasta` assemblies. You should see some that match expected structure of the plastome e.g. 01-02-03-02

If one good circular assembly is not found, it may be that the multiple isoforms described in intro were detected. The lengths of these assemblies will be the same, just different ordering of the contigs. In this case it is arbitrary which assembly is used. Sometimes minor variants between contigs result in multiple possibly haplotypes, in which case read mapping to each assembly should be used to determine the dominant haplotype.

#### Manual Assembly
NOVOPlasty may fail to identify a good candidate assembly, especially with lower coverage datasets. It may also be valuable to manually assembly the plastome as a check against the NOVOPlasty assembly.
###### Extract only chloroplast reads from total dataset. Requires a relatively close reference plastome (can use NOVOPlasty output if it seems complete).
1. Install tools
```
conda install -c bioconda bowtie2 spades
```
2. Map reads to reference plastome and create new FASTQ files of reads likely originating from the plastome.
```
bowtie2-build plastome.fasta
bowtie2 -p 8 -x plastome.fasta -1 fastp_output_R1.fastq -2 fastp_output_R2.fastq --very-sensitive-local --al-conc chloroplast_reads_R%.fastq
```
3. Assemble chloroplast reads into contigs
```
spades.py -t 8 -1 chloroplast_reads_R1.fastq -2 chloroplast_reads_R2.fastq -o spades --careful
```
4. Check spades assembly. Contig names will contain length and coverage information.
```
grep ">" spades/scaffolds.fasta
>NODE01_len_83209_cov_503.2035
$
$ >N
$ >
$ >
$ >
$ >
$ >
.
.
.
```
Ideally 3 contigs will be assembled (large subunit, inverted repeat, small subunit). Usually there will be many short contigs (<1000 bp) that can be discarded. There should be order of magnitude or greater difference between good chloroplast contigs and spurious ones (and inverted repeat contig should have ~2x coverage of the LSC/SSC).

### Step 4.
If assemblies look reasonable, create preliminary annotations of NOVOPlasty and spades results. Upload files to [GeSeq](https://chlorobox.mpimp-golm.mpg.de/geseq.html). Select linear sequences unless using NOVOPlasty `Circularized_assembly_...fasta`. Add any NCBI RefSeqs that are fairly close (e.g. same family). Check box to include MPI-MP chloroplast (land plant) references. Leave other options on default. After jobs are complete, download GenBank file and open in Geneious. Load read files too.

### Step 5.
Try to find correct orientation of your contigs to stitch into a draft. Align contigs to reference plastome. Depending on similar the reference is to the sample, this may take some trial and error with different Geneious settings and plugins (minimap2, Mauve, LASTZ). It is possible that there have been rearrangements. Can also use annotations to guide piecing contigs together.  





## Troubleshooting
1. Many contigs
