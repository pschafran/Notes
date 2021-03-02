#! /bin/bash

SPECIES="Anthoceros_punctatus" # No spaces!
CORES=12
AUGUSTUS_CONFIG_PATH="/home/ps997/bin/Augustus/config" # No trailing slash /
AUGUSTUS_SCRIPTS_PATH="/home/ps997/bin/Augustus/scripts"
AUGUSTUS_BIN_PATH="/home/ps997/miniconda3/envs/braker/bin"
GENEMARK_PATH="/home/ps997/bin/gmes_linux_64"
ASSEMBLER="flye miniasm masurca-flye masurca-cabog raven haslr redbean" # spaces separated list 
NANOPORE="/home/ps997/Anthoceros_punctatus/Anthoceros_punctatus_20180502_20180817_20180911_nanopore_guppy_porechop.fastq.gz"
ILLUMINA1="/home/ps997/Anthoceros_punctatus/illumina/150806_I262_FCC7F06ACXX_L7_RSZABPI005851-13_1.fq.gz"
ILLUMINA2="/home/ps997/Anthoceros_punctatus/illumina/150806_I262_FCC7F06ACXX_L7_RSZABPI005851-13_2.fq.gz"
RNA1="/home/ps997/Anthoceros_punctatus/rna/Anthoceros_punctatus/replicate1_RNA_R1.fq.gz"
RNA2="/home/ps997/Anthoceros_punctatus/rna/Anthoceros_punctatus/replicate1_RNA_R2.fq.gz"
REFPROT="/home/ps997/hornwort-orthologs/Hornwort_orthogroups.faa"
PLASTOME=""
CHONDROME=""

######### DON'T EDIT BELOW THIS LINE #########

echo "SPECIES: $SPECIES"
echo "AUGUSTUS_CONFIG_PATH: $AUGUSTUS_CONFIG_PATH"
echo "AUGUSTUS_SCRIPTS_PATH: $AUGUSTUS_SCRIPTS_PATH"
echo "AUGUSTUS_BIN_PATH: $AUGUSTUS_BIN_PATH"
echo "GENEMARK_PATH: $GENEMARK_PATH"
echo "RNA files: $RNA1 $RNA2"
echo "Reference protein file: $REFPROT"

# Start script in ~/$SPECIES/assemblies/ directory
# Expects a directory structure like this:
#
# Creates this directory structure:
# ./assemblies/
#	|./flye/
#	|./miniasm/
#	|./masurca-flye/
#	|./masurca-cabog/
#	|./raven/
#	|
#	|./repeat-masking/
#	|	|./flye/
#	|	|	|./0_RemoveOrganelles/
#	|	|	|./1a_RepeatModeler/
#	|	|	|./1b_EDTA/
#	|	|	|./2_ProteinExclude/
#	|	|	|./3_RepeatMasker
#	|	|
#	|	|./masurca-cabog/
#	|	|	|./0_RemoveOrganelles/
#	|	|	|./1a_RepeatModeler/
#	|	|	|./1b_EDTA/
#	|	|	|./2_ProteinExclude/
#	|	|	|./3_RepeatMasker
#	|	|
#	|	|./masurca-flye/
#	|	|	|./0_RemoveOrganelles/
#	|	|	|./1a_RepeatModeler/
#	|	|	|./1b_EDTA/
#	|	|	|./2_ProteinExclude/
#	|	|	|./3_RepeatMasker
#	|	|
#	|	|./miniasm/
#	|	|	|./0_RemoveOrganelles/
#	|	|	|./1a_RepeatModeler/
#	|	|	|./1b_EDTA/
#	|	|	|./2_ProteinExclude/
#	|	|	|./3_RepeatMasker
#	|
#	|./gene-prediction/
#	|	|./flye/
#	|	|	|./1a_RepeatModeler/
#	|	|	|	|./braker-RNA/augustus.hints.[aa|codingseq]  <-- These are probably the final files you want
#	|	|	|	|./braker-RNA+PROT/augustus.hints.[aa|codingseq]  <-- These are probably the final files you want
#	|	|	|./1b_EDTA/
#	|	|	|	|./braker-RNA/augustus.hints.[aa|codingseq]  <-- These are probably the final files you want
#	|	|	|	|./braker-RNA+PROT/augustus.hints.[aa|codingseq]  <-- These are probably the final files you want
#	|	|	|./1c_EDTA_noFilter/
#	|	|	|	|./braker-RNA/augustus.hints.[aa|codingseq]  <-- These are probably the final files you want
#	|	|	|	|./braker-RNA+PROT/augustus.hints.[aa|codingseq]  <-- These are probably the final files you want
#	|	|
#	|	|./masurca-cabog/
#	|	|	|./1a_RepeatModeler/
#	|	|	|	|./braker-RNA/augustus.hints.[aa|codingseq]  <-- These are probably the final files you want
#	|	|	|	|./braker-RNA+PROT/augustus.hints.[aa|codingseq]  <-- These are probably the final files you want
#	|	|	|./1b_EDTA/
#	|	|	|	|./braker-RNA/augustus.hints.[aa|codingseq]  <-- These are probably the final files you want
#	|	|	|	|./braker-RNA+PROT/augustus.hints.[aa|codingseq]  <-- These are probably the final files you want
#	|	|	|./1c_EDTA_noFilter/
#	|	|	|	|./braker-RNA/augustus.hints.[aa|codingseq]  <-- These are probably the final files you want
#	|	|	|	|./braker-RNA+PROT/augustus.hints.[aa|codingseq]  <-- These are probably the final files you want
#	|	|
#	|	|./masurca-flye/
#	|	|	|./1a_RepeatModeler/
#	|	|	|	|./braker-RNA/augustus.hints.[aa|codingseq]  <-- These are probably the final files you want
#	|	|	|	|./braker-RNA+PROT/augustus.hints.[aa|codingseq]  <-- These are probably the final files you want
#	|	|	|./1b_EDTA/
#	|	|	|	|./braker-RNA/augustus.hints.[aa|codingseq]  <-- These are probably the final files you want
#	|	|	|	|./braker-RNA+PROT/augustus.hints.[aa|codingseq]  <-- These are probably the final files you want
#	|	|	|./1c_EDTA_noFilter/
#	|	|	|	|./braker-RNA/augustus.hints.[aa|codingseq]  <-- These are probably the final files you want
#	|	|	|	|./braker-RNA+PROT/augustus.hints.[aa|codingseq]  <-- These are probably the final files you want
#	|	|
#	|	|./miniasm/
#	|	|	|./1a_RepeatModeler/
#	|	|	|	|./braker-RNA/augustus.hints.[aa|codingseq]  <-- These are probably the final files you want
#	|	|	|	|./braker-RNA+PROT/augustus.hints.[aa|codingseq]  <-- These are probably the final files you want
#	|	|	|./1b_EDTA/
#	|	|	|	|./braker-RNA/augustus.hints.[aa|codingseq]  <-- These are probably the final files you want
#	|	|	|	|./braker-RNA+PROT/augustus.hints.[aa|codingseq]  <-- These are probably the final files you want
#	|	|	|./1c_EDTA_noFilter/
#	|	|	|	|./braker-RNA/augustus.hints.[aa|codingseq]  <-- These are probably the final files you want
#	|	|	|	|./braker-RNA+PROT/augustus.hints.[aa|codingseq]  <-- These are probably the final files you want
#		|./$SPECIES_combined_proteins.faa
#		|./$SPECIES_clustered_proteins.faa <-- All proteins from all assemblies and repeat-maskings clustered with CD-HIT at 0.9 similarity
#



RM_CORES=$((CORES/3))

if [ ! -e $NANOPORE ]
	then echo "ERROR: Can't access Nanopore file"
	exit 1
fi
if [! -e $ILLUMINA1 ]
	then echo "ERROR: Can't access Illumina R1 file"
	exit 1
fi
if [! -e $ILLUMINA2 ]
	then echo "ERROR: Can't access Illumina R2 file"
	exit 1
fi
if [ ! -e $RNA1 ]
	then echo "ERROR: Can't access RNA1 file"
	exit 1
fi
if [ ! -e $RNA2 ]
	then echo "ERROR: Can't access RNA2 file"
	exit 1
fi
if [ ! -e $REFPROT ]
	then echo "WARNING: REFPROT file not provided"
fi
if [ ! -e $PLASTOME]
	then echo "WARNING: PLASTOME file not provided"
fi
if [ ! -e $CHONDROME ]
	then echo "WARNING: CHONDROME file not provided"
fi
if [ ! command -v flye &> /dev/null ]
	then echo "ERROR: Flye couldn't be executed"
	exit 1
fi
if [ ! command -v minimap2 &> /dev/null ]
	then echo "ERROR: Minimap2 couldn't be executed"
	exit 1
fi
if [ ! command -v miniasm &> /dev/null ]
	then echo "ERROR: Miniasm couldn't be executed"
	exit 1
fi
if [ ! command -v minipolish &> /dev/null ]
	then echo "ERROR: Minipolish couldn't be executed"
	exit 1
fi
if [ ! command -v raven &> /dev/null ]
	then echo "ERROR: Raven couldn't be executed"
	exit 1
fi
if [ ! command -v masurca &> /dev/null ]
	then echo "ERROR: Masurca couldn't be executed"
	exit 1
fi
if [ ! command -v RepeatMasker &> /dev/null ]
	then echo "ERROR: RepeatMasker couldn't be executed"
	exit 1
fi


for i in $ASSEMBLER
	do if [ ! -d $i ]
		then mkdir $i
	fi
done

if [ -d flye ]
	then cd flye
	flye --nano-raw $NANOPORE -t $CORES -o flye_output 1> flye.out 2> flye.err
	cd ..
fi
if [ -d miniasm ]
	then cd miniasm
	minimap2 -x ava-ont -t $CORES $NANOPORE $NANOPORE | gzip -1 > reads.paf.gz 2> minimap.err
	miniasm -f $NANOPORE reads.paf.gz > "$SPECIES".gfa 2> miniasm.err
	minipolish -t $CORES $NANOPORE "$SPECIES".gfa > "$SPECIES"_polished.gfa 2> minipolish.err
	awk '/^S/{print ">"$2"\n"$3}' "$SPECIES"_polished.gfa | fold > "$SPECIES"_polished.fa
	cd ..
fi
if [ -d raven ]
	then cd raven
	raven -t $CORES $NANOPORE > raven.fa 2> raven.err
	cd ..
fi

if [ -d masurca-flye ]
	then cd masurca-flye
	
fi


HOMEDIR=$(pwd)
if [ ! -d repeat-masking ]
	then mkdir repeat-masking
fi
cd repeat-masking

for i in $ASSEMBLER
	do if [ ! -d $i ]
		then mkdir $i
	fi
	if [ ! -d "$i"/0_RemoveOrganelles ]
		then mkdir "$i"/0_RemoveOrganelles
	fi
	if [ ! -d "$i"/1a_RepeatModeler ]
		then mkdir "$i"/1a_RepeatModeler
	fi
	if [ ! -d "$i"/1b_EDTA ]
		then mkdir "$i"/1b_EDTA
	fi
	if [ ! -d "$i"/2_ProteinExclude ]
		then mkdir "$i"/2_ProteinExclude
	fi
	if [ ! -d "$i"/3_RepeatMasker ]
		then mkdir "$i"/3_RepeatMasker
	fi

	cd "$i"/0_RemoveOrganelles
	CURRENT_PATH=$(pwd)
	echo $CURRENT_PATH
	if [ ! -e 0_assembly.fa ]
		then ln -sf ../../../"$i"/pilon-iter5/pilon-iter5.fasta 0_assembly.fa
	fi
	sed -i --follow-symlinks 's/_pilon_pilon_pilon_pilon_pilon//' 0_assembly.fa
	if [ ! -f organelles.blast ] || [ $(stat --printf="%s"  organelles.blast) -eq 0 ]
		then blastn -query 0_assembly.fa -db ../../organelle_genomes.fa -evalue 0.0001 -outfmt "6 qseqid sseqid pident length qcovs bitscore qstart qend sstart send" -out organelles.blast -num_threads "$CORES"
	fi
	awk -F"\t" '{if ($5 > 90) {print $1} }' organelles.blast | sort | uniq > organelle_contigs.txt
	awk -F"\t" '{if ($5 < 90 && $4 > 1500) {print $0} }' organelles.blast | sort -k1,1 -k7,7 | uniq > potential_misassemblies.txt
	/home/ps997/scripts/removeScaffoldsFromFasta.py 0_assembly.fa organelle_contigs.txt

	cd ../1a_RepeatModeler/
	CURRENT_PATH=$(pwd)
	echo $CURRENT_PATH
	if [ ! -e 1_assembly.fa ]
		then ln -sf ../0_RemoveOrganelles/0_assembly.fa_organelle_contigs.txt-removed.fasta 1_assembly.fa
	fi
	if [ ! -f 1_assembly.fa_db.nsq ]
		then ~/bin/RepeatModeler-2.0.1/BuildDatabase -name 1_assembly.fa_db 1_assembly.fa
	fi
	if [ ! -f 1_assembly.fa_db-families.fa ]
		then ~/bin/RepeatModeler-2.0.1/RepeatModeler -database 1_assembly.fa_db -pa $RM_CORES -LTRStruct >& run.out
	fi
	extract_unknownLTR.py 1_assembly.fa_db-families.fa
	if [ ! -f 1_assembly.fa_db-families.fa.LTRlib.unknown.fa.blastx ] || [ $(stat --printf="%s" 1_assembly.fa_db-families.fa.LTRlib.unknown.fa.blastx) -eq 0 ]
		then blastx -query 1_assembly.fa_db-families.fa.LTRlib.unknown.fa -db /home/fay-wei/bin/Tpases020812 -evalue 1e-10 -num_descriptions 10 -out 1_assembly.fa_db-families.fa.LTRlib.unknown.fa.blastx -num_threads "$CORES"
	fi
	perl /home/fay-wei/bin/Custom-Repeat-Library/transposon_blast_parse.pl --blastx *blastx --modelerunknown 1_assembly.fa_db-families.fa.LTRlib.unknown.fa
	cat 1_assembly.fa_db-families.fa.LTRlib.known.fa identified_elements.txt > 1_assembly.fa_db-families.fa.LTRlib.known.final.fa
	mv unknown_elements.txt 1_assembly.fa_db-families.fa.LTRlib.unknown.final.fa
	/home/ps997/scripts/rename_seq.py 1_assembly.fa_db-families.fa.LTRlib.known.final.fa
	cd $HOMEDIR/repeat-masking
done

conda activate EDTA
export PATH="/home/ps997/bin/NINJA-0.95-cluster_only/NINJA:/home/ps997/bin/cd-hit-v4.8.1-2019-0228:/home/ps997/bin/bbmap:/home/ps997/bin/BBMap:/home/ps997/scripts:/home/fay-wei/bin/racon/build/bin:/home/fay-wei/bin/jellyfish-2.2.5/bin:/home/ps997/miniconda3/envs/EDTA/bin:/home/ps997/miniconda3/condabin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin:/usr/games:/usr/local/games:/usr/lib/jvm/java-8-oracle/bin:/usr/lib/jvm/java-8-oracle/db/bin:/usr/lib/jvm/java-8-oracle/jre/bin:/opt/dell/srvadmin/bin:/home/ps997/edirect:/home/ps997/edirect"
for i in $ASSEMBLER
	do cd "$i"/1b_EDTA/
	CURRENT_PATH=$(pwd)
	echo $CURRENT_PATH
	if [ ! -e 1_assembly.fa ]
		then ln -sf ../0_RemoveOrganelles/0_assembly.fa_organelle_contigs.txt-removed.fasta 1_assembly.fa
	fi
	if [ ! -f 1_assembly.fa.mod.EDTA.TElib.fa ]
		then EDTA.pl --sensitive 1 --anno 1 --evaluate 1 -t $CORES --genome 1_assembly.fa
	fi
	extract_unknownLTR.py 1_assembly.fa.mod.EDTA.TElib.fa
	if [ ! -f 1_assembly.fa.mod.EDTA.TElib.fa.LTRlib.unknown.fa.blastx ] || [ $(stat --printf="%s" 1_assembly.fa.mod.EDTA.TElib.fa.LTRlib.unknown.fa.blastx) -eq 0 ]
		then blastx -query 1_assembly.fa.mod.EDTA.TElib.fa.LTRlib.unknown.fa -db /home/fay-wei/bin/Tpases020812 -evalue 1e-10 -num_descriptions 10 -out 1_assembly.fa.mod.EDTA.TElib.fa.LTRlib.unknown.fa.blastx -num_threads "$CORES"
	fi
	perl /home/fay-wei/bin/Custom-Repeat-Library/transposon_blast_parse.pl --blastx *blastx --modelerunknown 1_assembly.fa.mod.EDTA.TElib.fa.LTRlib.unknown.fa
	cat 1_assembly.fa.mod.EDTA.TElib.fa.LTRlib.known.fa identified_elements.txt > 1_assembly.fa.mod.EDTA.TElib.fa.LTRlib.known.final.fa
	mv unknown_elements.txt 1_assembly.fa.mod.EDTA.TElib.fa.LTRlib.unknown.final.fa
	/home/ps997/scripts/rename_seq.py 1_assembly.fa.mod.EDTA.TElib.fa.LTRlib.known.final.fa
	cd $HOMEDIR/repeat-masking
done
conda activate base

for i in $ASSEMBLER
	do cd "$i"/2_ProteinExclude/ 
	CURRENT_PATH=$(pwd)
	echo $CURRENT_PATH
	if [ ! -e 2_RepeatModeler.LTRlib.known.fa ]
		then ln -sf ../1a_RepeatModeler/1_assembly.fa_db-families.fa.LTRlib.known.final.fa.renamed.fa 2_RepeatModeler.LTRlib.known.fa
	fi
	if [ ! -e 2_EDTA.LTRLib.known.fa ]
		then ln -sf ../1b_EDTA/1_assembly.fa.mod.EDTA.TElib.fa.LTRlib.known.final.fa.renamed.fa 2_EDTA.LTRLib.known.fa
	fi
	if [ ! -f RepeatModeler2uniprot_plant_blast.out ] || [ $(stat --printf="%s" RepeatModeler2uniprot_plant_blast.out) -eq 0 ]
		then blastx -db /home/fay-wei/bin/uniprot/uniprot_sprot_plants.fasta -query 2_RepeatModeler.LTRlib.known.fa -out RepeatModeler2uniprot_plant_blast.out -num_threads $CORES
		perl /home/ps997/scripts/ProtExcluder1.1/ProtExcluder.pl RepeatModeler2uniprot_plant_blast.out 2_RepeatModeler.LTRlib.known.fa
	fi
	if [ ! -f EDTA2uniprot_plant_blast.out ] || [ $(stat --printf="%s" EDTA2uniprot_plant_blast.out) -eq 0 ]
		then blastx -db /home/fay-wei/bin/uniprot/uniprot_sprot_plants.fasta -query 2_EDTA.LTRLib.known.fa -out EDTA2uniprot_plant_blast.out -num_threads $CORES
		perl /home/ps997/scripts/ProtExcluder1.1/ProtExcluder.pl EDTA2uniprot_plant_blast.out 2_EDTA.LTRLib.known.fa
	fi

	cd $HOMEDIR/repeat-masking
done

for i in $ASSEMBLER
	do cd "$i"/3_RepeatMasker
	CURRENT_PATH=$(pwd)
	echo $CURRENT_PATH
	if [ ! -d EDTA_out ]
		then mkdir EDTA_out
	fi
	if [ ! -d RepeatModeler_out ]
		then mkdir RepeatModeler_out
	fi
	if [ ! -e 3_RepeatModeler_LTRLib ]
		then ln -sf ../2_ProteinExclude/2_RepeatModeler.LTRlib.known.fanoProtFinal 3_RepeatModeler_LTRLib
	fi
	if [ ! -e 3_EDTA_LTRLib ]
		then ln -sf ../2_ProteinExclude/2_EDTA.LTRLib.known.fanoProtFinal 3_EDTA_LTRLib
	fi
	if [ ! -e 3_assembly.fa ]
		then ln -sf ../1a_RepeatModeler/1_assembly.fa 3_assembly.fa
	fi

	if [ ! -e RepeatModeler_out/3_assembly.fa.masked ]
		then /home/ps997/bin/RepeatMasker/RepeatMasker -noisy -a -gff -u -pa $CORES -lib 3_RepeatModeler_LTRLib 3_assembly.fa
		rsync --progress --remove-source-files 3_assembly.fa.* RepeatModeler_out
	fi
	if [ ! -e EDTA_out/3_assembly.fa.masked ]
		then /home/ps997/bin/RepeatMasker/RepeatMasker -noisy -a -gff -u -pa $CORES -lib 3_EDTA_LTRLib 3_assembly.fa
		rsync --progress --remove-source-files 3_assembly.fa.* EDTA_out
	fi
	
	if [ ! -e EDTA_detailed_summary ]
		then /home/ps997/bin/RepeatMasker/util/buildSummary.pl -species "$SPECIES" -useAbsoluteGenomeSize EDTA_out/3_assembly.fa.out > EDTA_detailed_summary
	fi
	if [ ! -e RepeatModeler_detailed_summary ]
		then /home/ps997/bin/RepeatMasker/util/buildSummary.pl -species "$SPECIES" -useAbsoluteGenomeSize RepeatModeler_out/3_assembly.fa.out > RepeatModeler_detailed_summary
	fi
	cd $HOMEDIR/repeat-masking
done
cd $HOMEDIR

if [ ! -d gene-prediction ]
	then mkdir gene-prediction
fi
cd gene-prediction

for i in $ASSEMBLER
	do if [ ! -d $i ]
		then mkdir $i
	fi
done

for i in $ASSEMBLER
	do cd "$i"
	if [ ! -d 1a_RepeatModeler ]
		then mkdir 1a_RepeatModeler
	fi
	if [ ! -d 1b_EDTA ]
		then mkdir 1b_EDTA
	fi
	if [ ! -d 1c_EDTA_nofilter ]
		then mkdir 1c_EDTA_nofilter
	fi
	for j in 1a_RepeatModeler 1b_EDTA 1c_EDTA_nofilter
		do cd "$j"
		if [ ! -e RNA_R1.fq.gz ]
			then ln -sf $RNA1 RNA_R1.fq.gz
		fi
		if [ ! -e RNA_R2.fq.gz ]
			then ln -sf $RNA2 RNA_R2.fq.gz
		fi
		if [ ! -e Hornwort_orthogroups.faa ]
			then ln -sf $REFPROT Hornwort_orthogroups.faa
		fi
		cd ..
	done
	
	cd 1a_RepeatModeler/
	CURRENT_PATH=$(pwd)
	echo $CURRENT_PATH
	if [ ! -e 1_RepeatModeler_masked_assembly.fa ]
		then ln -sf ../../../repeat-masking/"$i"/3_RepeatMasker/RepeatModeler_out/3_assembly.fa.masked 1_RepeatModeler_masked_assembly.fa
	fi
	if [ ! -e RNA_mapped.bam ] || [ $(stat --printf="%s" RNA_mapped.bam) -lt 100 ]
		then hisat2-build -p "$CORES" 1_RepeatModeler_masked_assembly.fa 1_RepeatModeler_masked_assembly >& hisat-build.out
		hisat2 -p "$CORES" -x 1_RepeatModeler_masked_assembly -1 RNA_R1.fq.gz -2 RNA_R2.fq.gz 2> hisat-align.out | samtools view -b | samtools sort > RNA_mapped.bam
		samtools flagstat RNA_mapped.bam > RNA_mapped.bam.flagstats
		samtools stats RNA_mapped.bam > RNA_mapped.bam.samstats
	fi
	if [ ! -e RNA_mapped.bam.alignmentsummary ]
		then picard CollectAlignmentSummaryMetrics REFERENCE_SEQUENCE=1_RepeatModeler_masked_assembly.fa INPUT=RNA_mapped.bam OUTPUT=RNA_mapped.bam.alignmentsummary
	fi
	
	cd ../1b_EDTA/
	CURRENT_PATH=$(pwd)
	echo $CURRENT_PATH
	if [ ! -e 1_EDTA_masked_assembly.fa ]
		then ln -sf ../../../repeat-masking/"$i"/3_RepeatMasker/EDTA_out/3_assembly.fa.masked 1_EDTA_masked_assembly.fa
	fi
	if [ ! -e RNA_mapped.bam ] || [ $(stat --printf="%s" RNA_mapped.bam) -lt 100 ]
		then hisat2-build -p $CORES 1_EDTA_masked_assembly.fa 1_EDTA_masked_assembly >& hisat-build.out
		hisat2 -p "$CORES" -x 1_EDTA_masked_assembly -1 RNA_R1.fq.gz -2 RNA_R2.fq.gz 2> hisat-align.out | samtools view -b | samtools sort > RNA_mapped.bam
		samtools flagstat RNA_mapped.bam > RNA_mapped.bam.flagstats
		samtools stats RNA_mapped.bam > RNA_mapped.bam.samstats
	fi
	if [ ! -e RNA_mapped.bam.alignmentsummary ]
		then picard CollectAlignmentSummaryMetrics REFERENCE_SEQUENCE=1_EDTA_masked_assembly.fa INPUT=RNA_mapped.bam OUTPUT=RNA_mapped.bam.alignmentsummary
	fi
	
	cd ../1c_EDTA_nofilter
	CURRENT_PATH=$(pwd)
	echo $CURRENT_PATH
	
	if [ ! -e 1_assembly.fa.mod.MAKER.masked.fa ]
		then ln -sf ../../../repeat-masking/"$i"/1b_EDTA/1_assembly.fa.mod.MAKER.masked 1_assembly.fa.mod.MAKER.masked.fa
	fi
	if [ ! -e RNA_mapped.bam ] || [ $(stat --printf="%s" RNA_mapped.bam) -lt 100 ]
	then hisat2-build -p $CORES 1_assembly.fa.mod.MAKER.masked.fa 1_assembly.fa.mod.MAKER.masked >& hisat-build.out
		hisat2 -p "$CORES" -x 1_assembly.fa.mod.MAKER.masked -1 RNA_R1.fq.gz -2 RNA_R2.fq.gz 2> hisat-align.out | samtools view -b | samtools sort > RNA_mapped.bam
		samtools flagstat RNA_mapped.bam > RNA_mapped.bam.flagstats
		samtools stats RNA_mapped.bam > RNA_mapped.bam.samstats
	fi
	if [ ! -e RNA_mapped.bam.alignmentsummary ]
		then picard CollectAlignmentSummaryMetrics REFERENCE_SEQUENCE=1_assembly.fa.mod.MAKER.masked.fa INPUT=RNA_mapped.bam OUTPUT=RNA_mapped.bam.alignmentsummary
	fi
	
	cd ../1a_RepeatModeler/
	CURRENT_PATH=$(pwd)
	echo $CURRENT_PATH
	conda activate braker
	if [ ! -f braker-RNA/augustus.hints.aa ]
		then k=1
		for k in {1..10000} ; do if [ -d /home/ps997/bin/Augustus/config/species/Sp_"$k" ] ; then true ; else break ; fi ; done
		echo "Running Braker as Sp_$k"
		/home/ps997/miniconda3/envs/braker/bin/braker.pl --genome 1_RepeatModeler_masked_assembly.fa --species Sp_"$k" --bam RNA_mapped.bam --verbosity 3 --cores $CORES --GENEMARK_PATH=/home/ps997/bin/gmes_linux_64 --AUGUSTUS_CONFIG_PATH=/home/ps997/bin/Augustus/config --AUGUSTUS_BIN_PATH=/home/ps997/miniconda3/envs/braker/bin/ --AUGUSTUS_SCRIPTS_PATH=/home/ps997/bin/Augustus/scripts 1>braker.out 2>braker.err
		rsync -ar --remove-source-files braker/* braker-RNA/
	fi
	if [ ! -f braker-RNA+PROT/augustus.hints.aa ]
		then k=1
		for k in {1..10000} ; do if [ -d /home/ps997/bin/Augustus/config/species/Sp_"$k" ] ; then true ; else break ; fi ; done
		echo "Running Braker as Sp_$k"
		/home/ps997/miniconda3/envs/braker/bin/braker.pl --genome 1_RepeatModeler_masked_assembly.fa --species Sp_"$k" --bam RNA_mapped.bam --prot_seq Hornwort_orthogroups.faa --prg=gth --gth2traingenes --verbosity 3 --cores $CORES --GENEMARK_PATH=/home/ps997/bin/gmes_linux_64 --AUGUSTUS_CONFIG_PATH=/home/ps997/bin/Augustus/config --AUGUSTUS_BIN_PATH=/home/ps997/miniconda3/envs/braker/bin/ --AUGUSTUS_SCRIPTS_PATH=/home/ps997/bin/Augustus/scripts 1>braker.log 2>braker.err
		rsync -ar --remove-source-files braker/* braker-RNA+PROT/
	fi
	
	cd ../1b_EDTA/
	CURRENT_PATH=$(pwd)
	echo $CURRENT_PATH
	if [ ! -f braker-RNA/augustus.hints.aa ]
		then k=1
		for k in {1..10000} ; do if [ -d /home/ps997/bin/Augustus/config/species/Sp_"$k" ] ; then true ; else break ; fi ; done
		echo "Running Braker as Sp_$k"
		/home/ps997/miniconda3/envs/braker/bin/braker.pl --genome 1_EDTA_masked_assembly.fa --species Sp_"$k" --bam RNA_mapped.bam --verbosity 3 --cores $CORES --GENEMARK_PATH=/home/ps997/bin/gmes_linux_64 --AUGUSTUS_CONFIG_PATH=/home/ps997/bin/Augustus/config --AUGUSTUS_BIN_PATH=/home/ps997/miniconda3/envs/braker/bin/ --AUGUSTUS_SCRIPTS_PATH=/home/ps997/bin/Augustus/scripts 1>braker.out 2>braker.err
		rsync -ar --remove-source-files braker/* braker-RNA/
	fi
	if [ ! -f braker-RNA+PROT/augustus.hints.aa ]
		then k=1
		for k in {1..10000} ; do if [ -d /home/ps997/bin/Augustus/config/species/Sp_"$k" ] ; then true ; else break ; fi ; done
		echo "Running Braker as Sp_$k"
		/home/ps997/miniconda3/envs/braker/bin/braker.pl --genome 1_EDTA_masked_assembly.fa --species Sp_"$k" --bam RNA_mapped.bam --prot_seq Hornwort_orthogroups.faa --prg=gth --gth2traingenes --verbosity 3 --cores $CORES --GENEMARK_PATH=/home/ps997/bin/gmes_linux_64 --AUGUSTUS_CONFIG_PATH=/home/ps997/bin/Augustus/config --AUGUSTUS_BIN_PATH=/home/ps997/miniconda3/envs/braker/bin/ --AUGUSTUS_SCRIPTS_PATH=/home/ps997/bin/Augustus/scripts 1>braker.out 2>braker.err
		rsync -ar --remove-source-files braker/* braker-RNA+PROT/
	fi
	
	cd ../1c_EDTA_nofilter/
	CURRENT_PATH=$(pwd)
	echo $CURRENT_PATH
	if [ ! -f braker-RNA/augustus.hints.aa ]
		then k=1
		for k in {1..10000} ; do if [ -d /home/ps997/bin/Augustus/config/species/Sp_"$k" ] ; then true ; else break ; fi ; done
		echo "Running Braker as Sp_$k"
		/home/ps997/miniconda3/envs/braker/bin/braker.pl --genome 1_assembly.fa.mod.MAKER.masked.fa --species Sp_"$k" --bam RNA_mapped.bam --verbosity 3 --cores $CORES --GENEMARK_PATH=/home/ps997/bin/gmes_linux_64 --AUGUSTUS_CONFIG_PATH=/home/ps997/bin/Augustus/config --AUGUSTUS_BIN_PATH=/home/ps997/miniconda3/envs/braker/bin/ --AUGUSTUS_SCRIPTS_PATH=/home/ps997/bin/Augustus/scripts 1>braker.out 2>braker.err
		rsync -ar --remove-source-files braker/* braker-RNA/
	fi
	if [ ! -f braker-RNA+PROT/augustus.hints.aa ]
		then k=1
		for k in {1..10000} ; do if [ -d /home/ps997/bin/Augustus/config/species/Sp_"$k" ] ; then true ; else break ; fi ; done
		echo "Running Braker as Sp_$k"
		/home/ps997/miniconda3/envs/braker/bin/braker.pl --genome 1_assembly.fa.mod.MAKER.masked.fa --species Sp_"$k" --bam RNA_mapped.bam --prot_seq Hornwort_orthogroups.faa --prg=gth --gth2traingenes --verbosity 3 --cores $CORES --GENEMARK_PATH=/home/ps997/bin/gmes_linux_64 --AUGUSTUS_CONFIG_PATH=/home/ps997/bin/Augustus/config --AUGUSTUS_BIN_PATH=/home/ps997/miniconda3/envs/braker/bin/ --AUGUSTUS_SCRIPTS_PATH=/home/ps997/bin/Augustus/scripts 1>braker.out 2>braker.err
		rsync -ar --remove-source-files braker/* braker-RNA+PROT/
	fi
	conda activate base
	cd $HOMEDIR/gene-prediction
done
conda activate base


