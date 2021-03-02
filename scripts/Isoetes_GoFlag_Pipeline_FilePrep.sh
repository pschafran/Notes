#! /bin/bash

# Purpose:
# This script prepares files for running GoFlag_Pipeline_Hybpiper.sh
# Check for even number of sequence files, single probe file
# Make sample directories and move sequence files to directories
# Install some dependencies
# 

baseDir="$(pwd)"

# Check for presence of fastq.gz files
if [ -n "$(ls -A ./*.fastq.gz)" ] ; then
	filelist="$(ls -A ./*.fastq.gz)"
	len_filelist=$(wc -w <<< "$filelist")
	for item in $filelist ; do
		echo $item >> filelist.tmp
	done
else
	echo ""
	echo "No fastq.gz files present!"
	echo ""
	exit 1
fi

#Check for appropriate number of fastq.gz files
if [ $((len_filelist%2)) -eq 0 ] ;
then
	echo ""
	echo "Sequence files present..."
	echo ""
else
	echo ""
	echo "ERROR: ODD NUMBER OF FASTQ FILES"
	echo ""
	exit 1
fi

#Check for probe directory and probe file
if [ -d "./probes" ] ; then
	if [ $(find ./probes -name "*fasta" | wc -l) != 0 ] ; then
		echo ""
		echo "Probe file present..."
		echo ""
	else
		echo ""
		echo "ERROR: No probe file found"
		echo ""
		exit 1
	fi
else
	echo ""
	echo "ERROR: No probe directory present"
	echo ""
	exit 1
fi


### Maybe add later to allow different transcriptomes to be used for HybPiper
#Check for transcriptome directory and files
#if [ -d "./transcriptomes" ] ; then
#	if [ $(find ./transcriptomes -name "*fasta" | wc -l) != 0 ] ; then
#		echo ""
#		echo "Transcriptome present..."
#		echo ""
#		echo "Making HybPiper targetfile..."
#		echo ""
#		cd ./transctriptomes/
#		cat *.fasta > allTranscriptomes.fsa
#		mkblastdb -in allTranscriptomes.fsa -dbtype nucl
#		python $baseDir/scripts/makeTargetFile.py 
#	else
#		echo ""
#		echo "WARNING: No transcriptome found"
#		echo ""
#	fi
#else
#	echo ""
#	echo "WARNING: No transcriptome directory present"
#	echo ""
#fi


#Make sample/program directories, install sub-programs from github
echo "Making sample directories..."
mkdir Phyluce
git clone https://github.com/mossmatters/HybPiper.git
git clone https://github.com/tomas-fer/HybPhyloMaker.git
git clone https://github.com/PatrickKueck/FASconCAT-G.git
git clone https://github.com/smirarab/ASTRAL.git
unzip ./ASTRAL/*.zip
export PATH=$PATH:$workingDir/HybPiper/:$workingDir/HybPhyloMaker/:$workingDir/FASconCAT-G/
baseDir="$(pwd)"

#Use python script to split filenames
echo ""
echo "Making sample directories..."
echo ""
python ./scripts/filenamesplitter.py
sleep 10

#Check that files were copied then delete originals
while read line ; do
	if [ -n "$(ls -A ./$line/*.fastq.gz)" ] ; then
		filelist="$(ls -A ./$line/*.fastq.gz)"
		len_filelist=$(wc -w <<< "$filelist")
		if [ $len_filelist = 2 ] ; then
			echo "$line files in place..."
			rm "$line"*.fastq.gz
			echo ""
			echo "Unzipping fastq.gz files..."
			echo ""
			gunzip "./$line/$line"*.fastq.gz
		else
			echo ""
			echo "ERROR: Incorrect number of files in $line. $len_filelist files found"
			echo ""
		fi
	else
		echo ""
		echo "ERROR: Files not found in $line"
		echo ""
	fi
done < "samplelist.txt"

rm filelist.tmp
