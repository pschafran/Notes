#! /bin/bash



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
python ./scripts/filenamesplitter.py
sleep 5

#Check that files were copied then delete originals
while read line ; do
	if [ -n "$(ls -A ./$line/*.fastq.gz)" ] ; then
		filelist="$(ls -A ./$line/*.fastq.gz)"
		len_filelist=$(wc -w <<< "$filelist")
		if [ $len_filelist = 2 ] ; then
			echo "$line files in place..."
			#rm "$line"*.fastq.gz
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



### Run HybPiper and Cleanup ###
# Change into HybPiper directory
cp ./samplelist.txt ./HybPiper/
cd ./HybPiper

# Create shell scripts to run reads_first.py for each sample -- works better than trying to parallelize in same shell
while read i ; do
	cwd="$(pwd)"
	seqfiles="$(ls -A ../$i/*.fastq)"
	seqfiles=$(echo $seqfiles|tr '\n' ' ')
	echo -e "cd $cwd\npython ./reads_first.py -r $seqfiles -b *.fasta --prefix $i --bwa && python ./cleanup.py $i" > HybPiperReadsFirst_$i.sh ;
	chmod +x HybPiperReadsFirst_$i.sh
	open -a Terminal HybPiperReadsFirst_$i.sh
done < "samplelist.txt"


# Hold until all HybPiper scripts complete
hybseqProcs="$(ps aux | grep -c 'HybPiper')"

until [ $hybseqProcs > 1 ] ; do 
	sleep 15
	hybseqProcs="$(ps aux | grep -c 'HybPiper')"
done

echo ""
echo "All samples processed by HybPiper..."
echo ""
echo "Retrieving sequence statistics..."
echo ""

cd ./HybPiper

# Get sequence lengths and statistics from HybPiper results
python ./get_seq_lengths.py *.fasta samplelist.txt dna > ./sequence_lengths.txt
python ./hybpiper_stats.py ./sequence_lengths.txt samplelist.txt > ./statistics.txt

# Identify potential introns/non-coding sequences for each sample
while read i ; do
	python ./intronerate.py --prefix $i
done < "samplelist.txt"

mkdir ./dna
mkdir ./introns
mkdir ./supercontigs
mkdir ./paralogs

# Extract probe, intron, supercontig, and paralog sequences from individual samples to multi-fasta
python ./retrieve_sequences.py ./*.fasta . dna
mv *.FNA ./dna
python ./retrieve_sequences.py ./*.fasta . intron
mv *introns.fasta ./introns
python ./retrieve_sequences.py ./*.fasta . supercontig
mv *supercontig.fasta ./supercontigs

echo "" > paralog_loci.txt
for i in ./*/genes_with_paralog_warnings.txt ; do
	cat $i >> paralog_loci.txt
	echo "" >> paralog_loci.txt
done

while read i; 
	do python ./paralog_investigator.py ./$i; 
done < "samplelist.txt"
while read i; 
	do python ./paralog_retriever.py samplelist.txt ./$i > ./$i.paralogs.fasta; 
done < "paralog_loci.txt"
mv *.paralogs.fasta ./paralogs/

# Align all loci with MAFFT and trim alignments with trimAL
baseDir="$(pwd)"
cd ./HybPiper/dna/
for i in *.FNA ; do
	fileSize="$(stat -c%s $i)"
	if [ $fileSize = 0 ] ; then
		echo "$i empty..."
		rm $i
	else
		mafft --auto $i > MAFFT.dna.$i.fasta
		trimal -in MAFFT.dna.$i.fasta -out AUTO.TRIM.dna.$i.fasta -automated1
		trimal -in MAFFT.dna.$i.fasta -out STRICT.TRIM.dna.$i.fasta -strict
		trimal -in MAFFT.dna.$i.fasta -out SUPERSTRICT.TRIM.dna.$i.fasta -gt 0.8 -st 0.01 -resoverlap 0.9 -seqoverlap 80
	fi
done ;
mkdir mafft && mv MAFFT.dna* ./mafft/
mkdir auto && mv AUTO.TRIM.dna* ./auto/
mkdir strict && mv STRICT.TRIM.dna* ./strict/
mkdir superstrict && mv SUPERSTRICT.TRIM.dna* ./superstrict/


cd ../introns/
for i in *introns.fasta ; do
	fileSize="$(stat -c%s $i)"
	if [ $fileSize = 0 ] ; then
		echo "$i empty..."
		rm $i
	else
		mafft --auto $i > MAFFT.introns.$i
		trimal -in MAFFT.introns.$i -out AUTO.TRIM.introns.$i -automated1
		trimal -in MAFFT.introns.$i -out STRICT.TRIM.introns.$i -strict
		trimal -in MAFFT.introns.$i -out SUPERSTRICT.TRIM.introns.$i -gt 0.8 -st 0.01 -resoverlap 0.9 -seqoverlap 80
	fi
done ;
mkdir mafft && mv MAFFT.introns* ./mafft/
mkdir auto && mv AUTO.TRIM.introns* ./auto/
mkdir strict && mv STRICT.TRIM.introns* ./strict/
mkdir superstrict && mv SUPERSTRICT.TRIM.introns* ./superstrict/

cd ../supercontigs/
for i in *supercontig.fasta ; do
	fileSize="$(stat -c%s $i)"
	if [ $fileSize = 0 ] ; then
		echo "$i empty..."
		rm $i
	else
		mafft --auto $i > MAFFT.supercontig.$i
		trimal -in MAFFT.supercontig.$i -out AUTO.TRIM.supercontig.$i -automated1
		trimal -in MAFFT.supercontig.$i -out STRICT.TRIM.supercontig.$i -strict
		trimal -in MAFFT.supercontig.$i -out SUPERSTRICT.TRIM.supercontig.$i -gt 0.8 -st 0.01 -resoverlap 0.9 -seqoverlap 80
	fi
done ;
mkdir mafft && mv MAFFT.supercontig* ./mafft/
mkdir auto && mv AUTO.TRIM.supercontig* ./auto/
mkdir strict && mv STRICT.TRIM.supercontig* ./strict/
mkdir superstrict && mv SUPERSTRICT.TRIM.supercontig* ./superstrict/

cd ../paralogs/
for i in *paralogs.fasta ; do
	fileSize="$(stat -c%s $i)"
	if [ $fileSize = 0 ] ; then
		echo "$i empty..."
		rm $i
	else
		mafft --auto $i > MAFFT.paralogs.$i
		trimal -in MAFFT.paralogs.$i -out AUTO.TRIM.paralogs.$i -automated1
		trimal -in MAFFT.paralogs.$i -out STRICT.TRIM.paralogs.$i -strict
		trimal -in MAFFT.paralogs.$i -out SUPERSTRICT.TRIM.paralogs.$i -gt 0.8 -st 0.01 -resoverlap 0.9 -seqoverlap 80
	fi
done ;
mkdir mafft && mv MAFFT.paralogs* ./mafft/
mkdir auto && mv AUTO.TRIM.paralogs* ./auto/
mkdir strict && mv STRICT.TRIM.paralogs* ./strict/
mkdir superstrict && mv SUPERSTRICT.TRIM.paralogs* ./superstrict/

cd ../

#Set variables for 100%, 75%, and 50% of taxa present
totalSamples="$(grep -c "" samplelist.txt)"
fiftyPctSamples=$(awk -v totalSamples="${totalSamples}" -v multiplier=5 'BEGIN{print (totalSamples*multiplier/10)}')
seventyfivePctSamples=$(awk -v totalSamples="${totalSamples}" -v multiplier=7 'BEGIN{print (totalSamples*multiplier/10)}')

# Nested loops to parse through all markers and all trimming parameters. 
# Make trees with IQ-TREE and RAxML, then sort based on how many taxa are present 
for i in dna introns supercontigs paralogs ; do
	cd ./$i/
	for j in mafft auto strict superstrict ; do
		cd ./$j/
		mkdir 50pct
		mkdir 75pct
		mkdir 100pct
		
		# Concatenate all loci using FASconCAT-g and make trees with IQ-TREE and RAxML
		perl $baseDir/FASconCAT-G/FASconCAT-G_v1.04.pl -s
		iqtree -nt AUTO -s FcC_supermatrix.fas -m MFP -bb 5000 -alrt 5000 -nm 5000
		raxmlHPC -s FcC_supermatrix.fas -n FcC_supermatrix -m GTRGAMMA -f a -x 1234 -N 100 -d -p 1234
		
		for k in MAFFT*fasta AUTO*fasta STRICT*fasta SUPERSTRICT*fasta ; do
			iqtree -nt AUTO -s $k -m MFP -bb 5000 -alrt 5000 -nm 5000
			raxmlHPC -s $k -n $k -m GTRGAMMA -f a -x 1234 -N 100 -d -p 1234
			num="$(grep -c ">" $k )"
			if [ $num = $totalSamples ] ; then
				cp *$k* ./100pct && cp *$k* ./75pct && cp *$k* ./50pct && rm *$k*
				echo "$k" >> 100pctTaxa.txt && echo "$k" >> 75pctTaxa.txt && echo "$k" >> 50pctTaxa.txt
			elif [ $num -ge $seventyfivePctSamples ] ; then
				cp *$k* ./75pct && cp *$k* ./50pct && rm *$k*
				echo "$k" >> 75pctTaxa.txt && echo "$k" >> 50pctTaxa.txt
			elif [ $num -ge $fiftyPctSamples ] ; then
				cp *$k* ./50pct && rm *$k*
				echo "$k" >> 50pctTaxa.txt
			else
				echo "$k has less than 50% of taxa"
				echo "$k" >> LessThan50pctTaxa.txt
			fi
		done
		for l in 50pct 75pct 100pct ; do
			cd ./$l
			cat RAxML_bestTree* > RAxML.$i.$j.$l.trees
			cat *.treefile > IQTREE.$i.$j.$l.trees
			java -jar $baseDir/ASTRAL/astral.*.jar -i RAxML.$i.$j.$l.trees -o RAxML.ASTRAL.$i.$j.$l.trees
			java -jar $baseDir/ASTRAL/astral.*.jar -i IQTREE.$i.$j.$l.trees -o IQTREE.ASTRAL.$i.$j.$l.trees
			cd ../
		done
		cd ../
	done
	cd ../
done