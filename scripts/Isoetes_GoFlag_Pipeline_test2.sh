#! /bin/bash
baseDir="$(pwd)"
cd ./HybPiper/

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