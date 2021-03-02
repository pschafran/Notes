## Find HGT candidate genes
1. Do BLAST search of proteins, use taxonomy of best hits to identify genes with with possible/likely horizontal gene transfer origin.
```
# Do DIAMOND search of C. richardii proteins
diamond blastp -d /home/fay-wei/bin/NCBI_db/nr_20200613.dmnd -q Crichardii.v2.1.allTrs.pep.fa -o Crichardii.v2.1.allTrs.pep.nr.blastp --max-target-seqs 100 --query-cover 25 --subject-cover 25 --evalue 1e-10 --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore staxids sskingdoms skingdoms sphylums sscinames full_sseq -b12 -c1 -p 24
```
2. Measure alien index (AI) and ingroup/outgroup hit proportions of results.
```
/home/ps997/scripts/alienIndex.py --ingroup Eukaryota --missing ingroup --file Crichardii.v2.1.allTrs.pep.nr.blastp --output Crichardii.v2.1.allTrs.pep.nr.blastp.AI
```
3. Get proteins with AI > 0 indicating BLAST evalue equal or greater in the outgroup, and at least 25% of all hits in the outgroup, to avoid spurious nr entries.
```
awk -F"\t" '{ if ($8 >= 0 && $5 > 25) print $0 }' Crichardii.v2.1.allTrs.pep.nr.blastp.AI > Crichardii.v2.1.allTrs.pep.nr.blastp.AI.candidateHGT
cut -f 1 Crichardii.v2.1.allTrs.pep.nr.blastp.AI.candidateHGT > Crichardii.v2.1.allTrs.pep.nr.blastp.AI.candidateHGT.names.txt
# getFromFasta.py Crichardii.v2.1.allTrs.pep.fa Crichardii.v2.1.allTrs.pep.nr.blastp.AI.candidateHGT.names.txt > Crichardii.v2.1.allTrs.pep.nr.blastp.AI.candidateHGT.fasta
```

4. Get orthogroups (from Viridiplantae Orthofinder analysis, not described here) containing HGT candidates.
```
mkdir HGT_orthogroups
for i in ~/hornwort-orthologs/orthofinder_in/OrthoFinder/Results_Jan21_1/Orthogroup_Sequences/*.fa ; do grep -f Crichardii.v2.1.allTrs.pep.nr.cov33.blastp.AI.candidateHGT.names.txt $i && cp $i HGT_orthogroups/ ; done
cd HGT_orthogroups
```

5. Search 1KP database for additional members of each orthogroup
```
for i in OG*.fa ; do diamond blastp -d /home/ps997/1kp/1kp-translated-protein.dmnd -q $i -o "$i".1kp.blasthits --max-target-seqs 100 --query-cover 25 --subject- cover 25 -b12 -c1 --evalue 1e-10 -p 24 --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore full_sseq ; done
```

6. Retrieve BLAST hits for Ceratopteris in each OG from original DIAMOND Results
```
for i in OG*.fa ; do grep "Ceric" $i | sed 's/>//' | while read line ; do grep "$line" ../Crichardii.v2.1.allTrs.pep.nr.blastp ; done > "$i".nr.blasthits ; done
```

7. Combine BLAST results with orthogroups, remove/combine duplicate sequences
```
# Extract accession IDs from blast results
for i in OG*.fa ; do cut -f 2 "$i".nr.blasthits | sort | uniq > "$i".nr.accessions ; done
for i in OG*.fa ; do cut -f 2 "$i".1kp.blasthits | sort | uniq > "$i".1kp.accessions ; done

# Make backup copies of OGs
for i in OG*.fa ; do cp $i "$i".bak ; done

# Extract each unique accession and sequence present in blast results, append to orthogroup fasta file
for i in OG*.fa ; do while read line ; do grep -m 1 $line "$i".nr.blasthits ; done < "$i".nr.accessions | awk -F"\t" '{print ">"$17"-"$2"\n"$NF}' | sed 's/ /_/g' | awk -F";|-" '{ if ( $1 ~ /^>/ ) print $1"-"$NF ; else print $0 }' | cat >> "$i".nr.blasthits.fasta ; done

for i in OG*.fa ; do while read line ; do grep -m 1 $line "$i".1kp.blasthits ; done < "$i".1kp.accessions | awk -F"\t" '{print ">"$2"\n"$NF}' | sed 's/scaffold-//' | awk -F"-" '{ if ( $1 ~ /^>/ ) print ">"$3"-"$1"-"$2 ; else print $0 }' | sed 's/->/-/' | cat >> "$i".1kp.blasthits.fasta ; done

# Combine BLAST fastas to OG fasta
for i in OG*.fa ; do cat $i "$i".nr.blasthits.fasta "$i".1kp.blasthits.fasta > "$i".final.fasta ; done

# Remove duplicate sequences by name, combine entries with identical AA/DNA sequence (renamed with all combined entry names)
for i in OG*.final.fasta ; do removeDuplicateSeqFasta.py $i ; done
```

8. Use Pfam annotation to potentially erroneous groupings by looking for files with multiple clans (grouping of families believed to have a common evolutionary ancestor). If multiple clans found, each is split into a new fasta for downstream analysis. Two output TXT files list all input files that had a single clan, and all those with multiple clans.
```
for i in OG*.dupSeqsCombined.fasta ; do pfam_scan.pl -cpu 24 -fasta $i -dir /home/ps997/bin/PfamScanDB -outfile "$i".pfamscan ; done
/home/ps997/scripts/pfamscanSingleClan.py *.pfamscan
```

9. Continue with only groups/clans containing the candidate HGTs
```
grep -f ../Crichardii.v2.1.allTrs.pep.nr.cov25.blastp.AI.candidateHGT.names.txt *dupSeqsCombined.fasta*fasta | awk -F":" '{print $1}' | sort | uniq | while read line ; do clustalo --threads 12 -i "$line" -o "$line".clustal ; done
```

10. For single-clan orthogroups, align, trim, and build tree
```
# Align each orthogroup, trim sites with < 75% overlap and remove sequences with < 50% "good" sites, then infer phylogeny
while read i ; do clustalo --threads 12 -i $i -o "$i".clustal ; done < pfam_single_clan.txt
for i in *clustal ; do /home/ps997/bin/trimal/source/trimal -in $i -out "$i".trimal.fa -resoverlap 0.75 -seqoverlap 50 ; done
for i in *trimal.fa ; do iqtree -s $i -m MFP -B 5000 -T 8 ; done
```

11. Parse trees to find likely HGT based on eukaryotic sequences sister or nested in bacterial sequences
```
# Run traverse_tree script modified from Fay_Wei's script. Some options hardcoded in lines 14,15,16 need to be modified for each experiment
traverse_tree.py *contree

```


## Alien Index
Comparative measurement of BLAST scores to indicate strength of evidence for horizontal gene transfer, by comparing best ingroup evalue with best outgroup evalue. **Use script with DIAMOND BLASTp results generated with this specific format**:
```
diamond blastp -d ~/bin/NCBI_db/nr_20200613.dmnd -q query -o query_diamond2nr.out --max-target-seqs 100 --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore staxids sskingdoms skingdoms sphylums sscinames

/home/ps997/scripts/alienIndex.py --ingroup Eukaryota --file query_diamond2nr.out --output query_diamond2nr.alienIndex.txt
```
Output includes gene name, best ingroup and outgroup evalues, and alien index. AI ranges from about +/- 460, with scores > 0 showing better outgroup score than ingroup. AI of 10-20 has been used as cutoff for conclusive evidence of ingroup/outgroup gene origin. Output is not sorted by default, can be sorted by AI with:
```
sort -k 4,4nr  query_diamond2nr.alienIndex.txt
```
You can see what ingroups are available for your BLASTp file, try:
```
# Superkingdom
cut -f 14 query_diamond2nr.out | sort | uniq

# Kingdom
cut -f 15 query_diamond2nr.out | sort | uniq

# Phylum
cut -f 16 query_diamond2nr.out | sort | uniq

# Genus
cut -f 17 query_diamond2nr.out | awk -F";| " '{print $1}' | sort | uniq
```
In rare cases, a taxonomic group's name might be used at multiple taxonomic levels. In this case, the --taxon parameter can be used to specify (superkingdom, kingdom, phylum, genus). Run `alienIndex.py --help` for more details.



###### Citations
* Fan, X. et al. 2020. Phytoplankton pangenome reveals extensive prokaryotic horizontal gene transfer of diverse functions. <i>Sci. Adv.</i> 6 (18): eaba0111. [doi: 10.1126/sciadv.aba0111](https://advances.sciencemag.org/content/6/18/eaba0111)
* Gladyshev, E.A., M. Meselson, and A.R. Arkhipova. 2008. Massive horizontal gene transfer in bdelloid rotifers. <i>Science</i> 320 (5880): 1210-1213. [doi: 10.1126/science.1156407](https://science.sciencemag.org/content/320/5880/1210)
