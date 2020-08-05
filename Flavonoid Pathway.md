# Flavonoid Synthesis 

PAL (PHENYLALANINE AMMONIA LYASE) --> C4H (CINNAMATE 4-HYDROXYLASE) --> 4CL (4-COUMARATE-COA LIGASE) --> CHS (CHALCONE SYNTHASE) --> CHI (CHALCONE ISOMERASE)

## PAL
Selected Marchantia polymorpha PAL sequences based on "Transcriptome Analysis of Marchantin Biosynthesis from the Liverwort Marchantia polymorpha" Takahashi and Asakawa 2017. Natural Product Communications https://journals.sagepub.com/doi/pdf/10.1177/1934578X1701200831  
```
blastp -query Hornwort_proteins_all.faa -db ~/blastdbs/flavonoids/Marchantia_PAL -evalue 1e-10 -outfmt 6 -out Hornwort_proteins_all.PAL.blastp 
```

## C4H
Use CinnamateHydroxylases.fasta provided by Kevin  
```
blastp -query Hornwort_proteins_all.faa -db ~/blastdbs/flavonoids/CinnamateHydroxylases -evalue 1e-10 -outfmt 6 -out Hornwort_proteins_all.C4H.blastp
```

## 4CL
Use 4CL.fasta provided by Kevin  
```
blastp -query Hornwort_proteins_all.faa -db ~/blastdbs/flavonoids/4CL -evalue 1e-10 -outfmt 6 -out Hornwort_proteins_all.4CL.blastp
```

## CHS
Download chalcone synthase orthologs from eggNOG ENOG410IKG2, HMM from eggNOG, and chalcone synthase-like from Marchantia polymorpha (tr|Q68BL6|Q68BL6_MARPA Chalcone synthase-like OS=Marchantia paleacea subsp. diptera OX=93925 GN=MpCHSLK1 PE=2 SV=1)  
```
blastp -query Hornwort_proteins_all.faa -db ~/blastdbs/flavonoids/Chalcone_synthase_ENOG410IKG2 -evalue 1e-10 -outfmt 6 -out Hornwort_proteins_all.CHS.blastp 

blastp -query Hornwort_proteins_all.faa -db ~/blastdbs/flavonoids/Chalcone_synthase-like_Marchantia -evalue 1e-10 -outfmt 6 -out Hornwort_proteins_all.CHS-like.blastp 

hmmsearch ~/blastdbs/flavonoids/CHS_eggnog.hmm.txt Hornwort_proteins_all.faa > Hornwort_proteins_all.CHS.hmmsearch
```

## CHI and CHI-like
Use CHI.fasta provided by Kevin and CHI HMM from eggNOG ENOG410IJZH.
```
blastp -query Hornwort_proteins_all.faa -db ~/blastdbs/flavonoids/CHI -evalue 1e-10 -outfmt 6 -out Hornwort_proteins_all.CHI.blastp

hmmsearch ~/blastdbs/flavonoids/CHI_eggnog.hmm.txt Hornwort_proteins_all.faa > Hornwort_proteins_all.CHI.hmmsearch
```


# Results

## PAL
Conserved across plant lineages, all hornwort hits in single orthogroup (OG0001366 hornworts only, OG0000314 Charophyte to Angiosperm). Appears to be 2-copy in all hornworts but Notothylas. Hornwort clades relatively distant from other bryophytes and land plant phylogeny not well reconstructed from orthogroup.

<img src="/images/PAL_tree.png">

---------------------------

## C4H
Hits in ~260 orthogroups containing >7000 sequences!!! Most with low sequence similarity to references from Kevin, consistent with belowing to large Cytochrome P450 family.  

<img src="/images/C4H_percSimilarityDist.png">

Selecting just those >60% based on distribution below, which fall into 4 OGs with 41 sequences. One group of hornwort seqs nicely matches C3H refs. 

<img src="/images/C3H_C4H_alignment.png">

Two possibly misannotated seqs from A. punctatus and A. agrestis link one group of seqs matching C4H refs, and another group of hornwort seqs matching extra 5' bit of sequence. The non-C4H portion of these proteins has closest BLAST hits to fanconi-associated nuclease 1 homolog (FAN1, nuclease required for the repair of DNA interstrand cross-links) with 40-50% identity.

<img src="/images/C4H_alignment.png">

For C4H "core" consisting of only proteins with ~global alignment to refs, there still seems to be a outgroup of hornwort proteins to C4H refs. In C4H ingroup, Leiosporoceros protein very truncated and only ~20% of sequence aligns. In both hornwort clades, topology mostly matches species tree except incorrect position of Leio in both clades.

<img src="/images/C4H_core_alignment.png">

<img src="/images/C4H_core_tree.png">

-----------------------------

## 4CL
Hits in 44 orthogroups containing 930 sequences. Most with low sequence similarity (median 27%) to references from Kevin. 

<img src="/images/4CL_percSimilarityDist.png">

4 OGs with 86 sequences have similarity >50% to references. Alignment of refs to all seqs shows no clear alliance between refs and certain orthogroups.

<img src="/images/4CL_alignment.png">

Aligning refs to each OG, OG0009039 seems to be best match with 67% pairwise similarity. However, phylogenetic position makes this hornwort clade appear to be an outgroup to 4CL, with hornworts more similar to Arabidopsis than Marchantia. 

<img src="/images/4CL_OG0009039_alignment.png">
<img src="/images/4CL_OG0009039_tree.png">



```python

```
