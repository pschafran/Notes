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
Very conserved across plant lineages, all hornwort hits in single orthogroup (OG0001366 hornworts only, OG0000314 Charophyte to Angiosperm). Appears to be 2-copy in all hornworts but Notothylas.  

## C4H
Hits in ~260 orthogroups containing >7000 sequences!!! Most with low sequence similarity to references from Kevin. 

<img src="/pschafran/Notes/C4H_percSimilarityDist.png">





```python

```
