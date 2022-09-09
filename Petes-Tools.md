## Pete's Tools
**Various scripts I've written for little tasks. All located on server at /home/ps997/scripts/**

### getFromFasta.py
Use to get retrieve certain sequences from an FASTA file and return in FASTA format. Input is a string, file with one sequence name per line, or read from STDIN. Output goes to STDOUT, and can be redirected to a file. Any sequence names in the "retrieve" list that are not found in the source file are output to STDERR. Basic example returns Sequence1:
```
getFromFasta.py source.fasta "Sequence1"

>Sequence1
MAALAPVHLIGGNGGGDFYINGAHKGAILKRIRVWVGGWMIRGIEVQLSDGESQMFGSVDGGAREFTFQIGEKIAW
```

**Input file example**
An example from a source file, seqlist.txt, to get Sequence1, Sequence2, and Sequence3:
```
Sequence1
Sequence2
Sequence3
```
```
getFromFasta.py source.fasta seqlist.txt

>Sequence1
MAALAPVHLIGGNGGGDFYINGAHKGAILKRIRVWVGGWMIRGIEVQLSDGESQMFGSVDGGAREFTFQIGEKIAW
>Sequence2
MTSWPLKQEYAIEVASGLCVGLRGRAGADIDALGLTFLLPISHARLTNVRYPTLQLEAASIRPISIHEFYDENLSD
>Sequence3
MSIIYTPVHIIGSSSGGIVFNYDAAQSGGVLRRIGVWAGEWQLRGIRVWFTHTANPQTFGTANVGSYKEFEFTDGE
```

**STDIN pipeline example**
Can pipe a list of sequence names into the command, for example to get references from a BLAST tab-separated format output, blast.tsv:
```
Sequence1 Reference1     100.000 391     0       0       1       391     1       391     0.0     793
Sequence2 Reference2      100.000 391     0       0       1       391     1       391     0.0     793
Sequence3 Reference3     100.000 284     0       0       1       284     1       284     0.0     578
Sequence4 Reference4     100.000 391     0       0       1       391     1       391     0.0     793
Sequence5 Reference5     100.000 401     0       0       1       401     1       401     0.0     816
```
```
cut -f 2 blast.tsv | getFromFasta.py references.fasta - > references.blasthits.fasta
```


### parse_eggnog.sh
Use to retrieve eggNOG results for hornwort proteins. Input is a single transcript ID or a file with one transcript ID per line, and `ALL` to print whole line, or 1+ of these terms:
* SEED_ORTHOLOG
* SEED_EVALUE
* SEED_SCORE,
* BEST_TAX_LEVEL
* PREFERRED_NAME
* GO
* EC
* KEGG_KO
* KEGG_PATHWAY
* KEGG_MODULE
* KEGG_REACTION
* KEGG_RCLASS
* BRITE
* KEGG_TC
* CAZY
* BIGG
* TAXONOMY
* EGGNOG_OG
* BEST_EGGNOG_OG
* COG
* FREE_TEXT

Output goes to STDOUT and is tab-separated.


**Example 1**
Input:
```
/home/ps997/scripts/parse_eggnog.sh AagrBONN_evm.model.Sc2ySwM_1.5.1 ALL


```
Output:
```
AagrBONN_evm.model.Sc2ySwM_1.5.1	3067.XP_002950474.1	1.1e-49	203.4	Eukaryota													Eukaryota	COG2801@1,KOG0017@2759,KOG1075@1,KOG1075@2759	NA|NA|NA	E	Ribonuclease H protein


```
**Example 2**
Input:
```
/home/ps997/scripts/parse_eggnog.sh Leiosporoceros_dussii_v1_contig102_g113600.t1 PREFERRED_NAME FREE_TEXT


```
Output:
```
Leiosporoceros_dussii_v1_contig102_g113600.t1	RBR1	Regulator of biological processes that recruits a histone deacetylase to control gene transcription. May play a role in the entry into mitosis, negatively regulating the cell proliferation. Formation of stable complexes with geminiviridae replication-associated proteins may create a cellular environment which favors viral DNA replication


```

**Example 3**
Input from file:
```
parse_eggnog.sh names.txt PREFERRED_NAME FREE_TEXT

```
Output:
```
Phymatoceros_bulbiculosus_v1_contig3_g79550.t1		serine threonine-protein kinase
Phymatoceros_bulbiculosus_v1_contig3_g74910.t1	RAF1	Rubisco accumulation factor
Phymatoceros_bulbiculosus_v1_contig3_g79990.t1		ATP transmembrane transporter activity
Phymatoceros_bulbiculosus_v1_contig3_g69660.t1		HMGL-like
Phymatoceros_bulbiculosus_v1_contig3_g72050.t1		Pentatricopeptide repeat-containing protein
Phymatoceros_bulbiculosus_v1_contig3_g67680.t1		conserved protein UCP012943
```
