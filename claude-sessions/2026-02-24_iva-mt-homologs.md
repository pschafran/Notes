# Claude Code Session - 2026-02-24

**Project:** iva_phylogeny/analysis — Mitochondrial Homolog Analysis
**Location:** `/media/data/projects/iva_phylogeny/analysis`
**Duration:** ~2026-02-24 morning session

## Summary
Implemented a full pipeline to find homologous mitochondrial sub-sequences across 23 *Iva*
samples: collected and labelled all mt scaffold sequences, ran all-vs-all BLASTN, clustered
scaffolds into homolog groups via a BLAST graph, and aligned the top cluster with MAFFT.
Also surveyed all samples for missing mt and nr assemblies.

## Work Completed

### Files Created
- `phylogeny_paper/mt/collect_mt.sh` — collects mt scaffolds from 23 samples, labels headers `>SAMPLE_ID|original_header`, outputs `all_mt_scaffolds.fasta` and `mt_manifest.tsv`
- `phylogeny_paper/mt/analyze_mt_homologs.py` — BLAST graph clustering (Union-Find), per-cluster stats, representative extraction, MAFFT alignment, PI-site counting, summary report

### Files Generated (by scripts)
- `phylogeny_paper/mt/all_mt_scaffolds.fasta` — 317 sequences from 23 samples (5.3 MB)
- `phylogeny_paper/mt/mt_manifest.tsv` — per-sample sequence counts and total bp
- `phylogeny_paper/mt/mt_blastdb.*` — BLAST nucleotide database
- `phylogeny_paper/mt/mt_allvsall.blastn` — 69,223 BLAST hits (15 MB)
- `phylogeny_paper/mt/mt_clusters.tsv` — per-scaffold cluster assignment
- `phylogeny_paper/mt/mt_cluster_stats.tsv` — per-cluster statistics
- `phylogeny_paper/mt/mt_summary.txt` — printed summary report
- `phylogeny_paper/mt/cluster_01/cluster_01_seqs.fasta` — 23 representatives (one longest per sample)
- `phylogeny_paper/mt/cluster_01/cluster_01_seqs_for_align.fasta` — 20 representatives (3 repeat-path samples excluded)
- `phylogeny_paper/mt/cluster_01/cluster_01_aligned.fasta` — **MAFFT alignment: 20 taxa, 118,538 bp, 42.97% PI sites**

### Commands Executed
```bash
# Step 1 – Collect sequences
bash phylogeny_paper/mt/collect_mt.sh

# Step 2 – BLAST
makeblastdb -in all_mt_scaffolds.fasta -dbtype nucl -out mt_blastdb
blastn -query all_mt_scaffolds.fasta -db mt_blastdb \
  -outfmt "6 qseqid sseqid pident length qlen slen qstart qend sstart send evalue bitscore" \
  -perc_identity 70 -evalue 1e-10 -num_threads 32 -out mt_allvsall.blastn

# Step 3 – Cluster and align
python3 phylogeny_paper/mt/analyze_mt_homologs.py
```

## Key Results

### Clustering
- **1 core cluster (cluster_01):** all 265/317 mt scaffolds from all 23 samples — expected, since all scaffolds are fragments of the same organellar genome
- **52 singleton clusters:** mostly short sequences (87–7,465 bp) specific to one sample (likely NUMTs or short unique mt regions)
- At any coverage threshold (50%–100% of samples), there is exactly 1 cross-sample cluster

### Alignment — `cluster_01_aligned.fasta`
| Property | Value |
|---|---|
| Taxa | 20 (3 repeat-path samples excluded) |
| Alignment length | 118,538 bp |
| Parsimony-informative sites | 50,937 (42.97%) |
| Excluded samples | axillaris_89, frutescens_36, frutescens_49 |

The 3 excluded samples have repeat-path assemblies where GetOrganelle produced a single large
"path" sequence (236–343 kb) instead of individual scaffolds — these can't be globally aligned
with individual scaffolds from other samples.

## Key Decisions & Rationale

- **Single path representative for repeat-path samples:** For axillaris_89, frutescens_36,
  frutescens_49 (which have only `repeat_pattern*.path_sequence.fasta` files), used
  `repeat_pattern1.1` as the representative for clustering but excluded these from MAFFT
  alignment (their sequences are 236–343 kb vs. ~10–73 kb for normal scaffolds).

- **BLAST filters:** aln_len ≥ 500 bp AND ≥ 50% coverage of shorter sequence. This yielded
  6,490 graph edges from 66,182 cross-sample hits.

- **Union-Find clustering:** Simple connected-component clustering. All mt scaffolds connect
  transitively (scaffold A overlaps B, B overlaps C → all in one cluster), which is the
  correct biological result.

## Mt Assembly Survey (all samples)

### Missing mt (empty dir — GetOrganelle ran but failed):
- frutescens_37, imbricata_54, imbricata_61, imbricata_69
- microcephala_25, microcephala_70, texensis_50, xanthifolia_88

### No mt directory at all (reads exist, not run):
- cultigen_13-3-53

### Repeat-path only (no clean single-path assembly):
- axillaris_89 (216 paths), frutescens_36 (360 paths), frutescens_49 (40 paths)

### Samples with mt but NOT included in the 23-sample analysis:
- annua_51, axillaris_91, frutescens_66, frutescens_83

## Nr Assembly Survey (all samples)

### Missing nr entirely (reads exist, not run):
- cultigen_13-3-53

### Repeat-path only (effectively unusable):
- hayesiana_93 — 1,000 paths / 72 patterns (severe repeat issue, known)
- angustifolia_77 — 4 paths / 4 patterns (mild; may be workable by picking one path)

### Samples with elevated nr path counts (multiple rDNA copies):
- frutescens_86: 8 path files
- axillaris_91: 6 path files
- microcephala_57: 4 path files

## Next Steps
- [ ] Run IQ-TREE on `phylogeny_paper/mt/cluster_01/cluster_01_aligned.fasta` (20-sample mt phylogeny)
- [ ] Decide whether to add annua_51, axillaris_91, frutescens_66, frutescens_83 to the paper sample set
- [ ] Retry mt assembly for failed samples (frutescens_37, imbricata_54/61/69, microcephala_25/70, texensis_50)
- [ ] Run GetOrganelle mt+nr for cultigen_13-3-53
- [ ] Populate `phylogeny_paper/cp/` and `phylogeny_paper/nr/` with paper sample sequences

## Technical Details

### Sample file path patterns (mt)
- EREMID batch: `<species>/<num>_EREMID_240811/embplant_mt/embplant_mt.K115.scaffolds.graph1.1.path_sequence.fasta`
- STURM batch (most): `<species>/<num>_STURM_251201adq30ft/embplant_mt/embplant_mt.K115.scaffolds.graph1.1.path_sequence.fasta`
- STURM batch (some): `<species>/<num>_STURM_251201adq30ft/<species>_<pop>_<num>_mt/embplant_mt.K115.scaffolds.graph1.1.path_sequence.fasta`
- Special: angustifolia_46 is at top level `46_EREMID_240811/`; xanthifolia_64 dir spelled `xanthiifolia_CO1_64_mt/` (double i)
- frutescens_58 uses K85 (not K115)

### BLAST parameters
```
-perc_identity 70 -evalue 1e-10
-outfmt "6 qseqid sseqid pident length qlen slen qstart qend sstart send evalue bitscore"
```

## Related Files
- Previous session: `2026-02-03_iva-chloroplast-phylogeny.md`
- Scripts: `~/Notes/claude-scripts/2026-02-24_iva-mt-homologs_scripts/`

## Tags
`#iva` `#phylogeny` `#mitochondria` `#BLAST` `#MAFFT` `#python`
