# Claude Code Session - 2026-01-30

**Project:** hornwort_sex_chromosomes
**Location:** `/media/data/projects/hornwort_sex_chromosomes/analysis/`
**Duration:** Morning session

## Summary

Completed BUSCO analysis on all `--missing ingroup` FASTA files (26 runs) and performed comprehensive comparison between the two alien index parameter settings (`--missing outgroup` vs `--missing ingroup`). Results demonstrate that `--missing ingroup` produces significantly better plant transcriptome quality and cleaner contamination/HGT candidate sets. Also renamed all sequences in final FASTA files with species/sample identifiers for cross-sample uniqueness.

## Work Completed

### BUSCO Analysis
- Ran BUSCO v6.0.0 on 26 FASTA files (13 samples × 2 AI categories)
- Used viridiplantae_odb12 lineage database (822 BUSCOs)
- 32 CPU cores per run
- All runs completed successfully (~8 minutes total)

### Files Created
For each of 13 samples:
- `busco_negative_ai_missing_ingroup/` - BUSCO results for negative AI sequences
- `busco_positive_ai_missing_ingroup/` - BUSCO results for positive AI sequences
- `transcripts.fasta.transdecoder.pep.negative_ai.missing_ingroup.renamed.fasta` - Renamed sequences
- `transcripts.fasta.transdecoder.pep.positive_ai.missing_ingroup.renamed.fasta` - Renamed sequences

### Sequence Renaming
Renamed all sequences in final `--missing ingroup` FASTA files with species/sample identifiers.

**Naming scheme:** 2-letter genus + 3-letter species + sample ID (e.g., `>Fofuc_A017|NODE_10000...`)

| Code | Full Species Name | Sample ID |
|------|-------------------|-----------|
| Anspe | Anthoceros sp | A012 |
| Dejav | Dendroceros javanicus | A016 |
| Fofuc | Folioceros fuciformis | A017 |
| Fokas | Folioceros kashyapii | A008 |
| Ledus | Leiosporoceros dussii | ANON |
| Noaen | Nothoceros aenigmaticus | DXOU |
| Novin | Nothoceros vincentianus | TCBC |
| Ntjav | Notothylas javanica | A002 |
| Ntorb | Notothylas orbicularis | A011 |
| Pahal | Paraphymatoceros hallii | FAJB |
| Phcar | Phaeoceros carolinianus | A001 |
| Pmcor | Phaeomegaceros coriaceus | AKXB |
| Pybul | Phymatoceros bulbiculosus | A014 |

**Total sequences renamed:**
- Negative AI: 379,752 sequences across 13 samples
- Positive AI: 77,153 sequences across 13 samples
- All sequence names verified unique across all samples

## Results

### Aggregated Summary

| Metric | --missing outgroup | --missing ingroup | Difference |
|--------|-------------------|-------------------|------------|
| **Negative AI avg BUSCO** | 51.6% | 72.8% | **+21.2%** |
| **Positive AI avg BUSCO** | 19.7% | 3.0% | **-16.7%** |

### Individual Sample Results: Negative AI (Plant-like Sequences)

| Sample | Seqs (outgroup) | Seqs (ingroup) | BUSCO (outgroup) | BUSCO (ingroup) | Change |
|--------|-----------------|----------------|------------------|-----------------|--------|
| Anthoceros_sp | 20,038 | 27,123 | 50.1% | 72.3% | +22.2% |
| Dendroceros_javanicus | 23,371 | 33,112 | 51.8% | 78.8% | +27.0% |
| Folioceros_fuciformis | 20,169 | 28,406 | 49.8% | 74.7% | +24.9% |
| Folioceros_kashyapii | 12,507 | 17,450 | 34.9% | 51.0% | +16.1% |
| Leiosporoceros | 18,920 | 27,471 | 54.7% | 80.0% | +25.3% |
| Nothoceros_aenigmaticus | 9,980 | 13,976 | 24.3% | 38.6% | +14.3% |
| Nothoceros_vincentianus | 17,277 | 24,337 | 50.2% | 74.2% | +24.0% |
| Notothylas_javanica | 19,455 | 26,820 | 70.4% | 81.5% | +11.1% |
| Notothylas_orbicularis | 30,453 | 41,705 | 54.0% | 79.3% | +25.3% |
| Paraphymatoceros_hallii | 20,647 | 28,742 | 54.4% | 80.4% | +26.0% |
| Phaeoceros_carolinianus | 28,656 | 39,776 | 59.6% | 86.1% | +26.5% |
| Phaeomegaceros_coriaceus | 16,400 | 23,406 | 47.9% | 74.5% | +26.6% |
| Phymatoceros_bulbiculosus | 33,470 | 47,428 | 68.2% | 75.8% | +7.6% |

### Individual Sample Results: Positive AI (Potential Contamination/HGT)

| Sample | Seqs (outgroup) | Seqs (ingroup) | BUSCO (outgroup) | BUSCO (ingroup) | Change |
|--------|-----------------|----------------|------------------|-----------------|--------|
| Anthoceros_sp | 8,210 | 2,379 | 20.3% | 0.7% | -19.6% |
| Dendroceros_javanicus | 10,291 | 2,579 | 22.1% | 0.6% | -21.5% |
| Folioceros_fuciformis | 9,857 | 3,265 | 19.2% | 1.3% | -17.9% |
| Folioceros_kashyapii | 9,516 | 5,231 | 18.1% | 2.4% | -15.7% |
| Leiosporoceros | 7,043 | 558 | 18.6% | 0.2% | -18.4% |
| Nothoceros_aenigmaticus | 3,687 | 366 | 11.2% | 0.0% | -11.2% |
| Nothoceros_vincentianus | 8,074 | 2,735 | 19.2% | 0.6% | -18.6% |
| Notothylas_javanica | 10,885 | 5,245 | 9.7% | 1.3% | -8.4% |
| Notothylas_orbicularis | 15,749 | 6,570 | 20.7% | 2.3% | -18.4% |
| Paraphymatoceros_hallii | 7,920 | 1,702 | 20.8% | 1.3% | -19.5% |
| Phaeoceros_carolinianus | 16,793 | 7,848 | 25.3% | 3.8% | -21.5% |
| Phaeomegaceros_coriaceus | 11,410 | 5,786 | 25.2% | 5.4% | -19.8% |
| Phymatoceros_bulbiculosus | 46,231 | 33,902 | 25.8% | 19.2% | -6.6% |

### Detailed BUSCO Metrics: Negative AI with --missing ingroup

| Sample | Complete | Single | Duplicated | Fragmented | Missing |
|--------|----------|--------|------------|------------|---------|
| Anthoceros_sp | 72.3% | 50.4% | 21.9% | 20.8% | 6.9% |
| Dendroceros_javanicus | 78.8% | 53.5% | 25.3% | 14.7% | 6.4% |
| Folioceros_fuciformis | 74.7% | 55.0% | 19.7% | 18.6% | 6.7% |
| Folioceros_kashyapii | 51.0% | 40.5% | 10.5% | 32.0% | 17.0% |
| Leiosporoceros | 80.0% | 50.1% | 29.9% | 13.6% | 6.3% |
| Nothoceros_aenigmaticus | 38.6% | 35.5% | 3.0% | 37.5% | 24.0% |
| Nothoceros_vincentianus | 74.2% | 42.8% | 31.4% | 18.4% | 7.4% |
| Notothylas_javanica | 81.5% | 59.4% | 22.1% | 12.5% | 6.0% |
| Notothylas_orbicularis | 79.3% | 45.7% | 33.6% | 14.5% | 6.2% |
| Paraphymatoceros_hallii | 80.4% | 57.5% | 22.9% | 13.6% | 6.0% |
| Phaeoceros_carolinianus | 86.1% | 54.9% | 31.3% | 8.3% | 5.6% |
| Phaeomegaceros_coriaceus | 74.5% | 58.8% | 15.7% | 18.0% | 7.5% |
| Phymatoceros_bulbiculosus | 75.8% | 45.4% | 30.4% | 16.5% | 7.7% |

## Key Findings

### 1. N/A Taxonomy Sequences are Predominantly Plant
- Sequences lacking taxonomy annotation in BLASTP output are mostly plant sequences
- When treated as ingroup (plants), BUSCO completeness improves by ~21% on average
- This is consistent across all 13 samples

### 2. --missing ingroup Produces Cleaner Results
- **Negative AI (plant-like)**: Higher BUSCO = more complete transcriptome
- **Positive AI (contamination)**: Lower BUSCO = fewer misclassified plant genes
- The positive AI set with <5% BUSCO likely represents genuine non-plant sequences

### 3. Notable Outliers
- **Phymatoceros_bulbiculosus**:
  - Smallest improvement in negative AI (+7.6%)
  - Highest BUSCO retention in positive AI (19.2%)
  - Has unusually high positive AI sequence count (33,902 with --missing ingroup)
  - May indicate genuine HGT events or unique contamination pattern

- **Nothoceros_aenigmaticus**:
  - Lowest overall BUSCO completeness (38.6%)
  - Highest fragmented + missing (61.5%)
  - May indicate lower quality transcriptome assembly

- **Folioceros_kashyapii**:
  - Second lowest BUSCO (51.0%)
  - High fragmented rate (32.0%)
  - Similar quality concerns

## Recommendation

**Use `--missing ingroup` for hornwort alien index analysis**

Rationale:
1. Produces higher quality plant transcriptome set (avg 72.8% vs 51.6% BUSCO)
2. Creates cleaner contamination/HGT candidate set (avg 3.0% vs 19.7% BUSCO)
3. Correctly classifies N/A taxonomy sequences as plant (ingroup)

## Commands Executed

```bash
# BUSCO runs (26 total)
source ~/miniconda3/etc/profile.d/conda.sh
conda activate busco
busco -m prot -l viridiplantae_odb12 -c 32 \
  -i transcripts.fasta.transdecoder.pep.negative_ai.missing_ingroup.fasta \
  -o busco_negative_ai_missing_ingroup --offline

busco -m prot -l viridiplantae_odb12 -c 32 \
  -i transcripts.fasta.transdecoder.pep.positive_ai.missing_ingroup.fasta \
  -o busco_positive_ai_missing_ingroup --offline

# Sequence renaming (for each sample)
sed "s/^>/>Fofuc_A017|/" input.fasta > output.renamed.fasta
```

## Discussion: Protein Set Improvement Methods

Discussed potential methods to further improve transcriptome-derived protein sets:

1. **Redundancy reduction** - CD-HIT/MMseqs2 clustering at 95% identity
2. **ORF quality filtering** - Keep complete ORFs, minimum length, coding potential scores
3. **Expression-based filtering** - Remove low-expressed transcripts
4. **Isoform selection** - Longest ORF or highest-expressed per gene
5. **Functional annotation filtering** - Keep sequences with Pfam domains or SwissProt hits
6. **Homology-based validation** - OrthoFinder for cross-species ortholog detection
7. **Assembly improvement** - Long-read data, error correction

Not implemented this session - reserved for future work.

## Next Steps

- [x] ~~Rename sequences with species/sample identifiers~~ (completed)
- [ ] Analyze which sequences differ between the two parameter settings
- [ ] Investigate Phymatoceros_bulbiculosus positive AI sequences (19.2% BUSCO - potential HGT)
- [ ] Characterize the positive AI sequences taxonomically (what organisms?)
- [ ] Consider re-running samples with low BUSCO (Nothoceros_aenigmaticus, Folioceros_kashyapii)
- [ ] Redundancy reduction with CD-HIT
- [ ] Ortholog detection with OrthoFinder

## Related Files

- Test suite: `/media/data/projects/hornwort_sex_chromosomes/analysis/test_pipeline.py`
- Previous session: `~/Notes/claude-sessions/2026-01-29_hornwort-sex-chromosomes-alien-index.md`

## Tags

`#hornwort` `#transcriptomics` `#alien-index` `#BUSCO` `#contamination-detection` `#HGT` `#parameter-comparison`
