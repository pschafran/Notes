# Claude Code Session - 2026-02-01

**Project:** iva_phylogeny/analysis
**Location:** `/media/data/projects/iva_phylogeny/analysis`

## Summary
Analyzed GetOrganelle chloroplast and nuclear ribosomal DNA assemblies from 12 *Iva* samples. Classified chloroplast genomes by rbcL/ndhD orientation, reverse-complemented sample 58 to correct LSC inversion, then performed nrDNA alignment and phylogenetic analysis revealing unexpected long branch in sample 49.

## Work Completed

### Files Created
- `nrDNA_combined.fasta` - Combined nrDNA sequences from 12 samples
- `nrDNA_aligned.fasta` - MAFFT alignment with auto-reorientation (13,922 positions)
- `nrDNA_aligned_trimmed.fasta` - End-trimmed alignment (6,893 positions)
- `nrDNA_tree.*` - IQ-TREE output files (treefile, contree, iqtree, splits.nex, log, model.gz, mldist)
- `frutescens/58_STURM_251201adq30ft/embplant_pt/embplant_pt.K115.complete.graph1.1.path_sequence.revcomp.fasta`
- `frutescens/58_STURM_251201adq30ft/embplant_pt/embplant_pt.K115.complete.graph1.2.path_sequence.revcomp.fasta`

### Commands Executed
```bash
# Classify chloroplast genomes by rbcL and ndhD orientation
tblastn -query rbcL.faa -subject [fasta] -outfmt "6 sstart send"
tblastn -query ndhD.faa -subject [fasta] -outfmt "6 sstart send"

# Reverse complement sample 58 chloroplast sequences
seqtk seq -r embplant_pt.K115.complete.graph1.1.path_sequence.fasta > embplant_pt.K115.complete.graph1.1.path_sequence.revcomp.fasta
seqtk seq -r embplant_pt.K115.complete.graph1.2.path_sequence.fasta > embplant_pt.K115.complete.graph1.2.path_sequence.revcomp.fasta

# Align nrDNA with automatic reorientation
mafft --auto --adjustdirection --thread -1 nrDNA_combined.fasta > nrDNA_aligned.fasta

# Build phylogenetic tree
iqtree -s nrDNA_aligned_trimmed.fasta -m MFP -bb 1000 -nt AUTO -pre nrDNA_tree
```

## Key Decisions & Rationale
- **Decision:** Used both rbcL and ndhD to classify chloroplast orientation
  - **Rationale:** rbcL is in LSC region, ndhD in SSC; dual classification reveals structural variants like sample 58's LSC inversion
- **Decision:** Used MAFFT --adjustdirection for nrDNA alignment
  - **Rationale:** Automatically detects and reverse-complements sequences in wrong orientation
- **Decision:** Trimmed alignment ends where any sequence had gaps
  - **Rationale:** Removes overhanging regions from length variation while preserving internal gaps
- **Decision:** Used graph1.1 as representative for samples with multiple nrDNA sequences
  - **Rationale:** Samples 86 (4 seqs) and 91 (2 seqs) have intragenomic variation; graph1.1 is consistent choice

## Technical Details

### Chloroplast rbcL/ndhD Orientation Classification

| Sample | graph1.1 | graph1.2 |
|--------|----------|----------|
| 89_STURM | R-R | R-F |
| 91_STURM | R-R | R-F |
| 36_frutescens | R-F | R-R |
| 37_frutescens | R-F | R-R |
| 49_frutescens | R-R | R-F |
| **58_frutescens** | **F-R** | **F-F** |
| 63_frutescens | R-F | R-R |
| 68_frutescens | R-R | R-F |
| 71_frutescens | R-F | R-R |
| 86_frutescens | R-R | R-F |
| 94_frutescens | R-F | R-R |

Sample 58 has Forward rbcL (LSC inversion), confirmed from previous session. Reverse-complemented sequences now show R-F/R-R pattern matching other samples.

### nrDNA Sequence Summary

| Sample | Length (bp) | Assembly Type | Notes |
|--------|-------------|---------------|-------|
| 01_EREMID | 10,265 | K115 complete | Single sequence |
| 89_STURM | 10,311 | K115 complete | Single sequence |
| 91_STURM | 10,752 | K115 complete | 2 sequences (14 SNPs in ITS) |
| 36-94_frutescens | 6,806-8,936 | K115 scaffolds/complete | Most single sequence |
| 86_frutescens | 8,928 | K55 complete | 4 sequences (high intragenomic variation) |

### nrDNA Alignment
- 5 sequences reverse-complemented by MAFFT: 37, 63, 68, 71, 94
- Trimmed 3,512 positions from start, 3,517 from end
- Final alignment: 6,893 positions

### Phylogenetic Analysis
- **Best-fit model:** TIM2+F+I (selected by BIC)
- **Bootstrap:** 1000 ultrafast replicates
- **Key topology:**
  - 89_STURM + 91_STURM cluster (BS=98)
  - Most frutescens form tight cluster (BS=99)
  - 49_frutescens has extremely long branch (0.11 subs/site)

## Challenges & Solutions
**Problem:** Sample 58 chloroplast has inverted LSC region affecting rbcL orientation
**Solution:** Reverse-complemented entire sequences using seqtk; verified orientations now match other samples

**Problem:** Samples 86 and 91 have multiple nrDNA sequences
**Solution:** Used graph1.1 as representative; documented that variation is intragenomic rDNA heterogeneity

**Problem:** nrDNA sequences in different orientations
**Solution:** MAFFT --adjustdirection automatically detected and fixed 5 sequences

## Next Steps
- [ ] Investigate sample 49's long branch in nrDNA tree (contamination? assembly error? genuine divergence?)
- [ ] Consider excluding sample 49 and rebuilding tree
- [ ] Extract ITS1/ITS2 regions for more focused phylogenetic analysis
- [ ] Compare chloroplast and nrDNA phylogenies for congruence
- [ ] Complete analysis of 01_EREMID chloroplast (assembly still in progress)

## Related Files
- Previous session: `/home/peter/Notes/claude-sessions/2026-01-30_iva-frutescens-chloroplast.md`
- Reference proteins: `rbcL.faa`, `ndhD.faa`
- GetOrganelle environment: `/home/peter/miniconda3/envs/getorganelle`

## Tags
`#nrDNA` `#ITS` `#phylogenetics` `#getorganelle` `#iva` `#mafft` `#iqtree` `#chloroplast`
