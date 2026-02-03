# Claude Code Session - 2026-02-03

**Project:** iva_phylogeny/analysis
**Location:** `/media/data/projects/iva_phylogeny/analysis`

## Summary
Processed 34 chloroplast genomes from multiple *Iva* species assembled by GetOrganelle. Identified and corrected LSC inversions in 4 samples (01, 48, frutescens_58, hayesiana_92) and a circular rotation issue in imbricata_69. Generated a multi-species chloroplast phylogeny using IQ-TREE.

## Work Completed

### Files Created
- `chloroplast_forward_combined.fasta` - Initial combined file (34 samples, before fixes)
- `chloroplast_forward_combined_fixed.fasta` - Combined file with all orientation corrections
- `chloroplast_forward_aligned.fasta` - Initial MAFFT alignment (had misorientation artifacts)
- `chloroplast_forward_aligned_fixed.fasta` - Final corrected alignment (158,892 bp)
- `chloroplast_phylogeny.treefile` - Maximum likelihood tree (Newick format)
- `chloroplast_phylogeny.contree` - Bootstrap consensus tree
- `chloroplast_phylogeny.iqtree` - Full IQ-TREE analysis report

### Commands Executed
```bash
# Orientation analysis using tblastn
tblastn -query ndhD.faa -subject [fasta] -outfmt "6 sstart send"
tblastn -query rbcL.faa -subject [fasta] -outfmt "6 sstart send"

# Alignment
/home/peter/miniconda3/bin/mafft --auto --thread -1 chloroplast_forward_combined_fixed.fasta > chloroplast_forward_aligned_fixed.fasta

# Phylogenetic inference
/home/peter/miniconda3/bin/iqtree3 -s chloroplast_forward_aligned_fixed.fasta -m MFP -B 1000 -T AUTO -pre chloroplast_phylogeny
```

## Key Decisions & Rationale

- **Decision:** Used ndhD gene orientation (SSC region) as primary classifier for isoform selection
  - **Rationale:** ndhD flips between the two natural chloroplast isoforms; provides consistent reference orientation

- **Decision:** Used rbcL gene orientation (LSC region) to detect LSC inversions
  - **Rationale:** When ndhD is forward but rbcL is also forward (instead of reverse), indicates LSC inversion

- **Decision:** Fixed LSC inversions by reverse complementing the LSC region (~83-84kb)
  - **Rationale:** Maintains SSC/IR orientation while correcting the inverted LSC to match majority

- **Decision:** Fixed imbricata_69 by circular rotation rather than reverse complement
  - **Rationale:** Gene orientations were correct but genome was linearized at different start position

- **Decision:** Did NOT modify frutescens_49 despite elevated divergence
  - **Rationale:** Gene orientations correct; ~850bp size difference indicates genuine IR expansion variant

## Technical Details

### Sample Classification by Orientation

**Samples requiring LSC inversion fix (ndhD forward, rbcL forward):**
| Sample | LSC boundary | Fix applied |
|--------|-------------|-------------|
| 01 | 83,675 bp | Reverse complement LSC |
| 48 | 83,659 bp | Reverse complement LSC |
| frutescens_58 | 83,740 bp | Reverse complement LSC |
| hayesiana_92 | 83,708 bp | Reverse complement LSC |

**Sample requiring rotation fix:**
| Sample | Issue | Fix applied |
|--------|-------|-------------|
| imbricata_69 | Scaffolds assembly, different start position | Rotated 79,565 bp |

### Divergence Analysis Results (after fixes)
All samples now show <1% maximum divergence vs reference (frutescens_36), except:
- frutescens_49: 7.4% max (genuine IR expansion variant)
- xanthifolia samples: 3.2% max (inter-species divergence)

### Phylogenetic Analysis
- **Best-fit model:** K3Pu+F+I+G4
- **Invariant sites:** 63.5%
- **Alignment length:** 158,892 bp
- **Parsimony informative sites:** 3,128
- **Bootstrap replicates:** 1,000

### Key Phylogenetic Findings
1. frutescens_49 has extremely long branch (0.034 subs/site) due to IR expansion
2. *I. xanthifolia* (64, 88) forms distinct clade with ~0.9% divergence
3. *I. imbricata* samples cluster with sample 54
4. *I. hayesiana* (92, 93) forms well-supported clade (100% bootstrap)
5. Core *I. frutescens* samples nearly identical (36, 37, 58, 63, 83, 86)

## Challenges & Solutions

**Problem:** 4 samples (01, 48, 58, hayesiana_92) showed ~50% divergence in first ~110kb of alignment
**Solution:** Identified as LSC inversions using rbcL orientation; reverse complemented LSC region in each sample

**Problem:** imbricata_69 showed ~28% max divergence in first ~45kb
**Solution:** Identified as circular rotation issue (scaffolds assembly started at different position); rotated sequence to match reference start position

**Problem:** frutescens_49 showed moderate divergence (~7%)
**Solution:** Determined to be genuine biological variant with ~850bp IR expansion; no correction needed

## Scripts Created
- `analyze_chloroplast_orientations.py` - Analyzes ndhD/rbcL orientation for all samples
- `combine_chloroplast_forward.py` - Combines forward isoforms into single file
- `fix_lsc_inversion.py` - Detects and fixes LSC inversions
- `fix_rotation.py` - Fixes circular rotation in imbricata_69
- `regenerate_combined.py` - Regenerates combined file with all fixes
- `analyze_alignment_divergence.py` - Sliding window divergence analysis

## Next Steps
- [ ] Investigate frutescens_49 IR expansion in detail
- [ ] Consider removing frutescens_49 from phylogeny or using gene-based analysis
- [ ] Root the tree with appropriate outgroup
- [ ] Compare chloroplast phylogeny with nrDNA phylogeny (see 2026-02-01 session)
- [ ] Annotate chloroplast genomes to extract individual genes for gene tree analysis

## Related Files
- Previous session: `~/Notes/claude-sessions/2026-01-30_iva-frutescens-chloroplast.md`
- nrDNA session: `~/Notes/claude-sessions/2026-02-01_iva-nrDNA-phylogeny.md`
- Reference proteins: `ndhD.faa`, `rbcL.faa`

## Tags
`#chloroplast` `#phylogenetics` `#getorganelle` `#iva` `#genome-rearrangement` `#iqtree` `#mafft` `#lsc-inversion`
