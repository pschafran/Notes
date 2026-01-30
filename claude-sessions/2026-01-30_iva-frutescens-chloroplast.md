# Claude Code Session - 2026-01-30

**Project:** iva_phylogeny/analysis/frutescens
**Location:** `/media/data/projects/iva_phylogeny/analysis/frutescens`

## Summary
Analyzed GetOrganelle chloroplast genome assemblies from 9 *Iva frutescens* samples. Classified isoforms by ndhD gene orientation, generated a phylogenetic tree, and discovered a large-scale LSC inversion in sample 58 that was causing artifactual long-branch attraction.

## Work Completed

### Files Created
- `chloroplast_forward.fasta` - Multi-FASTA with 9 chloroplast genomes (ndhD forward orientation)
- `chloroplast_reverse.fasta` - Multi-FASTA with 9 chloroplast genomes (ndhD reverse orientation)
- `chloroplast_forward_aligned.fasta` - MAFFT alignment of forward-orientation genomes
- `chloroplast_forward_tree.*` - IQ-TREE output files (treefile, contree, iqtree, log, etc.)

### Commands Executed
```bash
# Find GetOrganelle output directories
find . -type d -name "*_pt"

# Determine ndhD orientation using tblastn
tblastn -query ndhD.faa -subject [fasta] -outfmt "6 sstart send"

# Align chloroplast genomes
mafft --auto --thread -1 chloroplast_forward.fasta > chloroplast_forward_aligned.fasta

# Build phylogenetic tree with model selection and bootstrap
iqtree -s chloroplast_forward_aligned.fasta -m MFP -bb 1000 -nt AUTO -pre chloroplast_forward_tree

# Investigate sample 58 rearrangement
blastn -query seq58_first100k.fasta -subject seq36_full.fasta -outfmt "6 qstart qend sstart send pident length"
```

## Key Decisions & Rationale
- **Decision:** Used ndhD gene orientation to classify isoforms as "forward" or "reverse"
  - **Rationale:** ndhD is located in the SSC region which flips between the two natural chloroplast isoforms; provides consistent reference orientation
- **Decision:** Used IQ-TREE with ModelFinder (MFP) for phylogenetic inference
  - **Rationale:** Automatic model selection via BIC ensures appropriate substitution model; ultrafast bootstrap (1000 replicates) provides branch support
- **Decision:** Investigated sample 58's long branch rather than excluding it
  - **Rationale:** Long branches can indicate assembly errors, contamination, or genuine biological differences worth understanding

## Technical Details

### Sample Isoform Classification
Based on ndhD orientation (tblastn sstart < send = forward):

| Sample | Forward Isoform | Reverse Isoform |
|--------|-----------------|-----------------|
| 36, 37, 63, 71, 94 | graph1.1 | graph1.2 |
| 49, 58, 68, 86 | graph1.2 | graph1.1 |

### Phylogenetic Analysis
- **Best-fit model:** K3Pu+F+G4 (selected by BIC)
- **Alignment length:** 170,779 bp
- **Notable topology:** Sample 58 showed extremely long branch (0.32 subs/site)

### LSC Inversion Discovery in Sample 58
Sliding window analysis of alignment revealed:
- Positions 0 - ~100 kb: ~47% divergence (POOR alignment)
- Positions ~102 kb - 170 kb: 0% divergence (GOOD alignment)

BLASTN confirmed ~84 kb LSC inversion:
```
Sample 58: 1-83,740 → Sample 36: 83,741-1 (REVERSE, 99.97% identity)
Sample 58: 83,741-152,182 → Sample 36: 83,742-152,173 (Forward, 100% identity)
```

### Chloroplast Genome Structure
```
Standard:    5'─[──────LSC──────→][IR-A][SSC→][IR-B]─3'
Sample 58:   5'─[←──────LSC──────][IR-A][SSC→][IR-B]─3'
```

## Challenges & Solutions
**Problem:** Sample 58 showed extremely long branch in phylogenetic tree
**Solution:** Sliding window analysis of alignment identified ~100 kb region with ~47% divergence; BLASTN confirmed this was due to LSC inversion, not true sequence divergence

**Problem:** GetOrganelle isoform numbering (1 vs 2) is arbitrary across samples
**Solution:** Used ndhD protein sequence as orientation reference via tblastn to consistently classify all samples

## Next Steps
- [ ] Reverse complement LSC region of sample 58 and rebuild phylogenetic tree
- [ ] Investigate whether LSC inversion is assembly artifact or genuine biological variant
- [ ] Consider using conserved gene sequences for phylogenetics instead of whole genomes
- [ ] Check if other samples have similar structural variants

## Related Files
- `/home/peter/miniconda3/envs/getorganelle` - GetOrganelle conda environment
- `ndhD.faa` - Reference ndhD amino acid sequence used for orientation

## Tags
`#chloroplast` `#phylogenetics` `#getorganelle` `#iva-frutescens` `#genome-rearrangement` `#iqtree` `#mafft`
