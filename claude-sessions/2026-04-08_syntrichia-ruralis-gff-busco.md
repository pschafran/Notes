# Claude Code Session - 2026-04-08

**Project:** hornwort_sex_chromosomes — Syntrichia ruralis GFF processing & BUSCO
**Location:** `/media/data/projects/hornwort_sex_chromosomes/analysis/moss_genomes/Syntrichia_ruralis/`
**Duration:** ~evening session (continued from prior context)

## Summary

Resolved a non-standard "linked-list exon" encoding in the CoGe GFF for *Syntrichia ruralis*, produced a correctly structured primary-transcript GFF3, extracted a primary protein FASTA, and ran BUSCO. Also completed CLAUDE.md documentation updates for this species and the parent moss_genomes directory. BUSCO completeness is 72.7%, judged insufficient for OrthoFinder — reannotation from scratch is planned.

## Work Completed

### Files Created
- `analysis/moss_genomes/Syntrichia_ruralis/flatten_v1a_chains.py` — Converts CoGe V1A chain-encoded genes to standard GFF3
- `analysis/moss_genomes/Syntrichia_ruralis/Syntrichia_ruralis_annos1-cds0-id_typename-nu1-upa1-add_chr0.gid66749_renamed.fixed.gff` — Corrected GFF3; 16,064 genes / 16,064 mRNAs
- `analysis/moss_genomes/Syntrichia_ruralis/Syntrichia_ruralis.PROT_primary_from_GFF.fa` — Primary protein FASTA from fixed GFF; 16,064 sequences
- `analysis/moss_genomes/Syntrichia_ruralis/CLAUDE.md` — New species-level documentation

### Files Modified
- `analysis/moss_genomes/CLAUDE.md` — Added note that Syrur CoGe annotation is problematic and reannotation is planned

### Commands Executed
```bash
# Flatten V1A chains
python3 flatten_v1a_chains.py

# Extract primary proteins
gffread -H -V -y Syntrichia_ruralis.PROT_primary_from_GFF.fa \
  -g Syntrichia_ruralis_renamed.faa \
  Syntrichia_ruralis_annos1-cds0-id_typename-nu1-upa1-add_chr0.gid66749_renamed.fixed.gff

# BUSCO
conda run -n busco busco -m proteins \
  -i Syntrichia_ruralis.PROT_primary_from_GFF.fa \
  -o Syntrichia_ruralis.PROT_primary_from_GFF.fa.busco \
  --out_path .../busco_analysis/ \
  -l .../viridiplantae_odb12/ --offline -c 16
```

## Key Decisions & Rationale

- **Did not add protein FASTA to OrthoFinder input:** BUSCO at 72.7% is below acceptable threshold; the CoGe annotation is structurally unusual and of uncertain quality. User plans to reannotate with BRAKER or similar.
- **Used `_renamed.fixed.gff` not `_renamed.primary.gff`:** The earlier `primary.gff` (gawk-filtered) only kept the first exon of V1A-chain genes, producing truncated proteins and 66.7% BUSCO. The fixed GFF recovers these to 72.7%.

## Technical Details

### CoGe GFF — Two Gene Classes
The CoGe GFF (genome ID 66749) encodes genes in two formats:

**Class 1 — Plain genes (10,450, ~65%):**
Standard GFF3: `gene → mRNA → exon/CDS`

**Class 2 — V1A chain genes (5,614, ~35%):**
Non-standard linked-list: `gene(Sr_gN_V1A) → mRNA1(exon1) → mRNA2(exon2) → ...`
- Each "mRNA" node contains exactly one exon/CDS
- Each node's Parent is the previous mRNA (not the gene)
- All nodes share the same `coge_fid`
- Chain is strictly linear (no branching)
- `flatten_v1a_chains.py` resolves these to standard GFF3

### BUSCO Comparison (viridiplantae_odb12)
| Input | Complete | S | D | Fragmented | Missing |
|---|---|---|---|---|---|
| All isoforms (CoGe original) | 81.6% | 77.5% | 4.1% | 12.3% | 6.1% |
| primary.gff (chain-broken) | 66.7% | 64.2% | 2.4% | 12.0% | 21.3% |
| fixed.gff (chains flattened) | 72.7% | 69.5% | 3.3% | 12.2% | 15.1% |

The all-isoform score is inflated. The 15.1% missing BUSCOs in the fixed primary set reflects genuine annotation incompleteness.

## Next Steps
- [ ] Reannotate *Syntrichia ruralis* gene models from scratch (BRAKER, using RNA-seq or protein evidence)
- [ ] After reannotation: run BUSCO, add to OrthoFinder viridiplantae input, run composition stats against new gene models
- [ ] Broader audit follow-up: BUSCO and renamed proteins still missing for ~10 moss/liverwort species (see session 2026-04-08 morning notes)

## Related Files
- `analysis/moss_genomes/Syntrichia_ruralis/CLAUDE.md`
- `analysis/moss_genomes/CLAUDE.md`
- Previous session notes: `2026-04-08_orthofinder-audit-gff-processing.md` (if it exists)

## Tags
`#moss-genomes` `#gff-processing` `#busco` `#syntrichia` `#python`
