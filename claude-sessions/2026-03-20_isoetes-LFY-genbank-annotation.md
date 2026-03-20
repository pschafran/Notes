# Claude Code Session - 2026-03-20

**Project:** isoetes_LFY_ONT — Kizzort Branch paper GenBank submission
**Location:** `/media/data/projects/isoetes_LFY_ONT/analysis/david_wickell_phylogeny/brunton_samples_20250308/kizzort_branch_paper/`

## Summary

Prepared GenBank submission annotations for 127 new *Isoetes* LFY sequences from the Kizzort Branch phylogenetic paper. Used `I_engelmannii_AY541781.1.gb` as the annotation template. Generated feature tables (exon 2 / intron 2 / CDS / mRNA / gene) for all sequences, with per-sequence `codon_start` determined by pairwise alignment to the reference. Also added CLAUDE.md files to four levels of the project directory tree.

## Work Completed

### Files Modified
- (none — all new files)

### Files Created
- `kizzort_branch_paper/annotate_lfy.py` — main annotation script (see scripts directory)
- `kizzort_branch_paper/new_sequences.fsa` — FASTA with source qualifier deflines for tbl2asn (127 seqs; placeholder fields for voucher/country/date)
- `kizzort_branch_paper/new_sequences.tbl` — tbl2asn feature table (gene, mRNA, CDS, exon, intron)
- `kizzort_branch_paper/new_sequences_annotation.tsv` — per-sequence summary: exon_end, intron_start, codon_start, translated protein, pass/fail
- `kizzort_branch_paper/annotation_issues.tsv` — 2 sequences that could not receive full annotations
- `kizzort_branch_paper/CLAUDE.md`
- `brunton_samples_20250308/CLAUDE.md`
- `brunton_samples_20250308/combined_runs/CLAUDE.md`
- `isoetes_LFY_ONT/CLAUDE.md`

### Commands Executed
```bash
cd .../kizzort_branch_paper
python3 annotate_lfy.py
```

## Key Decisions & Rationale

- **Identifying new sequences**: Filtered `Kizzort_Branch_paper.fasta` (179 total) for sequences lacking a GenBank accession pattern (`[A-Z]{2}\d{6}\.\d`) → 127 new sequences.
- **Exon/intron boundary**: Located by pairwise alignment to AY541781.1 reference (exon 2 ends at ref pos 61, intron 2 starts at pos 62 with canonical GT). Each new sequence aligned individually; GT splice donor searched ±8 bp of mapped position.
- **codon_start derivation**: Mapped reference position 2 (first codon start) into each query via alignment. Formula: `((first_codon_pos - 1) % 3) + 1`. Fallback: `((intron_start - 62) % 3) + 2` when ref pos 2 falls outside the query.
- **Translation verification**: All 125 full-annotation sequences checked for C-terminal tail `CPTK`; all pass. Proteins vary at N-terminus depending on how much of the upstream codon(s) the amplicon captures.
- **MS-08 and NC-05**: These two sequences start entirely within intron 2 (confirmed from full alignment; first non-gap column 124 vs exon end column 80). Annotated with `gene` + `intron` only; no CDS/exon. User should confirm whether to submit.

## Technical Details

**LFY gene structure** (AY541781.1 reference):
- Exon 2: `<1..61` (partial at 5′; encodes partial LFY SAM-domain protein)
- Intron 2: `62..>1041` (GT..AG; partial at 3′)
- CDS codon_start=2; protein = `RFLEEVQHICRERGEKCPTK`

**codon_start distribution** across 125 annotated sequences:
- `3` — 102 sequences (amplicon starts 2 bp after reference; protein = `FLEEVQHICRERGEKCPTK`)
- `1` — 15 sequences (Brunton samples; amplicon starts 1 bp before reference; protein = `GRFLEEVQHICRERGEKCPTK`)
- `2` — 8 sequences (same start as reference; protein = `RFLEEVQHICRERGEKCPTK` or `GRFLEEVQHICRERGEKCPTK`)

The extra G (GGG codon) and/or missing R at the N-terminus across groups simply reflects slightly different 5′ amplicon boundaries — all are biologically correct partial CDS annotations.

## Remaining TODO for GenBank Submission
- [ ] Fill in `[specimen-voucher]`, `[country]`, `[collection-date]` in `new_sequences.fsa`
- [ ] Confirm species IDs for Brunton21341 / 21342 / 21343 / 21345
- [ ] Confirm whether MS-08 and NC-05 (intron-only) should be submitted
- [ ] Clarify "Isoetes Leary" (Schafran114) species identity
- [ ] Clarify *Isoetes* sp. (Schafran103_1 and _2) species identity
- [ ] Decide on inclusion of `_maybecontam` and `_low` labeled sequences
- [ ] Complete tbl2asn template (.sbt file) with submitter info and citation

## Related Files
- Reference GenBank record: `kizzort_branch_paper/I_engelmannii_AY541781.1.gb`
- Metadata template: `kizzort_branch_paper/genbank_metadata_template.csv`
- Full dataset FASTA: `kizzort_branch_paper/Kizzort_Branch_paper.fasta`
- IQ-TREE result: `kizzort_branch_paper/Kizzort_Branch_paper.alignment.trim.fasta.contree`

## Tags
`#isoetes` `#phylogenetics` `#genbank` `#ONT` `#LFY` `#annotation`
