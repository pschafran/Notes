# Claude Code Session - 2026-03-24

**Project:** Lunularia_cruciata — four_assembly_cross
**Location:** `/media/data/projects/hornwort_sex_chromosomes/analysis/liverwort_genomes/Lunularia_cruciata/four_assembly_cross`
**Duration:** 2026-03-24

## Summary
Classified contigs from four *Lunularia cruciata* assemblies (2 female, 2 male) using pairwise alignments and reference genome mapping, with EDTA repeat annotation cross-referencing. Identified high-confidence V chromosome candidate contigs from the male HiFi assembly (M2) and produced a marked FASTA for gene prediction.

## Work Completed

### Scripts created/updated

- `four_assembly_cross/classify_contigs.py` — Fixed repeat flagging: default `--repeat-mapq` changed from 1 → 0 (minimap2 self-alignments are predominantly mapq=0); added `min_aln_len=500` bp filter to exclude noise. Now correctly identifies repeat-flagged contigs in all four assemblies.

- `four_assembly_cross/classify_hifi_contigs.py` — New script. Pairwise classification of the two HiFi assemblies (F2 vs M2) only, using reciprocal PAF coverage. Categories: `autosome`, `sex_specific`, `partial`, `repeat_ambiguous`. Outputs per-contig TSVs and sex_specific contig name lists.

- `four_assembly_cross/classify_by_reference.py` — New script. Classifies contigs from all four assemblies by mapping to the chromosome-scale reference (GCA_948567375.1), with EDTA repeat annotation cross-referencing. For each alignment, computes overlap of the *target* (reference) interval with repeat annotations to distinguish genuine hits from repeat-driven ones. Categories: `autosome`, `U_chromosome`, `organellar`, `unplaced`, `no_hit`. Outputs per-contig TSVs with `repeat_class` column (`non_repeat`, `mixed`, `repeat_driven`).

### Output files (in `four_assembly_cross/`)

- `contig_classification/` — Four-assembly pairwise classification results
- `hifi_classification/` — HiFi-only pairwise classification results
- `ref_classification/` — Reference-based classification results (current best):
  - `{F1,F2,M1,M2}_ref_classification.tsv` — per-contig tables
  - `ref_classification_summary.tsv` — summary counts/Mb per category
- `ERR10480608.bp.p_ctg.v_candidates_marked.fasta` — Full M2 HiFi assembly with V candidate contigs marked in header (`ptg*_V_candidate`), for use in gene prediction

### New PAF files used
Reference alignments for all four assemblies:
- `{asm}_mapped_to_GCA_948567375.1_cmLunCruc20.1_genome.fa.paf`

EDTA repeat GFF:
- `../EDTA/GCA_948567375.1_cmLunCruc20.1_genomic_renamed.fasta.mod.EDTA.TEanno_renamed.gff`

## Key Findings

### U chromosome
- **F1 (Linde female):** 99 contigs, 4.25 Mb map to U reference sequences (OX419763.1 + unloc1/2/4) — matches expected ~4.3 Mb ✓
- **F2 (HiFi female):** 2 contigs, 4.91 Mb — the HiFi assembly resolves the U chromosome into just 2 large contigs ✓
- Both F2 U contigs are `repeat_driven` — the U chromosome itself is >80% repeat-annotated (TE-rich), consistent with sex chromosome degeneration
- M2 contigs mapping to U reference (31 contigs, 1.93 Mb) are predominantly `repeat_driven` — shared TEs, not genuine U sequence

### V chromosome candidates
- **M2 no_hit (49 contigs, 2.21 Mb):** contigs with no alignment to the female reference at mapq≥1
- Cross-checked against F1 and F2 assemblies directly:
  - 9 contigs (0.32 Mb) present in female assemblies at high mapq → excluded (autosomal reference gaps)
  - 4 contigs (0.19 Mb) with weak female hit (10–50% coverage in F2) → borderline, possible PAR
  - **36 contigs, 1.70 Mb: confirmed V candidates** (<10% coverage in both F1 and F2)
- M1 (Linde male) no_hit contigs (2,078 contigs, 4.29 Mb) are overwhelmingly tiny (<5 kb, 96% <5 kb); do not add meaningful V sequence beyond what M2 captures

### F2 assembly artefacts
- F2 has 3,397 contigs (201.7 Mb) with no M2 alignment — these are secondary/haplotig contigs and repeat-rich regions, NOT U chromosome
- Reference mapping cleanly resolves this: only 2 contigs (4.91 Mb) map to U; 278 organellar contigs (9.0 Mb), 1,384 unplaced contigs (127.8 Mb), 1,887 no-hit contigs (75.5 Mb) are repeat/haplotig artefacts

### Assembled V chromosome completeness
- Only ~1.70 Mb assembled vs ~9.5 Mb estimated from Hi-C (scaffold_12 in reference) — expected, as V chromosome is highly repetitive and assembles poorly with HiFi

## Output FASTA for gene prediction
`ERR10480608.bp.p_ctg.v_candidates_marked.fasta`
- All 3,845 M2 contigs present
- 36 V candidate contigs have `_V_candidate` appended to sequence ID (no whitespace)
- e.g. `>ptg000324l_V_candidate`

## Next Steps
- [ ] Run gene prediction on `ERR10480608.bp.p_ctg.v_candidates_marked.fasta`
- [ ] Filter gene predictions to V_candidate contigs for V chromosome gene content analysis
- [ ] Compare V chromosome gene content with U chromosome genes (F2 U contigs)
- [ ] Consider whether the 4 weak_female_hit M2 contigs (ptg002620l, ptg000498l, ptg003237l, ptg001856l) represent PAR sequence
- [ ] Update CLAUDE.md with V chromosome assembly findings

## Tags
`#liverwort` `#Lunularia` `#sex-chromosomes` `#genome-assembly` `#V-chromosome` `#contig-classification`
