# Claude Code Session - 2026-04-01

**Project:** Lunularia cruciata female HiFi final genome
**Location:** `/media/data/projects/hornwort_sex_chromosomes/analysis/liverwort_genomes/Lunularia_cruciata/female_hifi/final_genome`

## Summary
Created a CLAUDE.md for the female HiFi final genome directory documenting assembly provenance, file contents, sequence ID conventions, BUSCO stats, and relationships to other Lunularia assemblies. Then produced a BUSCO completeness vs. minimum contig size filter plot mirroring the equivalent analysis already done for the male assembly.

## Work Completed

### Files Created
- `female_hifi/final_genome/CLAUDE.md` — documents assembly provenance chain (hifiasm → rename → EDTA → BRAKER → final prep), file inventory with stats, sequence ID conventions (noting the `LucruF-cmLunCruc20` hyphen format), BUSCO summary, relationship to GCA_948567375.1, and downstream integration pointers
- `female_hifi/final_genome/plot_busco_size_filter.py` — plots BUSCO completeness vs. minimum contig size filter (adapted from male equivalent)
- `female_hifi/final_genome/LucruF_cmLunCruc20_busco_size_filter.pdf` — output figure
- `female_hifi/final_genome/LucruF_cmLunCruc20_busco_size_filter.png` — output figure

### Commands Executed
```bash
python plot_busco_size_filter.py
```

## Key Decisions & Rationale
- **N50 axis in Mb not kb:** The female assembly is far more contiguous than the male (N50 16–23 Mb vs. 266–398 kb for male), so the axis unit was changed to Mb to avoid cluttered tick labels.
- **Size filter range 0–350 kb:** Matched the existing female BUSCO runs (50k, 100k, 150k, 200k, 250k, 300k, 350k) rather than the male range (10k–200k).

## Technical Details

### Female BUSCO data across size filters (viridiplantae_odb12, n=822)

| Min size | Scaffolds | Total (Mb) | N50 (Mb) | C% | S | D | F | M |
|---|---|---|---|---|---|---|---|---|
| 0 (unfiltered) | 3,611 | 778 | 16 | 98.4% | 757 | 52 | 2 | 11 |
| 50 kb | 1,211 | 703 | 20 | 98.2% | 757 | 50 | 3 | 12 |
| 100 kb | 463 | 652 | 23 | 98.1% | 765 | 41 | 3 | 13 |
| 150 kb | 294 | 632 | 23 | 98.1% | 774 | 32 | 3 | 13 |
| 200 kb | 201 | 616 | 23 | 98.1% | 780 | 26 | 3 | 13 |
| 250 kb | 132 | 600 | 23 | 98.1% | 784 | 22 | 3 | 13 |
| 300 kb | 102 | 592 | 23 | 98.1% | 785 | 21 | 3 | 13 |
| 350 kb | 88 | 588 | 23 | 98.1% | 786 | 20 | 3 | 13 |

**Key observation:** BUSCO completeness is essentially flat across all filter thresholds (98.4% → 98.1%), confirming the small contigs in the unfiltered assembly contain negligible gene content. This contrasts sharply with the male assembly where BUSCO drops from 95.1% (unfiltered) to 65.2% at 200 kb. The female assembly is substantially more contiguous.

## Next Steps
- [ ] Decide on a final size-filter cutoff for the female assembly (100 kb appears a reasonable trade-off: 463 scaffolds, 652 Mb, 98.1% BUSCO)
- [ ] Update OrthoFinder input if a size-filtered version is preferred over the full assembly

## Related Files
- Male equivalent: `male_hifi/final_genome/plot_busco_size_filter.py`
- Male figures: `male_hifi/final_genome/LucruM_cmLunCruc9_busco_size_filter.pdf/.png`
- BUSCO run logs: `female_hifi/final_genome/busco_3814397.log`

## Tags
`#Lunularia-cruciata` `#genome-assembly` `#BUSCO` `#python`
