# Claude Code Session - 2026-04-09

**Project:** hornwort_sex_chromosomes
**Location:** `/media/data/projects/hornwort_sex_chromosomes/analysis/moss_genomes` and `analysis/liverwort_genomes`
**Model:** claude-sonnet-4-6

## Summary

Scanned moss and liverwort genome directories for new composition stats data, developed a statistical outlier detection pipeline (Rosner ESD + modified Z-score) to identify putative sex/accessory chromosomes from per-scaffold gene and repeat content, and updated all CLAUDE.md files with statistically confirmed results. Also corrected a data error in the liverwort CLAUDE.md (Riccia values swapped) and incorporated Fu et al. 2025 sexuality data for Riccia species.

## Work Completed

### Files Modified
- `analysis/moss_genomes/CLAUDE.md` — Updated 1kp species note (composition stats now exist for all 18); replaced ad-hoc outlier section with ESD+modZ statistical results including GR/G/R/~ notation and modified Z-scores
- `analysis/liverwort_genomes/CLAUDE.md` — Replaced ad-hoc composition table with ESD+modZ results; corrected Riccia value swap; added Fu et al. 2025 sexuality data
- `analysis/CLAUDE.md` — Added Sexuality column to both moss tables; updated Riccia entries in liverwort table

### Files Created
- `analysis/moss_genomes/run_esd_outlier_test.py` — ESD + modified Z-score outlier detection script for mosses
- `analysis/moss_genomes/esd_outlier_results.tsv` — Per-scaffold statistics for all moss species
- `analysis/moss_genomes/esd_outlier_summary.tsv` — Per-species summary of confident outlier calls
- `analysis/liverwort_genomes/run_esd_outlier_test.py` — Same script adapted for liverworts
- `analysis/liverwort_genomes/esd_outlier_results.tsv` — Per-scaffold statistics for liverwort species
- `analysis/liverwort_genomes/esd_outlier_summary.tsv` — Per-species summary for liverworts

### Key Commands
```bash
# Run outlier analysis
cd analysis/moss_genomes && python3 run_esd_outlier_test.py
cd analysis/liverwort_genomes && python3 run_esd_outlier_test.py
```

## Key Decisions & Rationale

- **Rosner ESD (r=5)**: r=5 chosen to accommodate sex chromosomes fragmented across up to 5 scaffolds in unpolished assemblies
- **Combined ESD + modified Z-score criterion**: ESD alone has two failure modes — masking (multiple co-occurring outliers inflate SD, e.g. Fissidens) and tight-cluster false positives (extreme outlier removal leaves tiny residual SD, e.g. Ceratodon autosomes). Modified Z-score (median/MAD-based) acts as a minimum effect-size gate
- **1 Mb size filter**: excludes small unscaffolded contigs with noisy estimates
- **Scaffolds ≥ 1 Mb only**: Fontinalis antipyretica (fragmented assembly) had 0 sequences above threshold — correctly untestable

## Statistical Method

```python
# ESD: one-sided lower (gene depletion) and upper (repeat enrichment)
# Critical value at step i with current_n points:
#   p = alpha / current_n  (one-sided Bonferroni)
#   df = current_n - 2
#   lambda_i = ((current_n-1) * t_crit) / sqrt((current_n-2 + t_crit^2) * current_n)
# n_outliers = largest i where R_i > lambda_i

# Modified Z-score (Iglewicz & Hoaglin 1993):
#   M_i = 0.6745 * (x_i - median) / MAD
# Confident outlier = ESD significant AND |M_i| > 3.5
```

## Notable Findings

### Mosses
| Species | Confident outliers | Notes |
|---|---|---|
| Fissidens javanicus | scaffold_19 [~GR] suggestive only | 78.5 Mb, 0.9% gene — ESD masked by co-occurring low-gene scaffolds |
| Racopilum cuspidigerum | HiC_scaffold_1 [GR] | 122 Mb, largest scaffold, 1.4% gene |
| Paraleucobryum enerve | HiC_scaffold_12+14-17 | Large scaffold gene-only; small scaffolds GR |
| Tetraphis pellucida | scaffold_10/11/13 | 74-89% repeat — extreme heterochromatin signal |
| Barbula amplexifolia | scaffold_14-18 (all 5) | ~37 Mb fragmented, all GR |
| Physcomitrium patens | none | Monoicous — correct null result |
| Leucobryum bowringii | none | High genome variance prevents detection despite visually low-gene contigs |
| Jungermannia (liverwort) | none confident | Autosomes 70-76% repeat — MAD too large for stat significance |

### Liverworts
- Marchantia polymorpha chrV is **gene-enriched** (39.1% vs ~28% autosome avg) — flagged [R] only; gene-rich V is known biology
- Riccia cavernosa chr5: monoicous species, but chr5 is **homologous to R. fluitans chr5** (U sex chr) — shared degenerate composition from common ancestry
- Riccia fluitans chr5: putative **U sex chromosome** (dioicous, female genome)
- **Data correction**: previous CLAUDE.md had gene%/repeat% values swapped between the two Riccia species

## Challenges & Solutions

**Problem:** Ceratodon GG1 autosomes falsely flagged after U removal (tight SD=0.6%)
**Solution:** Modified Z-score gate eliminates these — autosomes don't clear |Z|>3.5

**Problem:** Jungermannia sex chromosome candidates not statistically detectable
**Solution:** Documented as visually suggestive but statistically non-significant; root cause identified (high autosomal repeat inflates MAD)

## Next Steps
- [ ] Verify dioicy/monoicy for priority 1kp moss species from literature (Fissidens, Racopilum, Paraleucobryum, Barbula, Diphyscium, Pyrrhobryum)
- [ ] Consider whether Jungermannia warrants a different statistical approach (e.g., clustering rather than outlier detection)
- [ ] Reannotation of Syntrichia ruralis planned before adding to OrthoFinder

## Related Files
- `analysis/moss_genomes/CLAUDE.md` — full outlier catalogue with statistics
- `analysis/liverwort_genomes/CLAUDE.md` — liverwort outlier catalogue
- Previous session: `2026-04-08_syntrichia-ruralis-gff-busco.md`

## Tags
`#hornwort-sex-chromosomes` `#moss-genomes` `#liverwort-genomes` `#statistics` `#sex-chromosomes` `#outlier-detection`
