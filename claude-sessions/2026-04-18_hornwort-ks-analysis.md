# Claude Code Session - 2026-04-18

**Project:** hornwort_sex_chromosomes — Ks analysis (continued)
**Location:** `/media/data/projects/hornwort_sex_chromosomes/analysis/ks`
**Duration:** short session

## Summary
Fixed TE gene mapping in `plot_ks_vs_position.py` (scaffold name substitution was missing), regenerated the Ks vs. chromosomal position figure with correctly coloured TE-related gametologs, and discussed the biological interpretation of high-Ks TE-flagged pairs.

## Work Completed

### Files Modified
- `analysis/ks/plot_ks_vs_position.py` — fixed TE gene mapping (see below); added `scaf_f`/`scaf_m` scaffold rename tuples to SPECIES config

## Key Bug Fixed

`load_te_genes()` was doing `query.replace(genotype_prefix, merged_prefix)` (e.g. `PaproF` → `Papro`) but not renaming the scaffold component, so `PaproF.S5G000300` became `Papro.S5G000300` instead of the correct `Papro.UG000300`. Fix: build a fragment replacement string `{prefix}.{old_scaf}G` → `{merged}.{new_scaf}G` before falling back to prefix-only substitution for autosomal genes.

Scaffold mappings (from `analyze_sex_chr_genes.py`):
| Species | Female sex chr | Male sex chr |
|---------|---------------|--------------|
| Papro | `PaproF.S5` → `Papro.U` | `PaproM.S5` → `Papro.V` |
| Phchi | `PhchiF.S4` → `Phchi.U` | `PhchiM.S6` → `Phchi.V` |
| Phphy | `PhphyF.S5` → `Phphy.U` | `PhphyM.S5` → `Phphy.V` |
| Ledus | `LedusF.S3` → `Ledus.U` | `LedusM.S5` → `Ledus.V` |

## TE-flagged gametolog pair counts (after fix)
| Species | TE pairs | Total pairs | % |
|---------|---------|-------------|---|
| Papro   | 55      | 117         | 47% |
| Phchi   | 0       | 36          | 0% |
| Phphy   | 8       | 559         | 1% |
| Ledus   | 97      | 114         | 85% |

Phchi 0/36 is a genuine result — the 28 unique Phchi U genes in the Ks output don't overlap with the 8 TE-annotated U genes (the TE genes lack V-chromosome counterparts that passed Ks filtering).

## Key Discussion Points

**Why TE-related genes could show elevated Ks:** No strong biological reason — synonymous sites are neutral. High Ks in TE pairs is most likely methodological:
1. Alignment artifacts (repetitive/divergent sequences → poor MAFFT/PAL2NAL alignment → inflated Ks)
2. Incorrect pairing (multiple TE family members grouped by MCL; non-orthologous copies paired)
3. Saturation (very old TE-derived sequences; multiple hits at same synonymous sites)

**Two compact TE clusters on Papro V (~2-2.5 Mb and ~3.5-4 Mb, Ks > 1):** User suggested possible ancient tandem duplication. More parsimonious alternatives:
- Pericentromeric TE clustering (two flanks of an unassembled centromere)
- Two TE insertion hotspots with alignment-artifact Ks values
Testable: check if the same U gene partners appear in both V clusters (tandem dup prediction) vs. different U genes per cluster. Synteny analysis would help distinguish.

## Next Steps
- [ ] Optionally: replot Ks distributions excluding TE-flagged pairs to see if the high-Ks tail disappears
- [ ] Optionally: identify specific gene pairs in the two Papro V TE clusters and check for shared U partners
- [ ] Consider synteny analysis to test tandem duplication hypothesis on Papro V

## Tags
`#hornworts` `#ks-analysis` `#sex-chromosomes` `#transposable-elements` `#python`
