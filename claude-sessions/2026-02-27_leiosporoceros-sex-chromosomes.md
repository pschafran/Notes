# Claude Code Session - 2026-02-27

**Project:** Leiosporoceros dussii sex chromosome genomics
**Female genome:** `/media/data/projects/hornwort_sex_chromosomes/analysis/Leiosporoceros/LeiosporocerosJC2/final_genome_prep/`
**Male genome:** `/media/data/projects/hornwort_sex_chromosomes/analysis/Leiosporoceros/LeiosporocerosH23/final_genome_prep/`

## Summary

Completed the repeat + gene landscape analysis for the male genome (LedusM) and generated a five-panel male vs. female comparison figure with emphasis on sex chromosome differentiation. The key finding is that the female V chromosome (S3) and male U chromosome (S5) have dramatically different repeat landscapes despite both being gene-poor, consistent with different stages of sex chromosome degeneration.

## Work Completed

### This session (continuation)

The previous session had completed all female figures and written the male analysis scripts (`plot_repeats.py`, `plot_genome_landscape.py` in the male directory) but had not yet run them or written the comparison figure.

### Scripts run

- `LeiosporocerosH23/final_genome_prep/plot_repeats.py` — generated `LedusM_repeat_composition.pdf/png`
- `LeiosporocerosH23/final_genome_prep/plot_genome_landscape.py` — generated `LedusM_genome_landscape_tracks.pdf/png` and `LedusM_genome_landscape_summary.pdf/png`
- `LeiosporocerosJC2/final_genome_prep/plot_sex_chromosome_comparison.py` — generated `LedusF_LedusM_sex_chromosome_comparison.pdf/png`

### Files Created (this session)
- `LeiosporocerosJC2/final_genome_prep/plot_sex_chromosome_comparison.py` — new multi-genome comparison script
- `LeiosporocerosH23/final_genome_prep/LedusM_repeat_composition.caption.txt`
- `LeiosporocerosH23/final_genome_prep/LedusM_genome_landscape_tracks.caption.txt`
- `LeiosporocerosH23/final_genome_prep/LedusM_genome_landscape_summary.caption.txt`
- `LeiosporocerosJC2/final_genome_prep/LedusF_LedusM_sex_chromosome_comparison.caption.txt`

### Files Created (previous session, now confirmed working)
- `LeiosporocerosH23/final_genome_prep/plot_repeats.py`
- `LeiosporocerosH23/final_genome_prep/plot_genome_landscape.py`

## Key Findings

### Female genome (LedusF / JC2)
- 140.4 Mb, 5 scaffolds + contigs
- **S3 = V chromosome** (22.0 Mb): 98.8% repeat, 252 genes (11.4/Mb), dominated by Unknown satellite TE_00000009 (~380 bp, 35,868 copies, ~80% of S3)
- Autosomes: 49–56% repeat, 113–137 genes/Mb, LTR-dominated

### Male genome (LedusM / H23)
- 121.0 Mb, 5 scaffolds + contigs
- **S5 = U chromosome** (5.3 Mb): 86.5% repeat, 234 genes (44.2/Mb), LTR-dominated (67.4% LTR)
- Autosomes: 49–56% repeat, 113.6–125.1 genes/Mb, LTR-dominated
- New TE class vs. female: `TIR/Tc1_Mariner` (mapped to DNA major class)

### Sex chromosome comparison
| Feature | V chromosome (F-S3) | U chromosome (M-S5) |
|---------|---------------------|---------------------|
| Size | 22.0 Mb | 5.3 Mb |
| Repeat fraction | 98.8% | 86.5% |
| Gene density | 11.4/Mb | 44.2/Mb |
| Dominant repeat | Unknown satellite (~77%) | LTR retrotransposons (~67%) |
| Degeneration stage | Advanced (near gene desert) | Moderate (partial gene retention) |

The V chromosome shows a more advanced degeneration characterized by satellite repeat amplification and near-complete gene loss. The U chromosome is smaller, LTR-rich, and still retains ~4× the V gene density — consistent with UV sex determination theory where both chromosomes independently accumulate repeats and lose genes, but at different rates.

## Figure Inventory

### Female genome (LeiosporocerosJC2/final_genome_prep/)
| File | Description |
|------|-------------|
| `LedusF_repeat_composition.pdf/png` | 4-panel: pie, stacked bar per scaffold, heatmap, area tracks |
| `LedusF_genome_landscape_tracks.pdf/png` | 5-scaffold stacked area landscape (genes+unannotated+repeats) |
| `LedusF_genome_landscape_summary.pdf/png` | 2×2: scatter, violin, gene density bar, compartment bar |

### Male genome (LeiosporocerosH23/final_genome_prep/)
| File | Description |
|------|-------------|
| `LedusM_repeat_composition.pdf/png` | 4-panel: pie, stacked bar per scaffold, heatmap, area tracks |
| `LedusM_genome_landscape_tracks.pdf/png` | 5-scaffold stacked area landscape |
| `LedusM_genome_landscape_summary.pdf/png` | 2×2 summary |

### Comparison figure (LeiosporocerosJC2/final_genome_prep/)
| File | Description |
|------|-------------|
| `LedusF_LedusM_sex_chromosome_comparison.pdf/png` | 5-panel: gene density, repeat fraction, scatter, repeat class composition, compartment composition |

## Technical Details

### Three-way genome partition (used in all landscape tracks)
Every 200 kb window is partitioned into three non-overlapping, exhaustive fractions:
- `gene_only` = union(genes ∪ repeats) − repeat_total (gene bases not overlapping repeats)
- `unannotated` = 1 − union(genes ∪ repeats)
- `repeats` = per-major-class fractions, scaled proportionally to sum to `rep_total`

### get_major() — male-specific TIR remapping
The male genome contains `TIR/Tc1_Mariner` elements (absent from female). The `get_major()` function maps `TIR` → `DNA` so these are grouped under DNA transposons.

### Figure style
- Liberation Sans font, no top/right spines, muted pastel palette
- `pdf.fonttype=42` for editable text in Illustrator
- Color constants in `figure_style.md` in project memory

### Bug fixed in previous session (relevant for understanding scripts)
A multi-counting bug in the per-scaffold detailed repeat bar (Panel B of `plot_repeats.py`) was fixed by adding `rep_by_seq_cls` (keyed by full class name, e.g. `'LTR/Gypsy'`) in addition to `rep_by_seq_major` (keyed by major class, e.g. `'LTR'`). Using the major-keyed dict for each detailed class caused LTR to be counted 3× (once per LTR subtype).

## Scripts Archive
Scripts saved to: `~/Notes/claude-scripts/2026-02-27_leiosporoceros-sex-chromosomes_scripts/`
- `plot_repeats_female.py`
- `plot_genome_landscape_female.py`
- `plot_repeats_male.py`
- `plot_genome_landscape_male.py`
- `plot_sex_chromosome_comparison.py`

## Next Steps
- [ ] Investigate the Unknown satellite (TE_00000009) on the female V chromosome in more detail — is it present on any other scaffolds? What is its structure?
- [ ] Check if synteny/collinearity can be established between female S3 and male S5 (or other scaffolds) using whole-genome alignment
- [ ] Consider adding a 6th panel to the comparison figure showing scaffold sizes to scale
- [ ] Potentially analyze the gene-rich island on S3 (~14–16 Mb) in more detail — are these ancestral pseudoautosomal region (PAR) genes?
- [ ] Compare repeat landscapes with other bryophyte sex chromosomes (Marchantia, Ceratodon)

## Related Files
- Female memory: `/home/peter/.claude/projects/-media-data-projects-.../memory/MEMORY.md`
- Figure style: `/home/peter/.claude/projects/-media-data-projects-.../memory/figure_style.md`

## Tags
`#genomics` `#sex-chromosomes` `#bryophytes` `#Leiosporoceros` `#repeat-annotation` `#python` `#matplotlib`
