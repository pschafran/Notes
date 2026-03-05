# Claude Code Session - 2026-03-05

**Project:** hornwort-sex-chromosomes
**Location:** `/media/data/projects/hornwort_sex_chromosomes/analysis/`
**Species covered:** *Leiosporoceros dussii* (LedusF/LedusM), *Phymatoceros phymatodes* (PhphyF/PhphyM)

## Summary

Produced comparative methylation analyses (Miami Manhattan plot, correlation scatter, per-site coverage/context figure) for Leiosporoceros dussii and Phymatoceros phymatodes, following the same pipeline established for Paraphymatoceros in the previous session. Discovered and fixed a GFF scaffold naming mismatch in the Phymatoceros female genome (EDTA TEanno.gff3 used contig names; correct file is `PhphyF_repeat_annotations.gff` which contains scaffold-level coordinates).

---

## Work Completed

### Files Created
- `Leiosporoceros/sex_chromosome_analyses/plot_methylation_manhattan.py`
  - Full pipeline: WGA (minimap2 asm5), syntenic 200 kb windows, Miami Manhattan, correlation scatter, context figure
  - Three output figures: `LedusF_LedusM_methylation_manhattan.pdf/png`, `_correlation.pdf/png`, `_context.pdf/png`

- `Phymatoceros_phymatodes/sex_chromosome_analyses/plot_methylation_manhattan.py`
  - Same pipeline; handles mixed PacBio (female, 18-col jasmine bed) and ONT (male, standard modkit bedMethyl)
  - Three output figures: `PhphyF_PhphyM_methylation_manhattan.pdf/png`, `_correlation.pdf/png`, `_context.pdf/png`

### Files Modified
- `Phymatoceros_phymatodes/sex_chromosome_analyses/plot_methylation_manhattan.py`
  - `GFF_F`: `pilon.3_renamed.fasta.mod.EDTA.TEanno.gff3` → `PhphyF_repeat_annotations.gff`
  - Fix: female context had showed 0% in all repeat categories due to contig/scaffold name mismatch

---

## Key Results

### Leiosporoceros dussii methylation (200 kb windows, 555 syntenic autosomal windows)
| | Female | Male |
|---|---|---|
| Autosomal 5mC | 56.3% | **67.9%** |
| U/V chromosome 5mC | 78.5% (U) | 78.4% (V) |
| Pearson r (syntenic) | 0.499 | — |

- **Inverted pattern**: male has higher autosomal methylation than female (opposite of Paraphymatoceros)
- Male intermediate sites have median coverage 12× vs unmethylated 21× — some low-coverage noise contribution possible
- Sex chromosomes (U and V) are equally highly methylated despite very different sizes (U=22 Mb, V=5.3 Mb)

### Phymatoceros phymatodes methylation (200 kb windows, 814 syntenic autosomal windows)
| | Female (PacBio) | Male (ONT) |
|---|---|---|
| Autosomal 5mC | 52.6% | 52.8% |
| U/V chromosome 5mC | 80.1% (U) | 78.0% (V) |
| Pearson r (syntenic) | 0.808 | — |

- **No sex difference** in autosomal methylation — both sexes ~52.6–52.8%
- Female coverage ~50× (PacBio); male coverage ~10× (ONT) — methodological difference noted in figure
- Both sex chromosomes highly methylated; U chr (4.9 Mb) vs V chr (4.0 Mb) — both small
- Context (post-fix): female and male both show repeat-enriched methylated sites (46.7% and 51.6% in-repeat for methylated CpGs)

### Cross-species autosomal methylation summary
| Species | Female | Male | Sex difference |
|---|---|---|---|
| Paraphymatoceros proskaueri | 59.3% | 36.7% | F >> M |
| Leiosporoceros dussii | 56.3% | 67.9% | M > F |
| Phymatoceros phymatodes | 52.6% | 52.8% | None |

---

## Key Decisions & Rationale

- **Use `*_repeat_annotations.gff` not EDTA TEanno.gff3**: The EDTA annotation may have been run on a pre-scaffolding assembly with contig names (C*), while the final genome uses scaffold names (S*). The project-level `*_repeat_annotations.gff` (which the existing landscape scripts also use) contains scaffold-level entries. Always verify: `cut -f1 file.gff | grep -v '^#' | sort -u | grep '\.S'`

- **PacBio 18-col bed format compatible with modkit**: jasmine-derived beds have 18 columns but columns 0–10 match modkit format exactly. Filter `p[3] == 'm'` works for both; just note 5mC-only in legend.

- **F_AUTO for Leiosporoceros must exclude S3** (U chromosome): `['LedusF.S1','LedusF.S2','LedusF.S4','LedusF.S5']` — the U chromosome is S3 in this species, not S5 as in others.

---

## Technical Details

### Methylation data formats
- **ONT/modkit**: standard bedMethyl, 11+ columns; col 3 = mod_type (`m`=5mC, `h`=5hmC); col 9 = Nvalid_cov; col 10 = pct_modified
- **PacBio/jasmine**: 18-column variant; columns 0–10 identical to modkit. Example:
  `PhphyF.S1  606  607  m  1  +  606  607  255,0,0  1  100.00  1  0  0  0  0  0  10`
  No 5hmC data (jasmine only writes `MM:Z:C+m?` tag, not `C+h`)

### GFF repeat annotation pitfall
The `pilon.3_renamed.fasta` is the polished pre-scaffolding assembly; EDTA was run on it using contig names. When contigs were scaffolded and the final genome renamed S1-S5, the EDTA GFF3 was not updated. The `*_repeat_annotations.gff` file was regenerated on the final scaffold assembly and has correct names. Check:
```bash
cut -f1 file.gff | grep -v '^#' | sort -u | grep '\.S'
```

### Per-site context analysis
Bisect-based interval labeling (from `papro_methyl_context.py` pattern):
- Load gene intervals from GTF (`gene` feature)
- Load repeat intervals from GFF3 (`Classification=` attribute), grouped by major class
- Use `bisect.bisect_left/right` to label each of N sites into context categories
- Categories: intergenic / gene_only / gene+repeat / repeat_only (broken down by major class)

### scipy unavailable
`from scipy import stats` raises `TypeError` in `_fitpack_impl.py`. All statistics (pearsonr, spearmanr, linregress) implemented with pure NumPy + `math.erfc` for p-values.

---

## Challenges & Solutions

**Problem:** Phymatoceros female context showed 0% in all repeat categories.
**Solution:** The EDTA GFF3 used contig names (PhphyF.C1...) from the pre-scaffolding assembly. Switched `GFF_F` to `PhphyF_repeat_annotations.gff`, which has scaffold-level (S1-S5) coordinates. Female repeat context now shows expected 14.5% → 46.7% repeat fraction from unmethylated → methylated sites.

---

## Next Steps
- [ ] Investigate whether male intermediate LTR sites in Leiosporoceros are noise (low coverage 12×) or biological (hemi-methylation) — compare coverage distribution vs unmethylated
- [ ] Check whether Leiosporoceros inverted F/M methylation difference holds when restricting to high-coverage sites
- [ ] Produce same analysis for Phaeomegaceros fimbriatus if methylation data exists
- [ ] Consider producing a combined cross-species methylation summary figure

---

## Related Files
- Previous session: `2026-03-04_hornwort-sex-chromosomes.md`
- Leiosporoceros script: `Leiosporoceros/sex_chromosome_analyses/plot_methylation_manhattan.py`
- Phymatoceros script: `Phymatoceros_phymatodes/sex_chromosome_analyses/plot_methylation_manhattan.py`

## Tags
`#hornwort` `#sex-chromosomes` `#methylation` `#transposable-elements` `#genome-landscape` `#python`
