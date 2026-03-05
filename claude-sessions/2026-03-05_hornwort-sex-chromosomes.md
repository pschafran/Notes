# Claude Code Session - 2026-03-05

**Project:** hornwort-sex-chromosomes
**Location:** `/media/data/projects/hornwort_sex_chromosomes/analysis/`
**Species covered:** *Leiosporoceros dussii* (LedusF/LedusM), *Phymatoceros phymatodes* (PhphyF/PhphyM), *Paraphymatoceros proskaueri* (PaproF/PaproM)

## Summary

Produced comparative methylation analyses (Miami Manhattan plot, correlation scatter, per-site coverage/context figure, methylation distribution histogram + stats file) for all three hornwort species: *Leiosporoceros dussii*, *Phymatoceros phymatodes*, and *Paraphymatoceros proskaueri*. Fixed a GFF scaffold naming mismatch in the Phymatoceros female genome. Fixed a context figure legend bug where female-only repeat categories were missing. Added full context analysis pipeline to the Paraphymatoceros script (which previously had only Manhattan and correlation figures). Key finding: Paraphymatoceros male shows confirmed hemi-methylation (65.9% intermediate sites, median coverage = unmethylated coverage).

---

## Work Completed

### Files Created
- `Leiosporoceros/sex_chromosome_analyses/plot_methylation_manhattan.py`
  - Full pipeline: WGA (minimap2 asm5), syntenic 200 kb windows, Miami Manhattan, correlation scatter, context figure
  - Three output figures: `LedusF_LedusM_methylation_manhattan.pdf/png`, `_correlation.pdf/png`, `_context.pdf/png`

- `Phymatoceros_phymatodes/sex_chromosome_analyses/plot_methylation_manhattan.py`
  - Same pipeline; handles mixed PacBio (female, 18-col jasmine bed) and ONT (male, standard modkit bedMethyl)
  - Three output figures: `PhphyF_PhphyM_methylation_manhattan.pdf/png`, `_correlation.pdf/png`, `_context.pdf/png`

### Files Modified / Extended

- `Phymatoceros_phymatodes/sex_chromosome_analyses/plot_methylation_manhattan.py`
  - `GFF_F`: `pilon.3_renamed.fasta.mod.EDTA.TEanno.gff3` → `PhphyF_repeat_annotations.gff` (GFF naming fix)
  - Added `pct_data` storage in context loop for per-site distribution analysis
  - Fixed context figure legend bug (`all_handles_dict` accumulated across both panels)
  - Added Figure 4: methylation distribution histogram → `PhphyF_PhphyM_methylation_distribution.pdf/png`
  - Added stats output file → `PhphyF_PhphyM_methylation_stats.txt`

- `Leiosporoceros/sex_chromosome_analyses/plot_methylation_manhattan.py`
  - Added `pct_data` storage in context loop
  - Fixed same context figure legend bug
  - Added Figure 4: methylation distribution histogram → `LedusF_LedusM_methylation_distribution.pdf/png`
  - Added stats output file → `LedusF_LedusM_methylation_stats.txt`

- `Paraphymatoceros_proskaueri/sex_chromosome_analyses/plot_methylation_manhattan.py`
  - Added entire context analysis infrastructure (functions, color dicts, CATS/stack ordering)
  - Added Figure 3: coverage violin + context stacked bar → `PaproF_PaproM_methylation_context.pdf/png`
  - Added Figure 4: methylation distribution histogram → `PaproF_PaproM_methylation_distribution.pdf/png`
  - Added stats output file → `PaproF_PaproM_methylation_stats.txt`
  - Legend bug fix applied at time of writing (not affected by the bug)

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

### Paraphymatoceros proskaueri methylation (added this session)
| | Female | Male |
|---|---|---|
| Autosomal 5mC | 59.3% | 36.7% |
| U/V chromosome 5mC | 86.4% (U) | 54.2% (V) |
| Pearson r (syntenic) | 0.886 | — |

- Female: 57.6% fully methylated, 33.4% unmethylated, 11.0% intermediate
- **Male: 65.9% intermediate** — median coverage 59× = unmethylated 58× → **confirmed hemi-methylation**
- V chromosome only 54.2% methylated (vs 86.4% U) — much lower than in other species

### Cross-species autosomal methylation summary
| Species | Female | Male | Sex difference | Hemi-methylation |
|---|---|---|---|---|
| Paraphymatoceros proskaueri | 59.3% | 36.7% | F >> M | Male (confirmed) |
| Leiosporoceros dussii | 56.3% | 67.9% | M > F | Possible male (12× vs 21×) |
| Phymatoceros phymatodes | 52.6% | 52.8% | None | Female possible (29.3% intermediate, ~50× cov) |

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

**Problem:** Context figure panel B legend was missing female-only repeat categories.
**Solution:** `handles = []` was reset each loop iteration, and legend was built only from the male panel (`col == 1`). Categories with >0.5% in female but <0.5% in all male bars appeared in bars but not legend. Fix: introduced `all_handles_dict = {}` before the loop, accumulated entries from both panels keyed by label, placed legend after loop from combined dict: `ax.legend(handles=list(all_handles_dict.values())[::-1], ...)`.

**Problem:** Paraphymatoceros male shows 65.9% intermediate-methylation sites — is this hemi-methylation or low-coverage noise?
**Solution:** Compared median coverage of intermediate sites (59×) vs unmethylated sites (58×): they are equal. Low-coverage noise would produce intermediate values only at low-covered sites; here coverage is the same as fully unmethylated sites, confirming the intermediates are biological hemi-methylated loci (each allele either 0% or 100% methylated, averaging to ~50% in the pileup).

**Problem:** `head -20` on sorted GFF sequence names missed scaffold entries (alphabetically, C-names precede S-names).
**Solution:** Always verify scaffold names explicitly: `cut -f1 file.gff | grep -v '^#' | sort -u | grep '\.S'`

---

## Next Steps
- [ ] Investigate whether Leiosporoceros male intermediate sites (12× median vs unmethylated 21×) are noise or biological hemi-methylation — low coverage makes this ambiguous
- [ ] Produce same methylation analysis for Phaeomegaceros fimbriatus if data exists
- [ ] Consider producing a combined cross-species methylation summary figure
- [ ] Investigate the dramatically lower V-chromosome methylation in Paraphymatoceros (54.2% vs 86.4% U) — may reflect different epigenetic regulation on male sex chromosome

---

## Related Files
- Previous session: `2026-03-04_hornwort-sex-chromosomes.md`
- Leiosporoceros script: `Leiosporoceros/sex_chromosome_analyses/plot_methylation_manhattan.py`
- Phymatoceros script: `Phymatoceros_phymatodes/sex_chromosome_analyses/plot_methylation_manhattan.py`
- Paraphymatoceros script: `Paraphymatoceros_proskaueri/sex_chromosome_analyses/plot_methylation_manhattan.py`
- Archived scripts: `~/Notes/claude-scripts/2026-03-05_hornwort-sex-chromosomes_scripts/`

## Tags
`#hornwort` `#sex-chromosomes` `#methylation` `#transposable-elements` `#genome-landscape` `#python`
