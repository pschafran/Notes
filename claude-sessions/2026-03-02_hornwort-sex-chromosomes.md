# Claude Code Session - 2026-03-02

**Project:** Hornwort sex chromosome genomics (multi-species)
**Location:** `/media/data/projects/hornwort_sex_chromosomes/analysis/`
**Continuing from:** `2026-02-27_leiosporoceros-sex-chromosomes.md`

## Summary

Three main bodies of work: (1) completed the full analysis pipeline for *Phaeomegaceros chiloensis* (8 scripts, all figures generated); (2) fixed a critical bug in the *Phymatoceros phymatodes* female repeat annotation that caused zero repeat statistics for all scaffolds; (3) performed a detailed analysis of Nanopore CpG methylation data for *Paraphymatoceros proskaueri* and added methylation tracks to the genome landscape and comparison figures.

---

## Part 1 — *Phaeomegaceros chiloensis* full analysis

### Species / accessions
- Female: `Phaeomegaceros_14765-5` → genome prefix `PhchiF`; **U chromosome = PhchiF.S4** (16.3 Mb, 170 genes)
- Male: `Phaeomegaceros_14765-4` → genome prefix `PhchiM`; **V chromosome = PhchiM.S6** (6.2 Mb, 349 genes)
- Working directory: `Phaeomegaceros_fimbriatus/sex_chromosome_analyses/`

### Scripts created and run

| Script | Location | Output |
|--------|----------|--------|
| `plot_repeats.py` (male) | `Phaeomegaceros_14765-4/final_genome_prep/` | `PhchiM_repeat_composition.pdf/png` |
| `plot_repeats.py` (female) | `Phaeomegaceros_14765-5/final_genome_prep/` | `PhchiF_repeat_composition.pdf/png` |
| `plot_genome_landscape.py` (male) | `Phaeomegaceros_14765-4/final_genome_prep/` | `PhchiM_genome_landscape_tracks/summary.pdf/png` |
| `plot_genome_landscape.py` (female) | `Phaeomegaceros_14765-5/final_genome_prep/` | `PhchiF_genome_landscape_tracks/summary.pdf/png` |
| `plot_sex_chromosome_comparison.py` | `sex_chromosome_analyses/` | `PhchiF_PhchiM_sex_chromosome_comparison.pdf/png` |
| `run_RBH_Ks.py` | `sex_chromosome_analyses/` | `rbh_ks/PhchiF_PhchiM_RBH_Ks.tsv` |
| `plot_Ks.py` | `sex_chromosome_analyses/` | `rbh_ks/PhchiF_PhchiM_Ks_distribution.pdf/png` |
| `plot_sex_chromosome_synteny.py` | `sex_chromosome_analyses/` | `rbh_ks/PhchiF_PhchiM_sex_chromosome_synteny.pdf/png` |

### RBH + Ks results
- 13,244 RBH pairs total; **4,538 ok** (valid Ks); 8,702 zero_Ks; 4 short_aln
- Median Ks = **0.0045** (intraspecific, consistent with other hornworts)
- **17 U↔V shared gene pairs** (PhchiF.S4 ↔ PhchiM.S6, flag=ok)
- 9 U-chromosome genes with RBH to male autosomes (no V partner)
- 16 V-chromosome genes with RBH to female autosomes (no U partner)

### Key script parameters (Phchi-specific)
```python
# Male landscape
SCAFFOLDS = ['PhchiM.S1', ..., 'PhchiM.S6']
SEX_CHR   = 'PhchiM.S6'   # V chromosome
SCAF_PALETTE = ['#6aab98']*5 + ['#c47a4a']   # S6 burnt orange
MAJOR_ORDER = ['LTR', 'DNA', 'LINE', 'MITE', 'Unknown']   # no pararetrovirus
get_major: TIR/polinton → DNA, Penelope → LINE

# Female landscape
SCAFFOLDS = ['PhchiF.S1', ..., 'PhchiF.S7']   # 7 scaffolds
SEX_CHR   = 'PhchiF.S4'   # U chromosome, S4 (not S5 or S7!)
SCAF_PALETTE = ['#7396b8']*3 + ['#b87070'] + ['#7396b8']*3

# Synteny plot
LEN_U = 16_316_474   # PhchiF.S4
LEN_V =  6_185_290   # PhchiM.S6
GTF_F = PhchiF_gene_annotations.gtf
GTF_M = PhchiM_gene_annotations.gtf
# Mb ticks: step 2 for U chr, step 1 for V chr
```

---

## Part 2 — *Phymatoceros phymatodes* female repeat GFF bug fix

### Problem
The *Phymatoceros* female (PhphyF) genome showed zero repeat statistics for all scaffolds S1–S5 in all three affected figures. Root cause: the wrong EDTA output file was specified.

| File | Scaffold names | Status |
|------|---------------|--------|
| `pilon.3_renamed.fasta.mod.EDTA.TEanno.gff3` (old) | `PhphyF.C*` (contigs only) | **Wrong** — no S scaffolds present |
| `PhphyF_repeat_annotations.gff` (correct) | `PhphyF.S*` (scaffolds) | **Correct** |

### Additional complication
The correct GFF contains TE classes absent from the original `ORDERED_CLASSES`: `TIR/Tc1_Mariner`, `TIR/Sola2`, `TIR/PiggyBac`, `MITE/DTA`, `Penelope`. `Penelope` (no `/` separator) also failed the original `get_major()` logic and was not grouped under LINE.

### Files fixed

**`Phymato_6/final_genome_prep/plot_repeats.py`** — 5 edits:
1. GFF path → `PhphyF_repeat_annotations.gff`
2. `CLASS_COLORS`: added new TE classes
3. `ORDERED_CLASSES`: added `TIR/Tc1_Mariner`, `TIR/Sola2`, `TIR/PiggyBac`, `Penelope`, `MITE/DTA`
4. `get_major()`: added `if major == 'Penelope': return 'LINE'`
5. `group_map` legend: added TIR classes → DNA transposons, Penelope → LINE

**`Phymato_6/final_genome_prep/plot_genome_landscape.py`** — 2 edits:
1. GFF path → `PhphyF_repeat_annotations.gff`
2. `get_major()`: added `if major == 'Penelope': return 'LINE'`

**`sex_chromosome_analyses/plot_sex_chromosome_comparison.py`** — 2 edits:
1. Female GFF path → `PhphyF_repeat_annotations.gff`
2. `get_major()`: added `if major == 'Penelope': return 'LINE'`

Note: male GFF (`Phymato_ref_genome/pilon.3_renamed.fasta.mod.EDTA.TEanno.gff3`) is correct — it has `PhphyM.S*` names and was not changed.

### Updated `get_major()` pattern (Phymato)
```python
def get_major(cls):
    major = cls.split('/')[0] if '/' in cls else cls
    if major == 'TIR':
        return 'DNA'
    if major == 'Penelope':
        return 'LINE'
    return major
```

All three scripts were re-run successfully after fixes.

---

## Part 3 — *Paraphymatoceros proskaueri* CpG methylation analysis

### Data
| Sample | Sex | Accession | BED file |
|--------|-----|-----------|---------|
| Papro252-3 | Male | PaproM | `Papro252-3/methylation/PaproM_genome.5mCG_5hmCG.bed` |
| Papro252-9.2.2 | Female | PaproF | `Papro252-9.2.2/methylation/PaproF_genome.5mCG_5hmCG.bed` |

Format: 18-column modkit bedMethyl with both `m` (5mC) and `h` (5hmC) rows per CpG position. ~26.6M unique CpG positions per sample (~53M rows total = 2 rows per site).

Sex chromosomes: **PaproF.S5 = U chromosome** (8.35 Mb); **PaproM.S5 = V chromosome** (5.91 Mb).

### Key findings (Nvalid_cov ≥ 5, scaffolds only)

**Global methylation:**
| Modification | Male | Female |
|---|---|---|
| 5mC | 37.5% | 61.0% |
| 5hmC | 1.8% | 3.6% |

**Per-scaffold mean 5mC%:**
| Scaffold | Male | Female | Note |
|---|---|---|---|
| S1 | 36.1% | 58.4% | autosome |
| S2 | 36.4% | 59.3% | autosome |
| S3 | 36.4% | 58.4% | autosome |
| S4 | 39.9% | 64.2% | autosome |
| **S5** | **54.2%** | **86.4%** | **sex chromosome** |

Sex chromosomes are hypermethylated relative to autosomes by ~16–18 pp (V chr) and ~22–28 pp (U chr).

**Distribution shape:** Bimodal — large unmethylated peak (0–9%) and methylated hump. Male methylated peak broad (~40–70%); female peak strongly right-shifted (80–100%), especially on S5.

### Coverage threshold analysis
- Male: stable 37.5–37.8% across Nvalid_cov ≥5 to ≥30 (well covered, ~42× median)
- Female: stable 60–61% for ≥5 to ≥20, then **drops sharply to 50% at ≥30**
- Drop is a sampling artefact: female at ≥30× preferentially retains accessible, less-methylated regions while excluding repeat-dense, highly methylated loci that have lower unique-read coverage from multi-mapping exclusion
- **Nvalid_cov ≥ 5 is the correct threshold** for both samples

On female autosomes, methylation is **inversely correlated with coverage** (78.5% at 20–29×; 20.7% at 50–99×) — high-coverage sites are in uniquely mappable gene-body regions (less methylated); low-coverage sites are in repeat-dense heterochromatin (highly methylated). Pattern is reversed for male (positive correlation, consistent with lower overall methylation and different coverage landscape).

### Processing workflow (from pschafran.github.io documentation)
The project uses **Megalodon v2.5.0** for methylation calling:
```bash
megalodon ../nanopore_reads/ \
  --guppy-config dna_r9.4.1_450bps_hac.cfg \
  --remora-modified-bases dna_r9.4.1_e8 hac 0.0.0 5mc CG 0 \
  --outputs basecalls mappings mods \
  --reference ../genome/genome.fasta
```
**Critical note:** Megalodon output BED is unsorted — must `bedtools sort` before downstream use. The Paraphymatoceros files have `5mCG_5hmCG` in the name, indicating they were generated with a newer pipeline (likely modkit from Dorado-basecalled data) that calls both 5mC and 5hmC simultaneously. The 18-column modkit bedMethyl format is used (col 10 = Nvalid_cov, col 11 = percent_modified).

---

## Part 4 — Methylation added to genome landscape and comparison figures

### Changes to both `plot_genome_landscape.py` scripts (Papro male + female)

**New constants:**
```python
METHYL_BED     = f'{BASE}/../methylation/Papro[MF]_genome.5mCG_5hmCG.bed'
METHYL_MIN_COV = 5
METHYL_COLOR   = '#6a3d7a'   # muted violet for 5mC
HMETH_COLOR    = '#b090c0'   # soft mauve for 5hmC
```

**New functions:**
- `load_methylation(bed_file, scaffolds, min_cov)` — reads bedMethyl, stores per-scaffold sorted numpy arrays of (position, met%) for 5mC and 5hmC separately; returns per-scaffold means
- `win_meth(pos_arr, met_arr, ws, we)` — `bisect`-based O(log n) window mean

**Figure 1 (landscape tracks):** Twin y-axis on each scaffold panel shows 5mC% as a violet line + light fill. Per-scaffold mean 5mC% added to panel title.

**Figure 2 (summary):** Expanded from 2×2 to 3×2 grid. New **Panel E** (spanning full bottom row): per-scaffold 5mC% as bars (scaffold palette) + 5hmC% as dot-line series on secondary right y-axis (0–10%).

### Changes to `plot_sex_chromosome_comparison.py` (Papro)

**New function:**
- `load_methylation_means(bed_file, scaffolds, min_cov)` — lightweight streaming version (no position arrays); returns `{seqid: (mC_mean, hmC_mean)}`
- `auto_meth_avg()` — genome-size-weighted autosomal methylation mean

**Figure:** Height expanded 15→18 inches. New **Panel F** (full-width bottom): 4-group bar chart comparing autosomal mean vs. sex chromosome for both sexes. 5hmC shown as dot-line on secondary right y-axis with per-point value labels.

Panel F values:
| Group | 5mC | 5hmC |
|---|---|---|
| F autosomes | 59.7% | ~3.5% |
| U chr (F-S5) | **86.4%** | 5.1% |
| M autosomes | 37.0% | ~1.8% |
| V chr (M-S5) | **54.2%** | 2.5% |

---

## Files Modified / Created

### Phaeomegaceros chiloensis (new)
- `Phaeomegaceros_14765-4/final_genome_prep/plot_repeats.py`
- `Phaeomegaceros_14765-4/final_genome_prep/plot_genome_landscape.py`
- `Phaeomegaceros_14765-5/final_genome_prep/plot_repeats.py`
- `Phaeomegaceros_14765-5/final_genome_prep/plot_genome_landscape.py`
- `sex_chromosome_analyses/plot_sex_chromosome_comparison.py`
- `sex_chromosome_analyses/run_RBH_Ks.py`
- `sex_chromosome_analyses/plot_Ks.py`
- `sex_chromosome_analyses/plot_sex_chromosome_synteny.py`
- `sex_chromosome_analyses/rbh_ks/PhchiF_PhchiM_RBH_Ks.tsv` (13,244 pairs)
- All output PDF/PNG figures

### Phymatoceros phymatodes (fixed)
- `Phymato_6/final_genome_prep/plot_repeats.py` — GFF path + new TE classes + Penelope
- `Phymato_6/final_genome_prep/plot_genome_landscape.py` — GFF path + Penelope
- `sex_chromosome_analyses/plot_sex_chromosome_comparison.py` — female GFF path + Penelope

### Paraphymatoceros proskaueri (methylation added)
- `Papro252-3/final_genome_prep/plot_genome_landscape.py` — added methylation tracks
- `Papro252-9.2.2/final_genome_prep/plot_genome_landscape.py` — added methylation tracks
- `sex_chromosome_analyses/plot_sex_chromosome_comparison.py` — added panel F

---

## Key Technical Patterns

### Efficient bedMethyl loading with bisect
```python
import bisect

def load_methylation(bed_file, scaffolds, min_cov=5):
    # Returns sorted position arrays per scaffold for fast windowing
    ...

def win_meth(pos_arr, met_arr, ws, we):
    lo = bisect.bisect_left(pos_arr, ws)
    hi = bisect.bisect_right(pos_arr, we)
    if lo >= hi:
        return float('nan')
    return float(np.mean(met_arr[lo:hi]))
```

### Coverage threshold recommendation
Nvalid_cov ≥ 5 is the standard filter (matches author's own awk filter from processing blog). Avoid ≥30 for female Paraphymatoceros — biases estimate downward by ~10 pp.

### Twin-axis methylation in matplotlib
```python
ax2 = ax.twinx()
ax2.plot(pos[valid], mC_win[valid], color=METHYL_COLOR, lw=1.2, zorder=5)
ax2.fill_between(pos[valid], 0, mC_win[valid], color=METHYL_COLOR, alpha=0.08, zorder=4)
ax2.set_ylim(0, 100)
ax2.spines['right'].set_color(METHYL_COLOR)
ax2.spines['top'].set_visible(False)
```

---

## Next Steps
- [ ] Run methylation analysis on other species (Leiosporoceros, Phymatoceros, Phaeomegaceros)
- [ ] Add methylation to Phaeomegaceros chiloensis figures (same pattern as Papro)
- [ ] Compare methylation between sex chromosomes and autosomes across all species
- [ ] Investigate whether methylation correlates with repeat content on a per-window basis
- [ ] Consider gene body methylation analysis (promoter, exon, intron bins) as described in processing blog
- [ ] Verify bedMethyl files are sorted (required for bisect-based windowing)
- [ ] Check if Phaeomegaceros, Leiosporoceros methylation data exists

## Scripts Archive
`~/Notes/claude-scripts/2026-03-02_hornwort-sex-chromosomes_scripts/`
- `phchi_[MF]_plot_repeats.py` — Phaeomegaceros chiloensis repeat scripts
- `phchi_[MF]_plot_genome_landscape.py` — Phaeomegaceros chiloensis landscape scripts
- `phchi_plot_sex_chromosome_comparison.py`
- `phchi_run_RBH_Ks.py`
- `phchi_plot_Ks.py`
- `phchi_plot_sex_chromosome_synteny.py`
- `phphy_F_plot_repeats_fixed.py` — fixed Phymatoceros female repeat script
- `phphy_F_plot_genome_landscape_fixed.py` — fixed Phymatoceros female landscape script
- `phphy_plot_sex_chromosome_comparison_fixed.py` — fixed Phymatoceros comparison script
- `papro_[MF]_plot_genome_landscape_methylation.py` — Paraphymatoceros landscape + methylation
- `papro_plot_sex_chromosome_comparison_methylation.py` — Paraphymatoceros comparison + panel F

## Tags
`#genomics` `#sex-chromosomes` `#bryophytes` `#Phaeomegaceros` `#Phymatoceros` `#Paraphymatoceros` `#methylation` `#bedMethyl` `#modkit` `#python` `#matplotlib`
