# Claude Code Session - 2026-04-17

**Project:** hornwort_sex_chromosomes — Ks analysis
**Location:** `/media/data/projects/hornwort_sex_chromosomes/analysis/ks`
**Duration:** ~21:00 - 21:12 EDT

## Summary
Set up and ran a two-part Ks analysis using wgd v2 to contextualise U/V sex chromosome divergence in hornworts. Completed all `wgd dmd` and gametolog `wgd ksd` runs; the interspecific `wgd ksd` run is still in progress at session end.

## Work Completed

### Files Created
- `analysis/ks/CLAUDE.md` — full pipeline documentation for this analysis
- `analysis/ks/filter_wgd_families.py` — script to filter wgd dmd TSVs (two modes: gametologs, interspecific)
- `analysis/ks/Ks_tree_wgd.newick` — copy of species tree with `.fa` appended to tip names to match wgd TSV column headers
- `analysis/ks/gametologs_Papro.tsv` — filtered paranome for Papro U/V pairs (68 families)
- `analysis/ks/gametologs_Phchi.tsv` — filtered paranome for Phchi U/V pairs (33 families)
- `analysis/ks/gametologs_Phphy.tsv` — filtered paranome for Phphy U/V pairs (103 families)
- `analysis/ks/gametologs_Ledus.tsv` — filtered paranome for Ledus U/V pairs (21 families)
- `analysis/ks/gametologs_AnangRef2.tsv` — filtered paranome for Anang U/V pairs (28 families)
- `analysis/ks/interspecific_autosomal.tsv` — filtered MRBH orthologs (7290 families; 18 Ledus sex-linked excluded)

### Files Modified
- `analysis/ks/CLAUDE.md` — updated to use simplified filenames (Ledus.fa etc.) and added Noaen
- `analysis/synteny/CLAUDE.md` — added wgd v2 reference section

### Directories Created
- `wgd_dmd_Papro/`, `wgd_dmd_Phchi/`, `wgd_dmd_Phphy/`, `wgd_dmd_Ledus/`, `wgd_dmd_AnangRef2/` — paranome outputs
- `wgd_dmd_interspecific/` — MRBH ortholog output (Ledus focal)
- `wgd_ksd_Papro/`, `wgd_ksd_Phchi/`, `wgd_ksd_Phphy/`, `wgd_ksd_Ledus/`, `wgd_ksd_AnangRef2/` — gametolog Ks outputs (all completed)

### Commands Executed
```bash
# dmd runs (all completed)
conda run -n wgd2 wgd dmd {Species}.fa -o wgd_dmd_{Species} -n 4  # x5 species
conda run -n wgd2 wgd dmd Ledus.fa Papro.fa Phchi.fa Phphy.fa AnangRef2.fa Noaen.fa \
    --focus Ledus.fa -o wgd_dmd_interspecific -n 8

# filtering
python filter_wgd_families.py gametologs wgd_dmd_{Species}/{Species}.fa.tsv \
    {Species}.U {Species}.V gametologs_{Species}.tsv   # x5 species
python filter_wgd_families.py interspecific wgd_dmd_interspecific/merge_focus.tsv \
    Ledus Ledus.U Ledus.V -o interspecific_autosomal.tsv

# ksd runs (gametologs all completed; interspecific still running)
conda run -n wgd2 wgd ksd gametologs_{Species}.tsv {Species}.fa -o wgd_ksd_{Species}  # x5
conda run -n wgd2 wgd ksd interspecific_autosomal.tsv \
    Ledus.fa Papro.fa Phchi.fa Phphy.fa AnangRef2.fa Noaen.fa \
    --spair Ledus.fa Papro.fa --spair Ledus.fa Phchi.fa \
    --spair Ledus.fa Phphy.fa --spair Ledus.fa AnangRef2.fa --spair Ledus.fa Noaen.fa \
    --speciestree Ks_tree_wgd.newick -o wgd_ksd_interspecific
```

## Key Decisions & Rationale
- **Merged U+V genomes as input:** Each merged genome has one copy of autosomes + one U + one V chromosome, making gametologs appear as within-genome paralogs identifiable by scaffold name
- **Sex chromosome scaffold IDs confirmed from CDS files:** Papro/Phchi/Phphy/Ledus use `{Sp}.U` and `{Sp}.V`; AnangRef2 uses `chrU_` and `chrV_`
- **filter_wgd_families.py strips non-UV genes from gametolog families:** Ensures wgd ksd only computes U-V Ks, not autosomal-UV pairs
- **18 Ledus sex-linked genes excluded from interspecific calibration:** Prevents gametologs from inflating the autosomal divergence distribution
- **Ks_tree_wgd.newick:** wgd requires tree tip names to match TSV column headers exactly; wgd dmd uses filenames as species names, so tips needed `.fa` suffix

## Bug Fixed
`filter_wgd_families.py` interspecific mode had an off-by-one column index: the MRBH TSV header starts with species names at column 0, but data rows have GF ID at column 0 (species data offset by 1). Fixed by using `header_idx + 1` for the data column index.

## Current Status
- All 5 gametolog `wgd ksd` runs: **complete**
- Interspecific `wgd ksd` run: **still running** (background process, will survive session close)
- Monitor: `tail -f analysis/ks/wgd_ksd_interspecific.log`

## Next Steps
- [ ] Check interspecific ksd output once complete (`wgd_ksd_interspecific/`)
- [ ] Examine output format: `*.ks.tsv` columns — key ones are `dS` (Ks), `gene1`, `gene2`
- [ ] For gametolog outputs, may want to filter pairs to keep only U-V (not U-U or V-V) before plotting
- [ ] Build visualisation overlaying gametolog Ks distributions vs Ledus interspecific Ks distribution
- [ ] Interpret: if U/V Ks < Ledus interspecific Ks → sex chromosomes younger than Leiosporoceros divergence

## Input Files
| File | Species |
|---|---|
| `Ledus.fa` | *Leiosporoceros dussii* |
| `Papro.fa` | *Paraphymatoceros proskaueri* |
| `Phchi.fa` | *Phaeomegaceros chiloensis* |
| `Phphy.fa` | *Phymatoceros phymatodes* |
| `AnangRef2.fa` | *Anthoceros angustus* |
| `Noaen.fa` | *Nothoceros aenigmaticus* (interspecific only) |

## Tags
`#hornworts` `#ks-analysis` `#sex-chromosomes` `#wgd` `#python`
