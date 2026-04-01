# Claude Code Session - 2026-04-01

**Project:** odp-liverworts
**Location:** `/media/data/projects/hornwort_sex_chromosomes/analysis/synteny/`
**Duration:** ~2026-04-01 afternoon - 18:08 EDT

## Summary
Set up a new Oxford Dot Plots (ODP) analysis directory (`odp_liverworts_20260401/`) for 10 high-quality liverwort genomes, excluding the fragmented Lunularia cruciata and MaSuRCA-assembled Marchantia inflexa. Ran input file validation, then prepared and successfully ran a riparian plot config for all 10 species in phylogenetic order.

## Work Completed

### Files Created
- `odp_liverworts_20260401/config.yaml` — main ODP config for 10 liverwort species
- `odp_liverworts_20260401/chrom/Acrsan.chrom` — generated from primary GFF3 (AUGUSTUS)
- `odp_liverworts_20260401/chrom/Blpus.chrom` — generated from primary GFF3 (AUGUSTUS)
- `odp_liverworts_20260401/chrom/Concon.chrom` — generated from primary GFF3 (AUGUSTUS)
- `odp_liverworts_20260401/chrom/Junere.chrom` — generated from primary GFF3 (AUGUSTUS)
- `odp_liverworts_20260401/chrom/Radobc.chrom` — generated from primary GFF3 (Helixer)
- `odp_liverworts_20260401/odp/step3-riparian/liverworts_all.yaml` — riparian plot config
- `odp_liverworts_20260401/odp/step3-riparian/liverworts_all/config.yaml` — symlink to above
- `~/Notes/claude-scripts/2026-04-01_odp-liverworts-config-validation_scripts/check_odp_files_liverworts.sh` — validation script for liverwort species

### Files Symlinked
- `odp_liverworts_20260401/proteins/{Mapol,Maqua,Ricav,Riflu,Rinat}.fa` → odp_bryophytes_20260318 (already validated)
- `odp_liverworts_20260401/chrom/{Mapol,Maqua,Ricav,Riflu,Rinat}.chrom` → odp_bryophytes_20260318
- `odp_liverworts_20260401/genomes/*.fa` → all symlinked to source genomes in liverwort_genomes/
- `odp_liverworts_20260401/proteins/{Acrsan,Blpus,Concon,Junere,Radobc}.fa` → converted to real copies before validation (to avoid modifying originals)

### Commands Executed
```bash
# Chrom file generation for AUGUSTUS-annotated species (Acrsan, Blpus, Concon, Junere)
awk -F'\t' '$3=="mRNA" {
    split($9, a, ";"); gsub("ID=", "", a[1]);
    print a[1]"\t"$1"\t"$7"\t"$4"\t"$5
}' input.primary.gff3 > output.chrom

# Radula (Helixer GFF) had a bug: some mRNA ID= attributes had comma-delimited junk appended
# e.g., ID=Robc_5-13-25.v6.Chromosome01.g000030.t1,eggnog_free_text_desc;Parent=...
# Fix: gsub(",.*", "", a[1]) after splitting on ";"

# Riparian run
cd odp_liverworts_20260401/odp/step3-riparian/liverworts_all
conda run -n oxforddotplots snakemake --snakefile ~/bin/odp/scripts/odp_rbh_to_ribbon -j 4
```

## Species Included (10 total)

| Code | Species | Proteins | BUSCO | Assembly type |
|------|---------|----------|-------|---------------|
| Mapol | Marchantia polymorpha | 18,007 | — | chromosome-level |
| Maqua | Marchantia quadrata | 16,218 | — | chromosome-level |
| Ricav | Riccia cavernosa | 20,184 | — | chromosome-level |
| Riflu | Riccia fluitans | 20,307 | — | chromosome-level |
| Rinat | Ricciocarpos natans | 15,583 | — | chromosome-level |
| Concon | Conocephalum conicum | 17,677 | 96.0% | scaffold-level |
| Blpus | Blasia pusilla | 18,041 | 94.6% | scaffold-level |
| Junere | Jungermannia erecta | 30,742 | 91.5% | HiC scaffold-level |
| Radobc | Radula obconica | 24,917 | 83.3% (odb10) | chromosome-level |
| Acrsan | Acrolejeunea sandvicensis | 22,677 | 93.7% | scaffold-level |

## Key Decisions & Rationale
- **Symlinks vs copies for existing species**: Mapol/Maqua/Ricav/Riflu/Rinat already validated in odp_bryophytes_20260318 — symlinked directly to avoid duplication.
- **Real copies for new species proteins**: New species' protein files were symlinks to originals; converted to copies before running the IFS-fixing validation script to prevent modifying source files.
- **Radula GFF bug**: Some Helixer mRNA features had malformed attributes (`ID=gene.t1,eggnog_free_text_desc;Parent=...` with comma instead of semicolon before the eggnog tag). Fixed by stripping everything after the first comma in the ID value.
- **Junere scaffolds**: Excluded `HiC_scaffold_796` and `HiC_scaffold_797` from chromorder (2 Mb each, oddly numbered — likely assembly artifacts).
- **Validation script**: Existing `check_odp_files.sh` has a hardcoded bryophytes species list, so a liverwort-specific version was written to `~/Notes/claude-scripts/2026-04-01_odp-liverworts-config-validation_scripts/`. All 10 species passed with no IFS, no name mismatches, no chromosome ID issues.

## Phylogenetic Order (riparian plot, top = earliest-diverging)
Based on NCBI taxdump (`/media/data/resources/ncbi_taxdump/rankedlineage.dmp`):
```
Marchantiopsida:
  Marchantiales (Marchantiaceae):  Mapol → Maqua
  Marchantiales (Conocephalaceae): Concon
  Marchantiales (Ricciaceae):      Ricav → Riflu → Rinat
  Blasiales (Blasiaceae):          Blpus
Jungermanniopsida:
  Jungermanniales (Jungermanniaceae): Junere
  Porellales (Radulaceae):            Radobc
  Porellales (Lejeuneaceae):          Acrsan
```

## Technical Details
- Chromorder for scaffold-level assemblies: scaffolds ≥ 1 Mb only. ODP reordered these by size when the riparian ran (the yaml was modified in-place by ODP).
- RBH source: `odp_liverworts_20260401/odp/step2-figures/synteny_nocolor/` (*.plotted.rbh files)
- Validation script location: `~/Notes/claude-scripts/2026-04-01_odp-liverworts-config-validation_scripts/check_odp_files_liverworts.sh`

## Next Steps
- [ ] ODP pairwise run still completing remaining pairs (Ricav/Riflu/Rinat vs each other and vs Mapol/Maqua) — these are all covered in the riparian config
- [ ] Riparian plot output is `odp_liverworts_20260401/odp/step3-riparian/liverworts_all/output.pdf`
- [ ] Consider adding Lunularia cruciata (male + female HiFi) once assembly is less fragmented
- [ ] Consider adding Marchantia inflexa once a better assembly is available

## Related Files
- Previous session: `~/Notes/claude-sessions/2026-03-19_odp-bryophytes-config-validation.md`
- Bryophytes ODP run: `/media/data/projects/hornwort_sex_chromosomes/analysis/synteny/odp_bryophytes_20260318/`

## Tags
`#synteny` `#ODP` `#liverworts` `#riparian` `#comparative-genomics`
