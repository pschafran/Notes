# Claude Code Session - 2026-04-02

**Project:** odp-mosses-setup
**Location:** `/media/data/projects/hornwort_sex_chromosomes/analysis/synteny/`
**Duration:** 2026-04-02

## Summary
Updated the liverworts riparian config to include Lunularia cruciata (which had been added manually to the ODP run), corrected its phylogenetic placement using the NCBI taxdump, and ordered its chromosomes against Rinat using RBH hit counts. Then set up a new ODP analysis directory for 26 moss species (`odp_mosses_20260402/`) and updated both relevant CLAUDE.md files.

## Work Completed

### Files Modified
- `odp_liverworts_20260401/odp/step3-riparian/liverworts_all.yaml` — added Lucru to species_order (after Rinat) and chromorder (RBH-sorted against Rinat)
- `analysis/CLAUDE.md` — split moss species table into fully-processed vs 1kp sections; added 16 new code prefixes; added Fontinalis exclusion note
- `analysis/synteny/CLAUDE.md` — added entries for odp_liverworts_20260401 and odp_mosses_20260402

### Files Created
- `odp_mosses_20260402/config.yaml` — ODP config for 26 moss species
- `odp_mosses_20260402/chrom/*.chrom` — 16 new chrom files (from AUGUSTUS primary.gff3)
- `odp_mosses_20260402/proteins/*.fa` — 16 new protein files (copied from source)
- `odp_mosses_20260402/genomes/*.fa` — 26 genome symlinks (16 new → moss_genomes/, 10 → odp_bryophytes_20260318/)
- `odp_mosses_20260402/proteins/*.fa` and `chrom/*.chrom` — 10 symlinks → odp_bryophytes_20260318/

### Commands Executed
```bash
# Chrom generation for 1kp species (AUGUSTUS format)
awk -F'\t' '$3=="mRNA" {
    split($9, a, ";"); gsub("ID=", "", a[1]);
    print a[1]"\t"$1"\t"$7"\t"$4"\t"$5
}' input.primary.gff3 > output.chrom

# Validation (inline Python) — all 16 new species passed:
# no IFS, no name mismatches, protein counts == chrom gene counts
```

## Key Decisions & Rationale
- **Lucru phylogenetic position:** Lunulariales is sister to Marchantiales (not within it — corrected from initial wrong claim). Confirmed via NCBI taxdump. Placed after Rinat (last Marchantiales species) in riparian config.
- **Lucru chromorder:** Ordered by RBH hit count against Rinat. Best matches: chr1→chr1, chr3→chr2, chr2→chr3, chr6→chr4, chr7→chr5, chr5→chr6, chr4→chr7, chr8→chr8. chr9+unloc placed last (mirrors chrV position in Rinat).
- **10 existing moss species:** Symlinked directly from odp_bryophytes_20260318 (already validated).
- **16 new 1kp species:** Proteins copied (not symlinked) to protect source files during any future validation. Chrom generated from primary.gff3.
- **Fontinalis excluded:** Genome has raw numeric IDs (>1, >3...) that don't match MAKER GFF scaffold IDs (100381, 101507...).
- **Archidium included as-is:** 53,463 proteins (anomalously high vs 18–38k for other mosses), but protein count matches GFF mRNA count and BUSCO shows only modest duplication. Flagged in CLAUDE.md.

## Species Table — new 1kp codes

| Code | Species | BUSCO |
|---|---|---|
| Aesub | Aerobryopsis subdivergens | 96.6% |
| Aralt | Archidium alternifolium | 96.1% |
| Atang | Atrichum angustatum | 93.4% |
| Baamp | Barbula amplexifolia | 94.5% |
| Brnor | Bryoxiphium norvegicum | 95.9% |
| Diful | Diphyscium fulvifolium | 94.3% |
| Fijav | Fissidens javanicus | 95.7% |
| Grobt | Grimmia obtusifolia | 96.8% |
| Hecil | Hedwigia ciliata | 95.5% |
| Lebow | Leucobryum bowringii | 95.6% |
| Paene | Paraleucobryum enerve | 95.3% |
| Phtur | Philonotis turneriana | 95.5% |
| Ptwil | Ptychomitrium wilsonii | 95.6% |
| Pyspi | Pyrrhobryum spiniforme | 96.1% |
| Racus | Racopilum cuspidigerum | 96.6% |
| Tepel | Tetraphis pellucida | 96.6% |

## Next Steps
- [ ] Run ODP snakemake in `odp_mosses_20260402/` (325 pairwise comparisons — will take significantly longer than liverworts)
- [ ] Investigate Archidium alternifolium protein count anomaly
- [ ] Fix Fontinalis antipyretica genome/GFF ID mismatch to enable future inclusion
- [ ] Create riparian plot config for mosses once ODP run completes

## Related Files
- Previous session: `~/Notes/claude-sessions/2026-04-01_odp-liverworts.md`
- Liverworts ODP run: `odp_liverworts_20260401/`
- Moss genomes: `/media/data/projects/hornwort_sex_chromosomes/analysis/moss_genomes/`

## Tags
`#synteny` `#ODP` `#mosses` `#riparian` `#comparative-genomics`
