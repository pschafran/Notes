# Claude Code Session - 2026-04-01

**Project:** hornwort_sex_chromosomes — CLAUDE.md review
**Location:** `/media/data/projects/hornwort_sex_chromosomes/analysis`

## Summary

Reviewed all 11 CLAUDE.md files across the hornwort sex chromosomes project tree. Found and fixed several species naming inconsistencies and missing entries in species reference tables.

## Work Completed

### Files Modified
- `analysis/CLAUDE.md` — Added Noyun (*Notothylas yunnanensis*) and Phlae (*Phaeoceros laevis*) to Other Hornwort Species table; updated section header to note Phlae is dioicous (not hermaphroditic) with sex of genome individual unknown; added Mainf (*Marchantia inflexa*) to Liverwort Species table
- `analysis/liverwort_genomes/CLAUDE.md` — Fixed code prefix for *Marchantia inflexa* from incorrect `Minflexa` to `Mainf (MainfF/MainfM)` matching actual file naming convention
- `analysis/orthofinder/viridiplantae_20260307/genomes_only/CLAUDE.md` — Removed Noyun and Phlae from focal U/V species list; corrected them to "non-focal species"

## Key Decisions & Rationale

- **Noyun / Phlae classification:** Both appear in the viridiplantae OrthoFinder run but are not focal species. Noyun (*Notothylas yunnanensis*) is hermaphroditic; Phlae (*Phaeoceros laevis*) is dioicous but the sex of the genome individual is unknown. Neither was created as part of this project.
- **Marchantia inflexa prefix:** The actual PROT files use `MainfF_` / `MainfM_` prefixes; the `Minflexa` label in the liverwort_genomes CLAUDE.md was incorrect and has been standardised.

## CLAUDE.md Files Reviewed (all 11)

| File | Status |
|---|---|
| `/home/peter/.claude/CLAUDE.md` | No changes needed |
| `hornwort_sex_chromosomes/CLAUDE.md` | No changes needed |
| `analysis/CLAUDE.md` | Updated (species tables) |
| `analysis/orthofinder/CLAUDE.md` | No changes needed |
| `analysis/orthofinder/viridiplantae_20260307/genomes_only/CLAUDE.md` | Updated (species classification) |
| `analysis/synteny/CLAUDE.md` | No changes needed |
| `analysis/liverwort_genomes/CLAUDE.md` | Updated (Mainf prefix) |
| `analysis/liverwort_genomes/Lunularia_cruciata/CLAUDE.md` | No changes needed |
| `analysis/liverwort_genomes/Lunularia_cruciata/female_hifi/final_genome/CLAUDE.md` | No changes needed |
| `analysis/liverwort_genomes/Lunularia_cruciata/male_hifi/final_genome/CLAUDE.md` | No changes needed |
| `analysis/liverwort_genomes/Marchantia_inflexa/CLAUDE.md` | No changes needed |
| `analysis/moss_genomes/CLAUDE.md` | No changes needed |

## Next Steps
- [ ] Resolve the TODO in `analysis/CLAUDE.md`: confirm which moss species are dioicous vs monoicous from the literature
- [ ] Consider whether Phlae sex could be determined from read-depth analysis if data is available

## Tags
`#documentation` `#claude-md` `#species-tables` `#hornworts`
