# Claude Code Session - 2026-04-01

**Project:** liverwort_genomes
**Location:** `/media/data/projects/hornwort_sex_chromosomes/analysis/liverwort_genomes`
**Duration:** ~14:00 - 14:35

## Summary
Examined 5 new liverwort genomes added to the project today, created a top-level CLAUDE.md for the liverwort_genomes directory, and launched a sequential EDTA repeat annotation batch run on all 5 species with auto-debugging logic.

## Work Completed

### Files Created
- `analysis/liverwort_genomes/CLAUDE.md` — Top-level documentation for the liverwort_genomes directory (species table, pipeline summaries, GFF3 quirks, EDTA bug fixes, data file patterns, pointers to subdirectory CLAUDE.mds)
- `analysis/liverwort_genomes/run_edta_liverworts.sh` — EDTA batch script (see Scripts section)

## New Species Added 2026-04-01

All 5 originate from CoGe platform (indicated by `scaffold_NNN-SpeciesName` filename convention), annotated with AUGUSTUS. Files per species: genome `.fa`, `.cds.fa`, `.pep.fa`, `.primary.gff3`, `.fa.gz.md5`.

| Species | CoGe prefix | Scaffolds | Genome size | Proteins |
|---|---|---|---|---|
| *Blasia pusilla* | scaffold_890 | 26 | 390 Mb | 18,041 |
| *Plagiochasma appendiculatum* | scaffold_952 | 67 | 429 Mb | 20,623 |
| *Jungermannia erectum* | scaffold_955 | 1,284 | 772 Mb | 30,742 |
| *Conocephalum conicum* | scaffold_664 | 18 | 288 Mb | 17,677 |
| *Acrolejeunea sandvicensis* | scaffold_888 | 40 | 387 Mb | 22,677 |

Scaffold IDs: `scaffold_N` for most; `HiC_scaffold_N` for *Jungermannia*. *Conocephalum* has zero Ns/gaps (gap-free assembly).

### GFF3 structure (important for downstream tools)
- Feature types: **mRNA + CDS only — no `gene` feature lines**
- Gene ID embedded in mRNA attribute: `geneID=jgNNN`
- Protein/CDS IDs: `jgNNN.t1` format — **no species code prefix**
- These will need renaming before use in OrthoFinder

### Genome size calculation bug
Initial genome size estimates were wrong (inflated 2-3×). Root cause: glob pattern `scaffold_*.fa` matched both genome `.fa` and `.pep.fa`, so `wc -c` summed both. Fixed by using `grep -v 'pep\|cds'` to isolate genome file. Correct sizes obtained with `assembly-stats`.

## EDTA Batch Run

**Script:** `run_edta_liverworts.sh`
**Command per species:** `EDTA.pl --genome genome.fa --cds cds.fa --anno 1 --sensitive 1 --threads 24`
**Status at session end:** Running in background (PID 2591970), on first species (*Blasia pusilla*)
**Master log:** `analysis/liverwort_genomes/edta_run_liverworts.log`
**Per-species log:** `<species_dir>/edta.log`

Auto-fix logic in script:
1. **TIR-Learner/absl-py TypeError** → downgrades `absl-py==0.9.0` and retries
2. **LTR_retriever `*.mod.mod.*` naming bug** → copies files to `*.mod.*` paths and retries

Script is re-entrant: skips species where `*.fa.mod.EDTA.TEanno.gff3` already exists.

## Next Steps
- [ ] Check EDTA run completion: `tail -f analysis/liverwort_genomes/edta_run_liverworts.log`
- [ ] Verify EDTA outputs for each species: `*.EDTA.TElib.fa`, `*.EDTA.TEanno.gff3`, `*.EDTA.intact.gff3`
- [ ] Add species code prefixes to sequence/protein IDs before OrthoFinder integration
- [ ] Add new species to OrthoFinder run
- [ ] Update CLAUDE.md and project-level CLAUDE.md species tables with new entries and code prefixes

## Tags
`#liverworts` `#repeat-annotation` `#EDTA` `#genomics`
