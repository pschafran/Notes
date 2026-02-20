# Claude Code Session - 2026-02-20

**Project:** Phaeomegaceros BRAKER annotation — failure analysis
**Location:** `/media/data/projects/hornwort_sex_chromosomes/analysis/Phaeomegaceros_fimbriatus/Phaeomegaceros_14765-5/annotation/braker`
**Duration:** 2026-02-20 (evening)

## Summary

Investigated why the BRAKER3 (ETP mode) gene annotation run for *Phaeomegaceros* failed. The original Jan 28–29 run had **already completed successfully** and produced a high-quality annotation (`bbc/better.gtf`, 3.53% missing BUSCOs). A Feb 20 retry run failed due to a chromosome name mismatch in intermediate GeneMark-ETP files caused by a partial re-run on Feb 4.

## Key Finding: Two Paths, One Directory

`Phaeomegaceros_chiloensis` is a **symlink** → `Phaeomegaceros_fimbriatus`:
```
/media/data/projects/hornwort_sex_chromosomes/analysis/Phaeomegaceros_chiloensis
  -> Phaeomegaceros_fimbriatus
```
All files in both paths are physically identical. The BRAKER run annotated the pilon-polished *Phaeomegaceros chiloensis* genome (`pilon.3.masked.fasta`) using species model `Phaeomegaceros_chiloensis_14765-5`.

## Work Completed

### Files Inspected
- `braker/braker.err` — top-level BRAKER stderr (hardlink shared between fimbriatus and chiloensis paths)
- `braker/braker/braker.log` — BRAKER run log
- `braker/braker/errors/GeneMark-ETP.stderr` — GeneMark-ETP failure details
- `braker/braker/errors/GeneMark-ETP.stdout` — GeneMark-ETP stdout
- `braker/braker/GeneMark-ETP/proteins.fa/nonhc/nonhc.fasta` — empty (0 bytes); root cause
- `braker/braker/GeneMark-ETP/proteins.fa/hc_regions.gtf` — has SHORT chromosome names
- `braker/braker/GeneMark-ETP/data/genome.softmasked.fasta` — has LONG chromosome names
- `braker/braker/GeneMark-ETP/proteins.fa/nonhc/nonhc.trace` — only header, no regions
- `braker/braker/GeneMark-ETP/proteins.fa/nonhc/prothint/log` — ProtHint error: missing genemark.gtf
- `/home/peter/bin/gmetp.pl` — GeneMark-ETP pipeline script (analyzed CreateThis, JoinNonhs functions)

### Final Output Files (from successful Jan 28–29 run)
- `braker/braker/braker.gtf` — primary annotation (Jan 29 00:23)
- `braker/braker/braker.aa` — protein sequences (Jan 29 00:23)
- `braker/braker/braker.codingseq` — coding sequences (Jan 29 00:23)
- `braker/braker/bbc/better.gtf` — **BUSCO-optimized final annotation (3.53% missing BUSCOs)**
- `braker/braker/GeneMark-ETP/genemark.gtf` — GeneMark output (Jan 28 23:36, 27.6 MB)

## Root Cause of Feb 20 Failure

### Failure chain

```
Feb 20 BRAKER retry
  └─ gmetp.pl called
       └─ PrepareGenomeTraining()
            └─ CreateThis("../nonhc.fasta") → returns 0 (file exists)
                 └─ SKIPS probuild entirely
                      └─ nonhc.fasta remains EMPTY
                           └─ pred_m/genemark.gtf never created
                                └─ prothint/prothint.gff never created
                                     └─ printRnaAlternatives.py crashes
                                          └─ gmetp.pl fails → BRAKER fails
```

### Why nonhc.fasta is empty (0 bytes, created Feb 4)

1. **Feb 4 partial re-run**: gmetp.pl regenerated `data/genome.softmasked.fasta` from `pilon.3.masked.fasta`. This file has **LONG chromosome names** (`>PhchiF.S1_pilon_pilon_pilon`).

2. **Chromosome name mismatch**: `proteins.fa/hc_regions.gtf` (from Jan 28) uses **SHORT names** (`PhchiF.S1`, no `_pilon_pilon_pilon` suffix). These came from the original run environment (possibly `/media/lilab/ps997/`).

3. **probuild produces no output**: `probuild --cut nonhc --seq genome.softmasked.fasta --regions hc_regions.gtf` matches region names against FASTA headers — no matches found → no `nonhc_*` region files written to `nonhc/regions/`.

4. **JoinNonhs creates empty file**: The `JoinNonhs()` function in `gmetp.pl` opens `nonhc.fasta` for writing, finds no `nonhc_*` files to join → creates a 0-byte `nonhc.fasta`.

5. **Permanently stuck**: `CreateThis()` in gmetp.pl checks only `if (! -e $fname)` — if the file exists (even if empty), it returns 0 and skips regeneration. Every subsequent retry sees the existing empty file and skips the probuild step.

### Key code paths in gmetp.pl

```perl
# CreateThis - only checks existence, NOT content:
sub CreateThis {
    if ( $force or ( ! -e $fname )) { return 1; }  # create
    else { return 0; }                              # skip
}

# JoinNonhs - creates empty output if no nonhc_N files exist:
sub JoinNonhs {
    open( my $OUT, ">", $out );  # truncates/creates file
    foreach $f (@list) {
        if ( $f =~ /^nonhc_(\d+)$/ ) { ... }  # no matches → empty output
    }
}
```

## BRAKER Run Configuration

```bash
braker.pl \
  --genome pilon.3.masked.fasta \
  --species Phaeomegaceros_chiloensis_14765-5 \
  --prot_seq uniprot_viridiplantae_reference_proteomes.fasta \
  --bam pilon.3.masked.fasta.RNA.bam \
  --threads 24 \
  --busco_lineage viridiplantae_odb10 \
  --useexisting
```

### Genome details
- **Genome**: `pilon.3.masked.fasta` — soft-masked, 148.7 Mb, 319 scaffolds
- **Headers**: `>PhchiF.S1_pilon_pilon_pilon` (LONG format — 3 rounds of pilon polishing appended `_pilon`)
- **RNA-seq BAM**: `pilon.3.masked.fasta.RNA.bam` — aligned to pilon.3 genome (long names)
- **Protein DB**: `uniprot_viridiplantae_reference_proteomes.fasta`
- **BUSCO lineage**: `viridiplantae_odb10`

### GeneMark-ETP statistics (successful Jan 28 run)
- 5,237 HC regions identified
- HC sequence: ~11.2 Mb (7.5% of genome)
- Non-HC sequence: ~137.5 Mb (92.5% of genome)

## BUSCO Summary (from best_by_compleasm.log)
| Annotation | Missing BUSCOs |
|---|---|
| braker.gtf | 19.06% |
| genemark.gtf | 5.18% |
| augustus.hints.gtf | 4.24% |
| **bbc/better.gtf** (final) | **3.53%** |

The final annotation was produced by TSEBRA merging BRAKER + GeneMark + Augustus gene sets.

## How to Fix (if re-run is needed)

### Option A — Use existing output (recommended)
The Jan 29 results are complete. Use `braker/braker/bbc/better.gtf` as the final annotation.

### Option B — Fix the stale intermediate state
Delete the empty/stale files, then fix the chromosome name mismatch before re-running:

```bash
cd braker/braker/GeneMark-ETP/proteins.fa/nonhc/

# 1. Remove the empty/stale files that block regeneration
rm nonhc.fasta nonhc.trace
rm -rf regions/ pred_m/ for_prothint/ prothint/

# 2. Fix chromosome name mismatch in genome.softmasked.fasta
# The hc_regions.gtf uses short names (PhchiF.S1)
# but genome.softmasked.fasta has long names (PhchiF.S1_pilon_pilon_pilon)
# Strip the _pilon_pilon_pilon suffix from FASTA headers:
cd ../../data/
sed -i 's/_pilon_pilon_pilon//' genome.softmasked.fasta

# 3. Re-run gmetp.pl (from the GeneMark-ETP working directory)
cd ../../
perl /home/peter/bin/gmetp.pl \
  --cfg etp_config.yaml \
  --workdir . \
  --bam etp_data/ \
  --cores 24 \
  --softmask
```

### Option C — Clean re-run from scratch
```bash
# Delete the BRAKER output directory and re-run without --useexisting
rm -rf braker/braker/
braker.pl \
  --genome pilon.3.masked.fasta \
  --species Phaeomegaceros_chiloensis_14765-5 \
  --prot_seq uniprot_viridiplantae_reference_proteomes.fasta \
  --bam pilon.3.masked.fasta.RNA.bam \
  --threads 24 \
  --busco_lineage viridiplantae_odb10
```

## Technical Details

### Directory structure note
```
analysis/Phaeomegaceros_chiloensis -> Phaeomegaceros_fimbriatus  (symlink)
```
Both paths point to the same physical directory. The braker.err in both locations are hardlinks (same inode 456337364). The `braker/` subdirectory also shares the same inode (457181904).

### gmetp.pl location
`/home/peter/bin/gmetp.pl` (version supporting ETP mode)

### Key file timestamps
| File | Date | Status |
|---|---|---|
| `rnaseq/stringtie/transcripts_merged.fasta` | Jan 28 22:00 | ✓ Complete |
| `proteins.fa/hc_regions.gtf` | Jan 28 22:27 | ✓ Complete (SHORT names) |
| `GeneMark-ETP/genemark.gtf` | Jan 28 23:36 | ✓ Complete (27.6 MB) |
| `braker.gtf`, `braker.aa` | Jan 29 00:23 | ✓ Complete |
| `data/genome.softmasked.fasta` | Feb 4 13:47 | ⚠ Regenerated (LONG names) |
| `nonhc/nonhc.fasta` | Feb 4 13:49 | ✗ Empty (0 bytes) — root cause |
| `errors/GeneMark-ETP.stderr` (Feb 20 run) | Feb 20 16:48 | ✗ Failure logged |

## Next Steps
- [ ] Decide whether Jan 29 annotation (`bbc/better.gtf`) is sufficient for downstream analysis
- [ ] If re-run needed: apply Option B fix (name mismatch correction) before retrying
- [ ] Note that the `Phaeomegaceros_chiloensis` symlink means all fimbriatus/chiloensis paths are shared — any BRAKER run for fimbriatus genome will need its own separate directory

## Tags
`#genome-annotation` `#BRAKER3` `#GeneMark-ETP` `#Phaeomegaceros` `#hornwort`
