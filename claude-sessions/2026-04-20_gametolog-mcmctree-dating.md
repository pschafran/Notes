# Claude Code Session - 2026-04-20

**Project:** gametolog_ancestry_tests
**Location:** `/media/data/projects/hornwort_sex_chromosomes/analysis/gene_analyses/gametolog_ancestry_tests`
**Duration:** ~07:00 – 07:58

## Summary

Completed MCMCtree Bayesian molecular clock analysis dating the divergence of U/V sex chromosome gametolog pairs in hornworts. Fixed three PAML 4.9j-specific configuration errors (wrong model number, wrong usedata for proteins, missing rate priors), ran all 8 HOGs through step1 and step2, and confirmed convergence with independent replicate chains. The primary finding is that U/V sex chromosomes diverged ~450–490 Ma, substantially predating hornwort crown diversification (~340 Ma).

## Work Completed

### Files Modified
- `build_mcmctree_inputs.py` — Fixed LG model specification (model=2 + aaRatefile), fixed step2 to use usedata=2, added rgene_gamma/sigma2_gamma priors, added print=1
- `CLAUDE.md` — Rewrote with corrected MCMCtree run instructions, PAML gotchas, and full results table

### Files Generated (by script)
- `mcmctree/HOG_OG*/HOG_OG*_focal.phy` — focal-only PHYLIP alignments (already existed)
- `mcmctree/HOG_OG*/step1/out.BV` — gradient/Hessian from usedata=3 run
- `mcmctree/HOG_OG*/step2/mcmc.txt` — 100,000 posterior MCMC samples (rep1)
- `mcmctree/HOG_OG*/step2/step2_out` — summary statistics
- `mcmctree/HOG_OG*/step2_rep2/mcmc.txt` — 100,000 posterior MCMC samples (rep2, seed=-1)

### Commands Executed
```bash
# Step 1 — gradient/Hessian (already completed for OG0000460; run remaining 7)
cd mcmctree/HOG_OGxxxxxxx/step1 && mcmctree mcmctree.ctl

# Step 2 — MCMC dating (run all 8 in parallel)
for hog in ...; do
    cp mcmctree/HOG_$hog/step1/out.BV mcmctree/HOG_$hog/step2/in.BV
    cd mcmctree/HOG_$hog/step2 && mcmctree mcmctree.ctl &
done

# Replicate run (seed=-1 added to ctl)
for hog in ...; do
    cp step1/out.BV step2_rep2/in.BV
    cd step2_rep2 && mcmctree mcmctree.ctl &
done
```

## Key Decisions & Rationale

- **LG model**: `model=2` + `aaRatefile=/usr/lib/paml/data/dat/lg.dat`, NOT `model=7`. In PAML 4.9j, model=7 is cpREV (unavailable) and silently falls back to REVaa. Confirmed by reading `tmp0001.out` which prints the AAML model name.
- **usedata=2 for step2**: `usedata=1` is nucleotide-only in PAML 4.9j. Protein approximate likelihood requires `usedata=2` reading `in.BV` (copied from step1's `out.BV`), not `rst2`.
- **Rate priors required**: Without `rgene_gamma` and `sigma2_gamma` in the ctl file, MCMCtree defaults to zeros, producing a completely stuck chain (all parameters frozen at starting values, lnL constant). Added `rgene_gamma=2 4` (mean rate 0.5 sub/site/100Ma) and `sigma2_gamma=1 10`.
- **print=1**: Without this, MCMCtree does not write individual samples to `mcmc.txt`, making it impossible to compute HPD intervals. Only summary means are written to the outfile.
- **Replicate seeding**: Without `seed=-1` in the ctl, MCMCtree uses a fixed default seed — both replicates were identical. Fixed by prepending `seed=-1` to rep2 ctl files.

## Results: U/V Divergence Dates

Posterior means and 95% HPD intervals (both replicates converged, Δmean < 1 Ma in all cases):

| HOG | Category | Mean (Ma) | HPD lo | HPD hi | ESS |
|-----|----------|-----------|--------|--------|-----|
| OG0000460 | clean UV | 490 | 450 | 519 | 21,311 |
| OG0000985 | clean UV | 481 | 436 | 514 | 45,472 |
| OG0003154 | clean UV (–PhchiM) | 477 | 427 | 513 | 35,087 |
| OG0002735 | clean UV (–LedusM) | 367 | 341 | 413 | 65,899 |
| OG0001301 | messy UV (–Phphy) | 492 | 454 | 519 | 32,278 |
| OG0001152 | messy UV (–PhphyF) | 449 | 396 | 507 | 59,914 |
| OG0002583 | U-only | 340 | 331 | 349 | 89,855 |
| OG0008579 | U-only | 340 | 331 | 349 | 90,356 |

**Primary finding**: U/V sex chromosomes diverged ~450–490 Ma, predating hornwort crown (~340 Ma) by ~130–150 Ma.

**Caveats**:
- Most HPD intervals push against the upper soft bound (505 Ma); true divergence could be older
- OG0002735 outlier (367 Ma) likely reflects absent LedusM compressing root age
- U-only HOGs recover calibration A prior by construction; no independent dating power

## Technical Details

### PAML 4.9j protein model numbers (seqtype=2, for future reference)
- model=0: Poisson
- model=1: EqualInput
- model=2: Empirical (uses aaRatefile; defaults to JTT if empty)
- model=3: Empirical+F
- model=4,7,10,11,12: unavailable in this build
- model=8: REVaa_0
- model=9: REVaa (full, 189 free parameters)

### MCMCtree approximate likelihood workflow (proteins)
1. Step 1 ctl: `usedata=3`, outputs `out.BV`
2. `cp out.BV ../step2/in.BV`
3. Step 2 ctl: `usedata=2`, reads `in.BV`; requires `rgene_gamma`, `sigma2_gamma`, `print=1`, `seed=-1`

### Checking the model name
```bash
grep "^Model:" mcmctree/HOG_OG0000460/step1/tmp0001.out
# Model: Empirical,  (/usr/lib/paml/data/dat/lg.dat) dGamma (ncatG=4) ns =   9  ls = 797
```

## Next Steps
- [ ] Plot posterior distributions for all UV HOGs (violin or ridgeline plots)
- [ ] Consider relaxing upper root bound (currently 505 Ma) to see if data want an older date
- [ ] Write up methods and results section

## Related Files
- Previous session: `2026-04-17_hornwort-ks-analysis.md`
- CLAUDE.md: `analysis/gene_analyses/gametolog_ancestry_tests/CLAUDE.md`

## Tags
`#hornwort` `#sex-chromosomes` `#molecular-clock` `#mcmctree` `#paml` `#bayesian`
