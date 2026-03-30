# Claude Code Session - 2026-03-30

**Project:** hornwort_sex_chromosomes
**Location:** `/media/data/projects/hornwort_sex_chromosomes`

## Summary
Developed a Python script to calculate per-contig k-mer abundance, GC content, and N-fraction (hard-masked repeat content) from genome assemblies. Tested on the *Lunularia cruciata* female HiFi assembly using two levels of repeat masking (MAKER/EDTA and a stricter masker), and used the metrics to classify contigs into primary genome, plastid-derived, and spurious/contaminant categories.

## Work Completed

### Files Created
- `analysis/avg_kmer_abundance.py` — Per-contig k-mer abundance, GC content, N-base count and fraction

### Files Generated (Lunularia female HiFi)
- `analysis/liverwort_genomes/Lunularia_cruciata/female_hifi/kmers21.jf` — Jellyfish k-mer counts (MAKER-masked)
- `analysis/liverwort_genomes/Lunularia_cruciata/female_hifi/kmers21.tsv` — Dumped k-mer counts (MAKER-masked)
- `analysis/liverwort_genomes/Lunularia_cruciata/female_hifi/contig_kmer_stats.tsv` — Per-contig stats (MAKER-masked)
- `analysis/liverwort_genomes/Lunularia_cruciata/female_hifi/kmers21_strict.jf` — Jellyfish k-mer counts (strict masking)
- `analysis/liverwort_genomes/Lunularia_cruciata/female_hifi/kmers21_strict.tsv` — Dumped k-mer counts (strict masking)
- `analysis/liverwort_genomes/Lunularia_cruciata/female_hifi/contig_kmer_stats_strict.tsv` — Per-contig stats (strict masking)

### Commands Executed
```bash
# K-mer counting (MAKER-masked)
jellyfish count -m 21 -s 1G -t 16 -C ERR10480607.bp.p_ctg.fasta.mod.MAKER.masked -o kmers21.jf
jellyfish dump -c -t kmers21.jf > kmers21.tsv

# Run analysis script
python3 analysis/avg_kmer_abundance.py \
    -f ERR10480607.bp.p_ctg.fasta.mod.MAKER.masked \
    -k kmers21.tsv -s 21 \
    -o contig_kmer_stats.tsv

# Same for strict masking with ERR10480607.bp.p_ctg.fasta.masked
```

## Key Decisions & Rationale

- **K-mer counting on genome assembly (not reads):** The jellyfish count was run on the genome FASTA, not raw reads. This means avg_kmer_abundance reflects within-assembly sequence repetitiveness (how often each unmasked k-mer appears across all contigs), NOT read depth/coverage. This is important for interpretation.
- **Canonical k-mers (-C flag):** Jellyfish and the Python script both use canonical k-mers (lexicographically smaller of forward/reverse complement) to avoid double-counting.
- **N-base handling:** K-mers containing N are skipped for abundance calculation but counted separately as `n_kmers`. N-bases are counted directly for `n_fraction`.

## Technical Details

### Script: avg_kmer_abundance.py
Output columns: `contig`, `seq_len`, `num_kmers` (valid, non-N), `n_kmers` (masked/spanning-N), `avg_kmer_abundance`, `gc_content`, `n_bases`, `n_fraction`

### Lunularia female HiFi assembly stats
- Total: 3,611 contigs, 778.2 Mb assembled

### MAKER-masked results — three populations identified

| Population | Contigs | Size | avg_kmer | GC | N-frac |
|---|---|---|---|---|---|
| Primary genome | 44 | 564 Mb | 150–300× | 38–42% | ~36% |
| Plastid-derived | ~1,457 | 116 Mb | 5–20× | <30% | variable |
| Debris/artifacts | ~2,110 | ~98 Mb | <5× | variable | variable |

**Primary genome separation criterion (MAKER-masked):**
```
avg_kmer_abundance >= 150  AND  0.38 <= gc_content <= 0.42
```
This is a clean, unambiguous split with no borderline cases.

**Plastid signal:** GC ~28–29% matches liverwort chloroplast genomes (Marchantia plastid ~28.5% GC). The 116 Mb total reflects high plastid copy number in HiFi data producing many fragmented organellar contigs.

### Strict masking comparison
- N-fraction increased by a uniform **+0.21** on every primary contig → strict masker targets the same repeat class across all chromosomes
- avg_kmer dropped from ~220× to ~90× on primary contigs
- Primary cluster shifted to 50–150× range, overlapping more with non-primary contigs → **MAKER masking is better for contig classification** via this k-mer metric; strict masking is better for annotation/read-depth work

## Next Steps
- [ ] Apply the script to other species with paired U/V genomes to assess contig purity
- [ ] Consider running jellyfish on raw HiFi reads to get true read-depth coverage per contig (complementary to this assembly-based approach)
- [ ] Use primary contig list (44 contigs from MAKER-masked analysis) to filter the Lunularia female assembly for downstream sex chromosome analyses

## Related Files
- `analysis/liverwort_genomes/Lunularia_cruciata/CLAUDE.md` — Lunularia assembly details
- Previous sessions on sex chromosome gene content analysis

## Tags
`#liverwort` `#genome-assembly` `#kmer` `#repeat-masking` `#Lunularia` `#contig-classification`
