# Claude Code Session - 2026-02-28

**Project:** iva_phylogeny / analysis
**Location:** `/media/data/projects/iva_phylogeny/analysis`

## Summary

Populated `phylogeny_paper/cp/`, `phylogeny_paper/nr/`, and `phylogeny_paper/mt/` with
oriented, properly named sequences for all 23 paper-subset samples. Investigated
anomalous assemblies for hayesiana_93 (nrDNA) and texensis_48 (mitochondria).
Produced a cross-marker issue summary for all samples.

---

## Work Completed

### Files Created

- `phylogeny_paper/collect_cp_nr.py` — collects and orients cp/nr sequences
- `phylogeny_paper/cp/{sample_id}_cp.fasta` — 23 oriented cp sequences (~152 kb each)
- `phylogeny_paper/cp/all_cp.fasta` — combined cp FASTA
- `phylogeny_paper/cp/cp_manifest.tsv` — source file + orientation notes per sample
- `phylogeny_paper/nr/{sample_id}_nr.fasta` — 23 oriented nr sequences (6.8–11.5 kb)
- `phylogeny_paper/nr/all_nr.fasta` — combined nr FASTA
- `phylogeny_paper/nr/nr_manifest.tsv` — source file + orientation notes per sample
- `phylogeny_paper/mt/{sample_id}_mt.fasta` — 23 per-sample mt scaffold files
- `phylogeny_paper/mt/mt_sample_manifest.tsv` — scaffold counts + sizes per sample

### Key Commands / Steps

```bash
# CP + NR collection
python3 phylogeny_paper/collect_cp_nr.py

# MT individual files (split from existing all_mt_scaffolds.fasta)
# done inline with Python — see script below

# hayesiana_93 nr: extracted GFA segment 23486 (8243 bp conserved core)
# from: hayesiana/93_STURM_251201adq30ft/hayesiana_CA2_93_nr/
#       embplant_nr.K55.complete.graph1.selected_graph.gfa
```

---

## CP Collection Details

**Script:** `phylogeny_paper/collect_cp_nr.py`

**Orientation logic:**
1. tblastn with `ndhD.faa` → select isoform where ndhD is forward (sstart < send)
2. tblastn with `rbcL.faa` → if rbcL also forward in that isoform, LSC is inverted → reverse-complement
3. STURM samples: files already labeled `forward_isoform.fasta` / `reverse_isoform.fasta`,
   but tblastn was used uniformly regardless

**Special cases:**
- `imbricata_67`: K55 assembly, 10 paths / 3 repeat patterns → restricted candidates to repeat_pattern1
- 5 samples needed LSC-inversion revcomp: annua_01, texensis_48, hayesiana_92, frutescens_58 (and hayesiana_92)

**Results:** 23/23 sequences collected, 151,707–153,030 bp each.
frutescens_49 is largest (153,030 bp) due to known ~850 bp IR expansion.

---

## NR Collection Details

**Orientation logic:** barrnap `--kingdom euk` → select 18S strand;
if 18S on `-` strand → reverse-complement. Convention: 18S on `+` strand.

**Special cases:**
- `hayesiana_93`: 1,000 paths / 72 repeat patterns → probable allopolyploid with ≥2 divergent
  rDNA families (IGS 55–70% divergent between families). Extracted large single-copy GFA
  segment 23486 (8,243 bp = 18S–ITS1–5.8S–ITS2–26S, no IGS). Already 18S-forward.
- `microcephala_57`: graph1 = 63 bp (fragment), used graph2 (K55, 7,212 bp)
- `angustifolia_77`: 4 repeat-pattern files (high intragenomic variation), used repeat_pattern1
- `frutescens_86`: 4 sequences (K55), used graph1.1 (first non-repeat-pattern file)

**Results:** 23/23 sequences collected, 6,806–11,548 bp each.
hayesiana_93 has no IGS (only coding core). Several frutescens are scaffold assemblies.

---

## MT Collection Details

Per-sample files split from existing `all_mt_scaffolds.fasta` (produced by `collect_mt.sh`).
Headers preserved as `sample_id|scaffold_N`. No orientation applied (no standard convention).

| Outlier | Issue |
|---------|-------|
| texensis_48 | 26 scaffolds, only 102 kb total (others ~215–343 kb). GetOrganelle failure: multiple isolated components, 14.6× coverage, timeout. Fallback scaffolding. |
| axillaris_89, frutescens_36, frutescens_49 | Single-sequence repeat_pattern1.1 files (216/360/40 total paths). Excluded from cluster_01 alignment. |
| frutescens_58 | K85 assembly, 215 kb (lower than ~238 kb for other frutescens). |
| imbricata 67/72/85 | Largest scaffold only 36 kb (vs 51–73 kb in other species). Consistent within genus, likely structural. |

---

## hayesiana_93 NR Investigation

Full pairwise distance analysis of all 1,000 nrDNA paths:
- All 1,000 unique sequences, all 11,389 bp
- Positions 0–8,191: 100% identical (18S–ITS1–5.8S–ITS2–26S rRNA coding)
- Positions 8,191–11,389: ~90% variable (IGS region)
- Assembly graph: 1 large segment (8,243 bp) + 25 small IGS segments (57–1,008 bp), 35 links

Ward hierarchical clustering of IGS sequences → 3 clusters (n=317, 345, 338):
- Medoid pairwise distances: C1 vs C2 = 70.4%, C1 vs C3 = 56.5%, C2 vs C3 = 55.8%
- Quick IQ-TREE (GTR+G): families 2+3 sister to hayesiana_92 (native hayesiana rDNA);
  family 1 deeply divergent, associated with non-hayesiana Iva clade
- Interpretation: hayesiana_93 is likely a hybrid/allopolyploid carrying a native hayesiana
  rDNA family and a divergent family from another Iva lineage
- For phylogeny: using the conserved single-copy core (GFA seg 23486) avoids the IGS ambiguity
  while retaining full 18S–26S coding information

---

## Cross-Marker Quality Summary

### Major issues
- **texensis_48 MT**: fragmented/incomplete assembly
- **hayesiana_93 NR**: allopolyploid signal; IGS excluded from paper sequence
- **frutescens_49 CP + NR**: long branches in previous phylogenies; enlarged IR in cp

### Moderate issues
- **imbricata_67 CP**: K55 assembly, 10 paths / 3 repeat patterns
- **axillaris_89, frutescens_36, frutescens_49 MT**: repeat-path assemblies, excluded from cluster_01 alignment
- **angustifolia_77, frutescens_86 NR**: multiple rDNA copies; single representative chosen

### Minor / notable
- Several frutescens NR are scaffold (not complete circle) assemblies: angustifolia_46, cultigen_133, frutescens_36, frutescens_58, frutescens_63, frutescens_68, frutescens_71, frutescens_94, asperifolia_60
- frutescens_68 NR: shortest at 6,806 bp
- xanthifolia_64 + imbricata_72 NR: K85 assemblies
- frutescens_58 MT: K85 assembly

### Cleanest samples (no caveats)
annua_01, frutescens_55, hayesiana_92, imbricata_85

---

## Next Steps

- [ ] Align `phylogeny_paper/cp/all_cp.fasta` (MAFFT) and run IQ-TREE (23-sample cp phylogeny)
- [ ] Align `phylogeny_paper/nr/all_nr.fasta` (MAFFT) and run IQ-TREE (23-sample nr phylogeny)
- [ ] Run IQ-TREE on existing `phylogeny_paper/mt/cluster_01/cluster_01_aligned.fasta` (20-sample mt phylogeny)
- [ ] Root all three trees with appropriate outgroup
- [ ] Compare cp/nr/mt topologies for congruence; flag discordant samples
- [ ] Decide whether to include frutescens_49 or note it as outlier
- [ ] Consider flagging texensis_48 mt and hayesiana_93 nr in paper methods
- [ ] Investigate frutescens_49 long-branch cause (contamination? genuine divergence?)
- [ ] Optional: retry texensis_48 mt with --max-rounds 50 or --disentangle-time-limit 7200

## Related Files

- Previous session: check `~/Notes/claude-sessions/` for earlier Iva sessions
- Memory file: `/home/peter/.claude/projects/-media-data-projects-iva-phylogeny-analysis/memory/MEMORY.md`
- Collect script: `/media/data/projects/iva_phylogeny/analysis/phylogeny_paper/collect_cp_nr.py`
- MT collect script: `/media/data/projects/iva_phylogeny/analysis/phylogeny_paper/mt/collect_mt.sh`

## Tags
`#phylogenomics` `#GetOrganelle` `#Iva` `#chloroplast` `#nrDNA` `#mitochondria` `#allopolyploidy`
