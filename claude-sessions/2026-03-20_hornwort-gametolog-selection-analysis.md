# Claude Code Session - 2026-03-20

**Project:** hornwort-sex-chromosomes / gametolog selection analysis
**Location:** `/media/data/projects/hornwort_sex_chromosomes/analysis/synteny/circos_hornwort_sex_acc_chr/`
**Session:** Continuation of a large multi-session analysis (context was compacted from a previous run)

---

## Summary

Extended and finalised the gametolog Ka/Ks selection analysis for hornwort sex chromosomes, adding Anthoceros angustus intra-genomic U-V paralog pairs computed from ODP RBH data, improving transcript assignments for pairs where OrthoFinder primary transcripts were suboptimal, and producing a comprehensive gametolog Ka/Ks table. Finished with a literature review of reproductive roles and interactions among the 16 candidate genes under purifying selection, and a KEGG pathway co-membership query.

---

## Work Completed

### Files Created (by Claude)
- `scripts/10_anang_dof_ks.py` — Ka/Ks for Anthoceros angustus DOF zinc-finger U-V paralog pair (chrU_g0091 vs chrV_g0007); tests all 4 transcript combinations; uses mafft + yn00
- `scripts/11_gametolog_table.py` — Comprehensive gametolog Ka/Ks table across all 5 species; applies transcript overrides; notes where primary transcripts were replaced; writes to `output/gametolog_kaks_table.tsv`
- `scripts/12_anang_gametolog_ks.py` — Ka/Ks for all 14 Anang chrU-chrV ODP RBH pairs; tests all transcript combinations; maps pairs to OrthoFinder OGs; writes to `output/anang_dof_ks/anang_gametolog_results.tsv`
- `notes/20260320_gametolog_literature_review.md` — Full literature review of 16 gametolog genes (reproductive roles, bryophyte evidence, cited sources, interactions)

### Files Modified (by Claude)
- `scripts/08_selection_analysis.py` — Major updates:
  - Added transcript override logic (loads `transcript_check_results.tsv`; for "IMPROVED" pairs substitutes best_Ks/Ka/omega)
  - Added Anang gametologs via ODP RBH (loads `anang_gametolog_results.tsv` + DOF pair)
  - Changed panel C from median-only dot plot to **per-species dots** with species-specific colors/markers and grey median bar
  - Added `OG_DESCRIPTION_OVERRIDE` dict with curated functional descriptions for all 16 OGs
  - Added "A. angustus" to SP_LABELS

### Output Files Updated
- `output/selection_candidates.tsv` — 20 OGs (19 purifying, 1 neutral), improved descriptions
- `output/selection_analysis.png/.pdf` — Updated figure with per-species panel C
- `output/gametolog_kaks_table.tsv` — 89-row comprehensive table (59 included, 10 Ks>3.0, 2 Ks<0.05, 15 n_sp<2, 3 omega missing)
- `output/anang_dof_ks/anang_gametolog_results.tsv` — 14 Anang RBH pair results (9 passed Ks filter)
- `output/anang_dof_ks/anang_dof_results.tsv` — DOF pair all-combo results

---

## Key Technical Details

### Anang naming convention fix (critical)
Anang gene names have the form `AnangRef2_chrU_g0091.t1` — no inner dots in the gene base name. Using `.split(".")[:2]` to strip the transcript suffix returns the unchanged string; the correct method is:
```python
re.sub(r'\.t\d+$', '', transcript_name)
```
This bug caused all Anang OG lookups to fail silently (returning "-") until fixed.

### gffread output header stripping
gffread strips the `AnangRef2_` prefix from sequence names in output FASTA (e.g. outputs `chrU_g0022.t1` not `AnangRef2_chrU_g0022.t1`). Script 12 uses a `seq_lookup()` helper that tries both forms.

### Transcript override in script 08
Two IMPROVED pairs (Phphy G055300, Phchi G001400) had primary Ks > KS_MAX but best Ks ≤ KS_MAX. Fixed by loading `sex_all` (no Ks pre-filter) and applying the Ks filter per-row after override lookup. Background distribution panels A/B still use pre-filtered primary-transcript data.

### OG_DESCRIPTION_OVERRIDE dict
```python
OG_DESCRIPTION_OVERRIDE = {
    "OG0003907":            "USP7/HAUSP (ubiquitin-specific protease)",
    "OG0007255":            "AGO1 (Argonaute 1)",
    "OG0000643":            "Cytokinin receptor kinase (CHASE-type)",
    "OG0000564":            "DDX39/UAP56 (mRNA export helicase)",
    "OG0010209":            "SnRK2 (ABA-activated protein kinase)",
    "OG0009804":            "EDR1-like MAP3K",
    "OG0003055":            "bHLH transcription factor",
    "OG0004941":            "C2H2 zinc finger protein",
    "OG0002767":            "ClpB/HSP101 (heat shock disaggregase)",
    "OG0010581+OG0011317":  "TCP transcription factor",
    "OG0002662+OG0008508":  "HD-ZIP transcription factor",
    "OG0011295+OG0011420":  "DOF zinc finger transcription factor",
    "OG0000346":            "TOPLESS/TPL transcriptional co-repressor",
    "OG0002632":            "UBXN7 (ubiquitin regulatory X domain)",
    "OG0006195":            "Importin-5 (IPO5/RanBP5)",
    "OG0006678":            "UBE2D (ubiquitin-conjugating enzyme E2)",
}
```

### yn00 output parsing (key indices)
```
# Yang & Nielsen section columns:
# seq1 seq2 S N t kappa omega dN±SE dS±SE
# indices: 0  1  2 3 4 5     6     7  8   9  10
omega = float(nums[6])
Ka    = float(nums[7])
Ks    = float(nums[9])
```

### Anang DOF pair best result
- Best combination: chrU_g0091.t1 × chrV_g0007.t1 (Ks=2.573, Ka=0.700, ω=0.272)
- OrthoFinder primary (t1×t2): Ks=2.811

---

## Literature Review Findings

Full review in: `notes/20260320_gametolog_literature_review.md`

### Proteins with strongest bryophyte gametophyte evidence
| Protein | Key evidence | Species |
|---|---|---|
| Cytokinin receptor | Female gametophyte abortion in *cki1*; gametophore initiation | *P. patens*, *Marchantia* |
| bHLH | Gametangia differentiation (MpBNB); sperm specification (DUO1) | *Marchantia* |
| TOPLESS | Sperm formation (DAZ1/2); 2D→3D gametophore transition | *P. patens* |
| CDF/DOF | Seasonal sexual reproduction regulation (CDL genes) | *P. patens* |

### Key protein–protein interactions within the set
- **AGO1 → TCP** (miR319): direct, well documented
- **TOPLESS ↔ CDF/DOF**: CDF1 physically recruits TPL to FT/CO promoters
- **TOPLESS ↔ bHLH** (DAZ1/2/DUO1): required for sperm formation
- **TOPLESS ↔ JAZ ↔ bHLH/MYC2**: jasmonate pathway, male fertility
- **SnRK2 → 14-3-3**: phosphorylation-dependent scaffold interaction
- **UBP26 ↔ UBE2D**: H2Bub1 writer/eraser pair on same chromatin mark

### KEGG shared pathways
Only 8 of 16 genes had KO IDs from eggNOG; two shared pathways found:
| Pathway | Genes |
|---|---|
| map03083 — Polycomb repressive complex | USP7/UBP26 (K11838) + UBE2D (K06689) |
| map04075 — Plant hormone signal transduction | Cytokinin receptor (K14489) + SnRK2 (K14498) |

Remaining 8 genes (TCP, bHLH, TOPLESS, CDF/DOF, 14-3-3, EDR1 MAP3K, IPO5, UBXN7) not queried — no KO IDs in hand.

---

## Key Decisions & Rationale

- **Per-species dots in panel C** (not medians): with only 4–5 species per OG, showing individual species values is more informative and avoids false precision from median of a small N
- **OG_DESCRIPTION_OVERRIDE** rather than modifying eggNOG output: allows manually curated functional descriptions without altering source data; applied as a final post-processing step in script 08 section 4
- **Anang added via ODP RBH** (not OrthoFinder Ks TSV): Anang has a combined U+V genome with no separate male/female genotype, so there is no between-genotype Ks TSV; intra-genomic chrU vs chrV pairs identified via ODP reciprocal best hits are the correct substitute

---

## Next Steps
- [ ] Query KEGG pathways for remaining 8 genes (TCP, bHLH, TOPLESS, CDF/DOF, 14-3-3, EDR1 MAP3K, IPO5, UBXN7) — need KO IDs from eggNOG annotations
- [ ] Consider adding literature references to `output/selection_candidates.tsv`
- [ ] Assess whether to add KEGG pathway co-membership to the literature review document

## Related Files
- Previous session (large, context-compacted): scripts 08–12 in `circos_hornwort_sex_acc_chr/scripts/`
- eggNOG annotations: `analysis/functional_annotations/`
- ODP RBH source: `odp_bryophytes_20260318/odp/step2-figures/synteny_nocolor/`
- Transcript check results: `circos_hornwort_sex_acc_chr/output/transcript_check_results.tsv`

## Tags
`#hornworts` `#sex-chromosomes` `#gametologs` `#Ka-Ks` `#selection` `#python` `#literature-review`
