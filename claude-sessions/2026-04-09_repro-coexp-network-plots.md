# Claude Code Session - 2026-04-09

**Project:** repro_coexp
**Location:** `/media/data/resources/marpolbase_expression_data/repro_coexp`
**Model:** claude-sonnet-4-6
**Session ID:** fe7d88bc-7bcb-49e2-bc0e-ddd494404ced

## Summary

Audited genes of interest from a full OG table (`tmp`) against the Marchantia reproductive co-expression network, correcting an error in community assignments for OG0000031. Regenerated OG assignment files from the complete table, added three new genes (BPCV/MpVg00350, BPCU/MpUg00370, BONOBO/Mp3g23300), fixed PDF font embedding for Adobe Acrobat compatibility, and added U/V chromosome visual markers to all network plots. Also extracted the community 9 node list and investigated why certain genes of interest are absent from the network.

## Work Completed

### Files Modified
- `repro_coexp/plot_communities.py` — added `pdf.fonttype=42` (Acrobat fix), U/V chromosome node symbols and legend entries

### Files Created
- `repro_coexp/OG*_community_assignment.tsv` (12 files) — regenerated from full `tmp` table via `/tmp/regen_og_assignments.py`
- `repro_coexp/OG_BPCVBPCU_community_assignment.tsv` — BPCV (MpVg00350) + BPCU (MpUg00370)
- `repro_coexp/OG_BONOBO_community_assignment.tsv` — BONOBO (Mp3g23300)
- `repro_coexp/community_*_network.pdf/png` (9 communities × 2 formats) — updated plots
- `repro_coexp/community_9_plotted_nodes.tsv` — full node list for community 9 subgraph

### Commands Executed
```bash
python3 /tmp/regen_og_assignments.py   # regenerate OG files
python3 plot_communities.py            # re-run all plots
```

## Key Decisions & Rationale

- **TrueType fonts (pdf.fonttype=42):** Matplotlib defaults to Type 3 bitmap fonts which Adobe Acrobat rejects. Setting fonttype=42 embeds TrueType fonts and resolves the issue.
- **U/V markers as centred text, not offset labels:** Symbols sit inside the node circle so they remain readable at small node sizes without cluttering the layout.
- **OG files must start with "OG":** The glob pattern `OG*_community_assignment.tsv` requires this prefix. New gene files (BPCVBPCU, BONOBO) must follow this convention.

## Key Findings

### OG0000031 was undercalled
The `tmp` table listed "Yes (MpVg00010.1)" but 3 of 4 Mp paralogs are in enriched communities:
- Mp4g15710 → community 7 (pan_male)
- MpUg00340 → community 9 (pan_reproductive)
- MpVg00010 → community 7 (pan_male)

The paralog split across communities suggests sex-chromosome-driven subfunctionalization: V-linked gene joins male module, U-linked gene joins pan-reproductive/female module.

### Community 9 is a U-chromosome female module
Of 34 plotted nodes, 28 are MpUg* genes, almost all elevated specifically in archegonia and archegoniophore. The neighbourhood of the three genes of interest (MpUg00340, MpUg00040/CULT1, MpUg00370/BPCU) is dominated by contiguous U-chromosome genes with female-biased expression.

### Absent genes: biology not data quality
Genes missing from the network were investigated via correlation graph presence and best HRR to any seed gene:

| Gene | OG | In corr. graph | Best HRR to seed |
|---|---|---|---|
| Mp4g09475 | OG0000002 | **No** | — |
| Mp8g15630 | EDR1 | Yes (best HRR=16) | 336 |
| Mp1g18540 | CDF6 | Yes (best HRR=6) | 147 |
| Mp1g29640 | HK3 | Yes (best HRR=7) | 390 |
| Mp2g15090 | HK3 | Yes (best HRR=8) | 263 |
| Mp3g15430 | USP24 | Yes (best HRR=11) | 156 |

Mp4g09475 is not expressed. All others are expressed and well co-expressed with non-reproductive genes, but their expression patterns don't couple to the reproductive programme.

## Next Steps
- [ ] Correct OG0000031 entry in `tmp` table (should list Mp4g15710, MpUg00340, MpVg00010)
- [ ] Optionally investigate which gene modules the absent genes belong to (EDR1, CDF6, HK3 paralogs)
- [ ] Consider whether a relaxed HRR threshold could capture weakly reproductive-connected genes

## Related Files
- Lab notebook entry: `/media/data/projects/repro_coexp/notebook/entry_2026-04-09_session2.md`
- Previous session summary (network construction): `/home/peter/.claude/projects/-media-data-resources-marpolbase-expression-data/fe7d88bc-7bcb-49e2-bc0e-ddd494404ced.jsonl`
- `repro_coexp/tmp` — full gene-of-interest table (OG, Mp orthologs, community annotations)

## Tags
`#coexpression` `#marchantia` `#network` `#python` `#visualization` `#sex-chromosomes`
