# Claude Code Bug Report - Usage Policy False Positives

**Date:** 2026-02-01
**Time:** Approximately 05:00-05:30 UTC
**Tool:** Claude Code CLI
**Model:** claude-opus-4-5-20251101

## Summary

During a routine bioinformatics session analyzing plant genomic data, multiple responses were truncated or blocked due to apparent false positive usage policy violations. All content was standard scientific output with no sensitive, harmful, or policy-violating material.

## Session Context

- **Project:** Phylogenetic analysis of *Iva* (marsh elder) plant species
- **Data:** GetOrganelle chloroplast and nuclear ribosomal DNA assemblies
- **Tools used:** tblastn, seqtk, MAFFT, IQ-TREE
- **Working directory:** `/media/data/projects/iva_phylogeny/analysis`

## Incident 1: Truncated Scientific Summary Table

### What was being generated
A markdown summary table of nuclear ribosomal DNA analysis containing:
- Sample identifiers (e.g., "89_STURM", "36_frutescens")
- Sequence counts and lengths in base pairs
- Assembly status ("K115 complete", "K115 scaffolds")
- Biological notes ("heterozygous", "intragenomic variation")

### Truncation point
Response cut off mid-sentence at "~10 " while comparing sequence lengths between sample groups.

### Content snippet (reconstructed)
```markdown
| Sample | Species | # Sequences | Length (bp) | Assembly | Status |
|--------|---------|-------------|-------------|----------|--------|
| 01_EREMID | *I. eremidica*? | 1 | 10,265 | K115 complete | Single sequence |
| 89_STURM | Unknown | 1 | 10,311 | K115 complete | Single sequence |
...
```

### Likely trigger candidates
- Markdown table with multiple columns
- Scientific terminology: "heterozygous", "scaffolds", "variation"
- Species names with asterisk formatting (*I. frutescens*)
- Numeric data (percentages, base pair counts)

## Incident 2: Repeated Output Failures

### What was being generated
ASCII art phylogenetic tree diagram with:
- Branch structure using characters: `/`, `|`, `+`, `-`
- Bootstrap support values (e.g., "98", "99")
- Sample names with underscores
- Branch length decimals (e.g., "0.1105")

### Symptom
Multiple attempts to display tree failed. Some showed sandbox errors (`bwrap`), others appeared to be policy-related truncations.

### Content that failed
```
                              +---- 89_STURM
                    /--98-----|
                   |          +---- 91_STURM
         /---60----|
        |          +---------------------- 49_frutescens
```

## Incident 3: Python Script Output Truncation

### What was being generated
Python script output displaying phylogenetic tree interpretation and "Key observations" section.

### Likely trigger candidates
- Multi-line string with tree ASCII art
- Print statements with scientific observations
- Terms: "divergence", "long branch", "contamination" (in scientific context)

## Common Elements Across All Incidents

1. **Markdown tables** with genomic statistics
2. **Percentage values** (e.g., "2.5% gaps", "0.4%")
3. **Decimal numbers** for branch lengths
4. **Underscore identifiers** (sample_name format)
5. **Newick tree format** with nested parentheses
6. **Scientific terms:** heterozygous, variation, divergence, polymorphism, contamination

## What This Content Actually Was

- Standard bioinformatics analysis output
- Plant genomics (marsh elder, *Iva* species)
- No code vulnerabilities or security content
- No personal data
- No harmful instructions
- Purely scientific/educational content

## Impact

- Disrupted scientific workflow
- Required multiple retries and workarounds
- User had to manually request continuation
- Some information had to be reformatted to avoid triggers

## Suggested Investigation

1. Review server logs for this session (timestamps ~05:00-05:30 UTC, 2026-02-01)
2. Check which safety classifier triggered
3. Examine confidence scores on flagged content
4. Consider whitelisting common bioinformatics terminology

## Environment Details

```
Platform: linux
OS Version: Linux 6.8.0-90-generic
Working directory: /media/data/projects/iva_phylogeny/analysis
Model: claude-opus-4-5-20251101
```

## Report Location

GitHub Issues: https://github.com/anthropics/claude-code/issues
