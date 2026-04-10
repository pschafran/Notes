#!/usr/bin/env python3
"""
Rosner's Generalized ESD test for gene-depleted / repeat-enriched outlier
chromosomes in liverwort genome composition statistics.

Runs one-sided lower ESD on gene_pct (gene depletion) and one-sided upper
ESD on repeat_pct (repeat enrichment) for each species independently, using
only scaffolds >= 1 Mb. r=5 accommodates sex chromosomes that may be
fragmented across multiple scaffolds.

A secondary modified Z-score filter (|Zi| > MOD_Z_THRESH) is applied as a
minimum effect-size gate to address two known ESD failure modes:
  1. Masking: multiple marginal outliers inflate SD, preventing detection of
     the most extreme value (e.g. Fissidens scaffold_19 at 0.9% gene content).
  2. Tight-cluster false positives: after removing a very extreme outlier, the
     remaining autosomes have tiny SD, causing minor deviations to exceed the
     (now low) critical value at later steps (e.g. Ceratodon autosomes after
     removing the U chromosome).

A scaffold is flagged as a CONFIDENT outlier only if it passes both ESD and
the modified Z-score criterion. ESD-only and modZ-only calls are also reported
separately in the detail output.

Outputs:
  esd_outlier_results.tsv  — one row per tested scaffold with all statistics
  esd_outlier_summary.tsv  — one row per species summarising outlier calls

Usage:
  python3 run_esd_outlier_test.py

References:
  Rosner (1983) Technometrics 25(2):165-172
  https://www.itl.nist.gov/div898/handbook/eda/section3/eda35h3.htm
  Iglewicz & Hoaglin (1993) How to Detect and Handle Outliers, ASQC
"""

import os
import csv
import numpy as np
from scipy import stats

BASE           = "/media/data/projects/hornwort_sex_chromosomes/analysis/liverwort_genomes"
MIN_SIZE       = 1_000_000   # 1 Mb minimum scaffold size
R_MAX          = 5           # upper bound on number of outliers (r in Rosner 1983)
ALPHA          = 0.05
MOD_Z_THRESH   = 3.5         # Iglewicz & Hoaglin threshold for modified Z-score


# ---------------------------------------------------------------------------
# Core statistical functions
# ---------------------------------------------------------------------------

def rosner_esd(values, r, alpha, side='lower'):
    """
    Rosner's Generalized ESD test for up to r outliers.

    Parameters
    ----------
    values : array-like of float
    r      : int   — maximum number of outliers to test for (r in Rosner 1983)
    alpha  : float — experiment-wise significance level
    side   : 'lower' (detect unusually small values)
             'upper' (detect unusually large values)

    Returns
    -------
    List of dicts (one per removal step, ordered 1..r_eff), each with:
      step        — removal step (1 = most extreme)
      orig_idx    — index in the original values array
      value       — the value removed at this step
      Ri          — ESD test statistic
      lambda_i    — critical value at this step
      is_outlier  — True if declared outlier by backwards ESD inference

    Reference: Rosner (1983) Technometrics 25(2):165-172;
               https://www.itl.nist.gov/div898/handbook/eda/section3/eda35h3.htm
    """
    data = np.array(values, dtype=float)
    n = len(data)

    # Need at least r+2 points to have >= 2 remaining after r removals
    r = min(r, n - 3)
    if r < 1:
        return []

    working     = data.copy()
    working_idx = np.arange(n)
    step_results = []

    # Step 1: compute r test statistics iteratively
    for i in range(1, r + 1):
        mean = np.mean(working)
        sd   = np.std(working, ddof=1)
        if sd == 0:
            break

        if side == 'lower':
            pos = int(np.argmin(working))
            Ri  = (mean - working[pos]) / sd
        else:   # upper
            pos = int(np.argmax(working))
            Ri  = (working[pos] - mean) / sd

        current_n = len(working)   # = n - i + 1

        step_results.append({
            'step':      i,
            'orig_idx':  int(working_idx[pos]),
            'value':     float(working[pos]),
            'Ri':        float(Ri),
            'current_n': current_n,
        })

        working     = np.delete(working,     pos)
        working_idx = np.delete(working_idx, pos)

    # Step 2: compute critical values
    # At step i with current_n = n - i + 1 points:
    #   df       = current_n - 2
    #   p        = alpha / current_n    (one-sided Bonferroni)
    #   t_crit   = t_{p, df}
    #   lambda_i = ((current_n - 1) * t_crit)
    #              / sqrt((current_n - 2 + t_crit^2) * current_n)
    for sr in step_results:
        current_n = sr['current_n']
        df = current_n - 2
        if df <= 0:
            sr['lambda_i'] = np.nan
            continue
        p      = alpha / current_n
        t_crit = stats.t.ppf(1.0 - p, df)
        sr['lambda_i'] = ((current_n - 1) * t_crit) / \
                          np.sqrt((current_n - 2 + t_crit**2) * current_n)

    # Step 3: backwards inference — largest i where R_i > lambda_i
    n_outliers = 0
    for sr in step_results:
        li = sr.get('lambda_i', np.nan)
        if not np.isnan(li) and sr['Ri'] > li:
            n_outliers = sr['step']

    for sr in step_results:
        sr['is_outlier'] = (sr['step'] <= n_outliers)

    return step_results


def modified_zscore(values, side='lower'):
    """
    Iglewicz & Hoaglin modified Z-score using median and MAD.
    Returns array of scores; large positive values = upper-tail outliers,
    large negative = lower-tail outliers.

    M_i = 0.6745 * (x_i - median) / MAD
    """
    data = np.array(values, dtype=float)
    med  = np.median(data)
    mad  = np.median(np.abs(data - med))
    if mad == 0:
        # Fall back to mean/SD if MAD is zero (all values identical)
        sd = np.std(data, ddof=1)
        return (data - med) / sd if sd > 0 else np.zeros_like(data)
    return 0.6745 * (data - med) / mad


# ---------------------------------------------------------------------------
# I/O helpers
# ---------------------------------------------------------------------------

def load_composition(tsv_path):
    rows = []
    with open(tsv_path) as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            try:
                rows.append({
                    'scaffold':   row['scaffold'],
                    'size_bp':    int(row['size_bp']),
                    'n_genes':    int(row['n_genes']),
                    'gene_pct':   float(row['gene_pct']),
                    'repeat_pct': float(row['repeat_pct']),
                })
            except (ValueError, KeyError):
                continue
    return rows


def find_species_tsvs():
    """
    Returns dict: species_label -> list of TSV paths.
    Walks one level of subdirectories to find composition_composition_stats.tsv
    files. Skips composition_test_composition_stats.tsv (duplicate test files).
    Species without stats files (Lunularia, Marchantia inflexa, Radula) are
    silently omitted.
    """
    species = {}
    for d in sorted(os.listdir(BASE)):
        full = os.path.join(BASE, d)
        if not os.path.isdir(full):
            continue

        # Check top-level first
        std = os.path.join(full, "composition_composition_stats.tsv")
        if os.path.isfile(std):
            species[d] = [std]
            continue

        # Check one level of subdirectories (e.g. Marchantia_quadrata/Mquadrata_genome_annot/,
        # Ricciocarpos_natans/Fu_et_al_2025/)
        for sub in sorted(os.listdir(full)):
            subfull = os.path.join(full, sub)
            if not os.path.isdir(subfull):
                continue
            p = os.path.join(subfull, "composition_composition_stats.tsv")
            if os.path.isfile(p):
                species[d] = [p]
                break   # use first match per species

    return species


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    species_tsvs = find_species_tsvs()

    detail_rows  = []
    summary_rows = []

    for species, tsv_list in species_tsvs.items():
        all_rows = []
        for tsv in tsv_list:
            all_rows.extend(load_composition(tsv))

        # Filter to >= 1 Mb
        large = [r for r in all_rows if r['size_bp'] >= MIN_SIZE]
        n = len(large)

        if n < 4:
            summary_rows.append({
                'species': species, 'n_sequences_tested': n,
                'gene_confident': 'n/a (too few sequences)',
                'repeat_confident': 'n/a (too few sequences)',
                'both_confident': 'n/a',
            })
            continue

        gene_vals   = [r['gene_pct']   for r in large]
        repeat_vals = [r['repeat_pct'] for r in large]
        scaffolds   = [r['scaffold']   for r in large]

        # ESD tests
        gene_esd   = rosner_esd(gene_vals,   R_MAX, ALPHA, side='lower')
        repeat_esd = rosner_esd(repeat_vals, R_MAX, ALPHA, side='upper')

        gene_esd_map   = {res['orig_idx']: res for res in gene_esd}
        repeat_esd_map = {res['orig_idx']: res for res in repeat_esd}

        # Modified Z-scores
        gene_modz   = modified_zscore(gene_vals,   side='lower')   # negative = depleted
        repeat_modz = modified_zscore(repeat_vals, side='upper')   # positive = enriched

        def fmt(val, decimals=4):
            if val is None or (isinstance(val, float) and np.isnan(val)):
                return ''
            return round(float(val), decimals)

        # Collect outlier names for summary
        gene_confident_names   = []
        repeat_confident_names = []

        for i, row in enumerate(large):
            g  = gene_esd_map.get(i, {})
            rp = repeat_esd_map.get(i, {})

            gene_esd_out    = g.get('is_outlier', False)
            repeat_esd_out  = rp.get('is_outlier', False)
            gene_modz_out   = gene_modz[i]   < -MOD_Z_THRESH
            repeat_modz_out = repeat_modz[i] >  MOD_Z_THRESH

            gene_confident   = gene_esd_out   and gene_modz_out
            repeat_confident = repeat_esd_out and repeat_modz_out

            if gene_confident:
                gene_confident_names.append(row['scaffold'])
            if repeat_confident:
                repeat_confident_names.append(row['scaffold'])

            detail_rows.append({
                'species':          species,
                'scaffold':         row['scaffold'],
                'size_mb':          round(row['size_bp'] / 1e6, 3),
                'n_genes':          row['n_genes'],
                'gene_pct':         row['gene_pct'],
                'repeat_pct':       row['repeat_pct'],
                # ESD columns — gene
                'gene_esd_step':    g.get('step', ''),
                'gene_Ri':          fmt(g.get('Ri')),
                'gene_lambda':      fmt(g.get('lambda_i')),
                'gene_esd_outlier': gene_esd_out,
                # Modified Z — gene
                'gene_modZ':        fmt(gene_modz[i], 3),
                'gene_confident':   gene_confident,
                # ESD columns — repeat
                'repeat_esd_step':    rp.get('step', ''),
                'repeat_Ri':          fmt(rp.get('Ri')),
                'repeat_lambda':      fmt(rp.get('lambda_i')),
                'repeat_esd_outlier': repeat_esd_out,
                # Modified Z — repeat
                'repeat_modZ':        fmt(repeat_modz[i], 3),
                'repeat_confident':   repeat_confident,
            })

        both_confident = sorted(set(gene_confident_names) & set(repeat_confident_names))

        summary_rows.append({
            'species':            species,
            'n_sequences_tested': n,
            'gene_confident':     ', '.join(gene_confident_names)   or 'none',
            'repeat_confident':   ', '.join(repeat_confident_names) or 'none',
            'both_confident':     ', '.join(both_confident)         or 'none',
        })

    # --- Write outputs ---
    detail_path  = os.path.join(BASE, "esd_outlier_results.tsv")
    summary_path = os.path.join(BASE, "esd_outlier_summary.tsv")

    detail_fields = [
        'species', 'scaffold', 'size_mb', 'n_genes', 'gene_pct', 'repeat_pct',
        'gene_esd_step', 'gene_Ri', 'gene_lambda', 'gene_esd_outlier', 'gene_modZ', 'gene_confident',
        'repeat_esd_step', 'repeat_Ri', 'repeat_lambda', 'repeat_esd_outlier', 'repeat_modZ', 'repeat_confident',
    ]
    with open(detail_path, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=detail_fields, delimiter='\t')
        w.writeheader()
        w.writerows(detail_rows)

    summary_fields = ['species', 'n_sequences_tested',
                      'gene_confident', 'repeat_confident', 'both_confident']
    with open(summary_path, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=summary_fields, delimiter='\t')
        w.writeheader()
        w.writerows(summary_rows)

    # --- Print summary ---
    print(f"\n{'Species':<42} {'N':>3}  {'Gene-depleted (confident)':<45} {'Repeat-enriched (confident)':<45} {'Both'}")
    print("─" * 165)
    for row in summary_rows:
        print(f"{row['species']:<42} {row['n_sequences_tested']:>3}  "
              f"{row['gene_confident']:<45} {row['repeat_confident']:<45} {row['both_confident']}")

    print(f"\nDetailed results : {detail_path}")
    print(f"Summary          : {summary_path}")
    print(f"\nNote: 'confident' = significant by both Rosner ESD (alpha={ALPHA}, r={R_MAX})")
    print(f"      AND |modified Z-score| > {MOD_Z_THRESH} (Iglewicz & Hoaglin 1993)")


if __name__ == '__main__':
    main()
