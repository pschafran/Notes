#!/usr/bin/env python3
"""
analyze_genome_composition.py
─────────────────────────────
Compute and visualize gene and repeat content per chromosome/scaffold.

Inputs
------
  -g / --genes    Gene annotation (GFF3 or GTF). Handles:
                    • GFF3 with 'gene' feature rows
                    • GFF3 with only 'mRNA' rows (groups by geneID= attribute)
                    • GTF with 'gene' or 'transcript' rows
  -r / --repeats  Repeat annotation (GFF3, e.g. EDTA output).
                  Classification read from feature column (e.g. LTR/Gypsy)
                  or from a Classification= / Class= attribute.
  -f / --fai      FASTA index (.fai) for authoritative scaffold sizes.
                  If omitted, sizes are inferred from GFF coordinates (less
                  accurate — may underestimate scaffold length).

Outputs
-------
  <prefix>_composition_stats.tsv   Per-scaffold counts and percentages
  <prefix>_composition.pdf/.png    Stacked bar chart

Usage examples
--------------
  # Basic
  analyze_genome_composition.py -g genes.gff3 -r repeats.gff -f genome.fa.fai

  # Only selected scaffolds, custom output location
  analyze_genome_composition.py -g genes.gff3 -r repeats.gff -f genome.fa.fai \\
      -s chr1,chr2,chrU,chrV -o results/ -p MySpecies

  # Add pararetrovirus class, rename prefix
  analyze_genome_composition.py -g genes.gtf -r repeats.gff -f genome.fa.fai \\
      --extra-classes pararetrovirus -p Papro
"""

import argparse
import os
import re
import sys
from collections import defaultdict

import matplotlib
matplotlib.use('Agg')
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np


# ── Repeat classification ─────────────────────────────────────────────────────

DEFAULT_MAJOR_ORDER = ['LTR', 'DNA', 'LINE', 'MITE', 'Unknown']

# Map non-standard base names → canonical major class
BASE_REMAP = {
    'TIR':        'DNA',
    'polinton':   'DNA',
    'Penelope':   'LINE',
    'RC':         'DNA',
    'SINE':       'LINE',
    'Helitron':   'DNA',
    'Gypsy_LTR_retrotransposon': 'LTR',
    'hAT_TIR_transposon': 'DNA',
    'I_LINE_retrotransposon': 'LINE',
    'IS3EU_TIR_transposon': 'DNA',
    'L1_LINE_retrotransposon': 'LINE',
    'L2_LINE_retrotransposon': 'LINE',
    'long_terminal_repeat': 'LTR',
    'LTR_retrotransposon': 'LTR',
    'Mutator_TIR_transposon': 'DNA',
    'Penelope_retrotransposon': 'LINE',
    'PIF_Harbinger_TIR_transposon': 'DNA',
    'RTE_LINE_retrotransposon': 'LINE',
    'Sola_TIR_transposon': 'DNA',
    'Tad1_LINE_retrotransposon': 'LINE',
    'Tc1_Mariner_TIR_transposon': 'DNA',
    '5S_SINE_retrotransposon': 'LINE',
    'CACTA_TIR_transposon': 'DNA',
    'Copia_LTR_retrotransposon': 'LTR',
    'CRE_LINE_retrotransposon': 'LINE'
}


def get_major_class(raw_class, major_order, remap):
    """Return canonical major repeat class from a raw classification string."""
    base = raw_class.split('/')[0].strip()
    base = remap.get(base, base)
    return base if base in major_order else 'Unknown'


# ── Interval utilities ────────────────────────────────────────────────────────

def merge_intervals(intervals):
    """Merge overlapping (start, end) intervals. Coordinates are 0-based half-open."""
    if not intervals:
        return []
    srt = sorted(intervals)
    merged = [list(srt[0])]
    for s, e in srt[1:]:
        if s < merged[-1][1]:
            merged[-1][1] = max(merged[-1][1], e)
        else:
            merged.append([s, e])
    return [tuple(x) for x in merged]


def count_bp(intervals):
    return sum(e - s for s, e in intervals)


def intersect_bp(a_intervals, b_intervals):
    """Total bp that overlap between two sorted, merged interval lists."""
    total = 0
    j = 0
    for as_, ae in a_intervals:
        while j < len(b_intervals) and b_intervals[j][1] <= as_:
            j += 1
        k = j
        while k < len(b_intervals) and b_intervals[k][0] < ae:
            total += min(ae, b_intervals[k][1]) - max(as_, b_intervals[k][0])
            k += 1
    return total


# ── GFF / GTF attribute parsing ───────────────────────────────────────────────

def parse_attrs_gff3(attr_str):
    """Parse GFF3 key=value;... attributes."""
    attrs = {}
    for field in attr_str.split(';'):
        field = field.strip()
        if '=' in field:
            k, v = field.split('=', 1)
            attrs[k.strip()] = v.strip()
    return attrs


def parse_attrs_gtf(attr_str):
    """Parse GTF key "value"; ... attributes."""
    attrs = {}
    for m in re.finditer(r'(\w+)\s+"([^"]*)"', attr_str):
        attrs[m.group(1)] = m.group(2)
    return attrs


def detect_format(path):
    """Guess GFF3 vs GTF from file extension."""
    return 'gtf' if path.endswith('.gtf') else 'gff3'


# ── Gene loading ──────────────────────────────────────────────────────────────

def load_genes(path, scaffolds):
    """
    Load gene loci and exon intervals (0-based half-open coordinates).

    Gene loci (for n_genes count) use the full gene span from 'gene' / 'transcript'
    / 'mRNA' feature rows.  Exon intervals (for gene_bp calculation) use 'exon'
    feature rows so that intronic sequence is excluded.  If no 'exon' features are
    present, gene spans are used as a fallback (with a warning).

    Returns
    -------
    gene_loci     : {scaffold: [(start, end), ...]}  merged gene spans
    exon_intervals: {scaffold: [(start, end), ...]}  merged exon intervals
    target        : feature type used for gene counting
    """
    fmt = detect_format(path)
    parse_attrs = parse_attrs_gtf if fmt == 'gtf' else parse_attrs_gff3

    # First pass: collect feature types present
    feature_counts = defaultdict(int)
    with open(path) as fh:
        for line in fh:
            if line.startswith('#') or not line.strip():
                continue
            p = line.split('\t')
            if len(p) < 3:
                continue
            if p[0] not in scaffolds:
                continue
            feature_counts[p[2]] += 1

    # Choose target feature for gene counting
    for target in ('gene', 'transcript', 'mRNA'):
        if feature_counts.get(target, 0) > 0:
            break
    else:
        sys.exit(f'ERROR: no gene/transcript/mRNA features found in {path} '
                 f'for the specified scaffolds.\n'
                 f'  Features seen: {dict(feature_counts)}\n'
                 f'  Scaffolds queried: {sorted(scaffolds)[:10]}')

    has_exons = feature_counts.get('exon', 0) > 0
    if not has_exons:
        print('  WARNING: no exon features found — gene_bp will include introns '
              '(falling back to full gene span)')

    loci     = defaultdict(dict)   # {scaffold: {gene_id: [min_start, max_end]}}
    exon_raw = defaultdict(list)   # {scaffold: [(start, end), ...]}

    with open(path) as fh:
        for line in fh:
            if line.startswith('#') or not line.strip():
                continue
            p = line.strip().split('\t')
            if len(p) < 9:
                continue
            scaf = p[0]
            if scaf not in scaffolds:
                continue
            feat  = p[2]
            start = int(p[3]) - 1   # GFF is 1-based; convert to 0-based
            end   = int(p[4])

            if feat == target:
                attrs = parse_attrs(p[8])
                # Derive gene ID
                if fmt == 'gtf':
                    gene_id = attrs.get('gene_id', attrs.get('transcript_id', f'{scaf}:{start}'))
                else:
                    gene_id = (attrs.get('geneID') or attrs.get('gene_id') or
                               attrs.get('ID', f'{scaf}:{start}').split('.')[0])
                b = loci[scaf].setdefault(gene_id, [start, end])
                b[0] = min(b[0], start)
                b[1] = max(b[1], end)

            elif feat == 'exon' and has_exons:
                exon_raw[scaf].append((start, end))

    # Merged gene spans (used for n_genes count)
    gene_loci = {scaf: merge_intervals([(b[0], b[1]) for b in genes.values()])
                 for scaf, genes in loci.items()}

    # Merged exon intervals (used for gene_bp; excludes intronic sequence)
    if has_exons:
        exon_intervals = {scaf: merge_intervals(ivs) for scaf, ivs in exon_raw.items()}
    else:
        exon_intervals = gene_loci   # fallback: same as gene spans

    n_total  = sum(len(v) for v in loci.values())
    n_merged = sum(len(v) for v in gene_loci.values())
    print(f'  Genes loaded: {n_total} loci → {n_merged} merged on '
          f'{len(gene_loci)} scaffolds  [feature: {target}]')
    if has_exons:
        n_exon_ivs = sum(len(v) for v in exon_intervals.values())
        n_exon_bp  = sum(count_bp(v) for v in exon_intervals.values())
        print(f'  Exon intervals: {n_exon_ivs} merged intervals, {n_exon_bp:,} bp '
              f'(introns excluded from gene_bp)')

    missing = [s for s in scaffolds if s not in gene_loci]
    if missing:
        print(f'  WARNING: no genes found on: {missing}')

    return gene_loci, exon_intervals, target


# ── Repeat loading ────────────────────────────────────────────────────────────

def load_repeats(path, scaffolds, major_order, remap):
    """
    Load repeat intervals as {scaffold: {major_class: [(start, end), ...]}}.

    Classification source (tried in order):
      1. feature column (col 2): e.g. 'LTR/Gypsy'
      2. Classification= attribute
      3. Class= attribute
      4. Name= attribute (first token before '/')
    """
    data = defaultdict(lambda: defaultdict(list))
    n_total = 0
    n_unknown = 0

    with open(path) as fh:
        for line in fh:
            if line.startswith('#') or not line.strip():
                continue
            p = line.strip().split('\t')
            if len(p) < 9:
                continue
            scaf = p[0]
            if scaf not in scaffolds:
                continue

            start = int(p[3]) - 1
            end   = int(p[4])
            raw   = p[2].strip()

            # Try feature column first; fall back to attributes
            if '/' in raw or raw in major_order or raw in remap:
                cls = raw
            else:
                attrs = parse_attrs_gff3(p[8])
                cls = (attrs.get('Classification') or
                       attrs.get('Class') or
                       attrs.get('class') or
                       attrs.get('Name', 'Unknown').split('#')[-1])

            major = get_major_class(cls, major_order, remap)
            data[scaf][major].append((start, end))
            n_total += 1
            if major == 'Unknown':
                n_unknown += 1

    # Merge overlapping intervals per scaffold per class
    result = {}
    for scaf in data:
        result[scaf] = {cls: merge_intervals(ivs)
                        for cls, ivs in data[scaf].items()}

    print(f'  Repeats loaded: {n_total} features on {len(result)} scaffolds '
          f'({n_unknown} → Unknown)')

    missing = [s for s in scaffolds if s not in result]
    if missing:
        print(f'  WARNING: no repeats found on: {missing}')

    return result


# ── Per-scaffold statistics ───────────────────────────────────────────────────

def compute_stats(scaffolds, sizes, gene_loci, exon_intervals, repeats, major_order):
    """
    For each scaffold compute:
      gene_bp, repeat_bp (by class), gene_repeat_overlap_bp, unannotated_bp.

    gene_bp is computed from merged exon intervals (introns excluded).
    n_genes is the count of merged gene-span loci.

    Fraction breakdown (mutually exclusive, sums to scaffold size):
      gene_only    = exon_bp - overlap_bp
      repeat_only  = repeat_bp_total - overlap_bp
      overlap      = exon ∩ repeat
      unannotated  = size - gene_only - repeat_only - overlap
    """
    rows = []
    for scaf in scaffolds:
        size = sizes.get(scaf, 0)
        g_ivs = gene_loci.get(scaf, [])       # gene spans → for n_genes
        e_ivs = exon_intervals.get(scaf, [])  # exon intervals → for gene_bp
        r_by_class = repeats.get(scaf, {})
        all_r_ivs  = merge_intervals(
            [iv for ivs in r_by_class.values() for iv in ivs])

        gene_bp   = count_bp(e_ivs)
        repeat_bp = count_bp(all_r_ivs)
        overlap   = intersect_bp(e_ivs, all_r_ivs)

        gene_only  = gene_bp - overlap
        rep_only   = repeat_bp - overlap
        unanno     = max(0, size - gene_only - rep_only - overlap)

        class_bp = {cls: count_bp(r_by_class.get(cls, [])) for cls in major_order}

        n_genes = len(g_ivs)

        rows.append({
            'scaffold':       scaf,
            'size_bp':        size,
            'n_genes':        n_genes,
            'gene_bp':        gene_bp,
            'gene_pct':       100 * gene_bp / size if size else 0,
            'repeat_bp':      repeat_bp,
            'repeat_pct':     100 * repeat_bp / size if size else 0,
            'gene_only_bp':   gene_only,
            'gene_only_pct':  100 * gene_only / size if size else 0,
            'overlap_bp':     overlap,
            'overlap_pct':    100 * overlap / size if size else 0,
            'unannotated_bp': unanno,
            'unanno_pct':     100 * unanno / size if size else 0,
            **{f'{cls}_bp':  class_bp[cls]               for cls in major_order},
            **{f'{cls}_pct': 100 * class_bp[cls] / size if size else 0
               for cls in major_order},
        })
    return rows


# ── Stats TSV ─────────────────────────────────────────────────────────────────

def write_stats(rows, path, major_order):
    fixed = ['scaffold', 'size_bp', 'n_genes',
             'gene_bp', 'gene_pct', 'repeat_bp', 'repeat_pct',
             'gene_only_bp', 'gene_only_pct',
             'overlap_bp', 'overlap_pct',
             'unannotated_bp', 'unanno_pct']
    class_cols = [f'{c}_{s}' for c in major_order for s in ('bp', 'pct')]
    fieldnames = fixed + class_cols

    with open(path, 'w', newline='') as fh:
        import csv
        w = csv.DictWriter(fh, fieldnames=fieldnames, delimiter='\t',
                           extrasaction='ignore')
        w.writeheader()
        for row in rows:
            w.writerow({k: (f'{v:.1f}' if isinstance(v, float) else v)
                        for k, v in row.items()})
    print(f'  Stats written: {path}')


# ── Composition plot ──────────────────────────────────────────────────────────

# Default colour palette (muted pastels matching project style)
DEFAULT_COLORS = {
    'LTR':          '#c8a882',
    'DNA':          '#a8b8c8',
    'LINE':         '#b8a8c8',
    'MITE':         '#c8c8a0',
    'pararetrovirus': '#d4a0a0',
    'Unknown':      '#c8c8c8',
    'gene_only':    '#8aab7e',
    'overlap':      '#6a8b5e',
    'unannotated':  '#e8e8e8',
}


def plot_composition(rows, major_order, outbase, colors=None):
    if colors is None:
        colors = DEFAULT_COLORS

    plt.rcParams.update({
        'font.family':       'Liberation Sans',
        'pdf.fonttype':      42,
        'ps.fonttype':       42,
        'axes.spines.top':   False,
        'axes.spines.right': False,
        'font.size':         10,
    })

    scaffolds = [r['scaffold'] for r in rows]
    n = len(scaffolds)
    x = np.arange(n)
    width = 0.65

    # Stack order: gene_only, overlap, repeat classes (LTR first), unannotated
    stack = ['gene_only', 'overlap'] + major_order + ['unannotated']
    stack_label = {
        'gene_only':    'Gene (non-repeat)',
        'overlap':      'Gene + repeat overlap',
        'unannotated':  'Unannotated',
        **{cls: cls for cls in major_order},
    }

    def get_pct(row, key):
        if key == 'gene_only':
            return row['gene_only_pct']
        if key == 'overlap':
            return row['overlap_pct']
        if key == 'unannotated':
            return row['unanno_pct']
        return row.get(f'{key}_pct', 0)

    fig_w = max(7, n * 1.2 + 3)
    fig, ax = plt.subplots(figsize=(fig_w, 6))

    bottoms = np.zeros(n)
    patches = []
    for key in stack:
        vals = np.array([get_pct(r, key) for r in rows])
        col  = colors.get(key, '#cccccc')
        bars = ax.bar(x, vals, bottom=bottoms, width=width,
                      color=col, edgecolor='white', linewidth=0.4)
        bottoms += vals
        patches.append(mpatches.Patch(color=col, label=stack_label[key]))

    # Scaffold size labels below bars
    for xi, row in zip(x, rows):
        size_mb = row['size_bp'] / 1e6
        ax.text(xi, -2.5, f'{size_mb:.1f} Mb', ha='center', va='top',
                fontsize=8, color='#555555')

    ax.set_xticks(x)
    ax.set_xticklabels(scaffolds, rotation=30, ha='right', fontsize=9)
    ax.set_ylabel('Genome composition (%)', fontsize=10)
    ax.set_ylim(0, 105)
    ax.set_xlim(-0.6, n - 0.4)

    # Legend (reverse so gene is on top)
    ax.legend(handles=patches[::-1], loc='upper right', frameon=False,
              fontsize=8.5, ncol=1)

    plt.tight_layout()
    plt.subplots_adjust(bottom=0.18)

    for ext in ('pdf', 'png'):
        p = f'{outbase}.{ext}'
        fig.savefig(p, dpi=150, bbox_inches='tight')
        print(f'  Figure saved: {p}')
    plt.close(fig)


# ── CLI ───────────────────────────────────────────────────────────────────────

def parse_args():
    p = argparse.ArgumentParser(
        description='Compute and plot gene + repeat composition per scaffold.',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__)
    p.add_argument('-g', '--genes',   required=True,
                   help='Gene annotation file (GFF3 or GTF)')
    p.add_argument('-r', '--repeats', required=True,
                   help='Repeat annotation file (GFF3, e.g. EDTA output)')
    p.add_argument('-f', '--fai',     default=None,
                   help='FASTA index (.fai) for scaffold sizes. '
                        'If omitted, sizes are inferred from GFF max coordinates.')
    p.add_argument('-s', '--scaffolds', default=None,
                   help='Comma-separated scaffold names to include. '
                        'If omitted, all scaffolds present in the gene file are used.')
    p.add_argument('-o', '--outdir',  default='.',
                   help='Output directory (default: current directory)')
    p.add_argument('-p', '--prefix',  default='composition',
                   help='Output filename prefix (default: composition)')
    p.add_argument('--extra-classes', default=None,
                   help='Comma-separated additional repeat classes to track '
                        '(e.g. pararetrovirus). Added before Unknown in stack order.')
    p.add_argument('--no-plot', action='store_true',
                   help='Write stats TSV only, skip figure')
    return p.parse_args()


def main():
    args = parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    # Build major class order
    major_order = list(DEFAULT_MAJOR_ORDER)
    if args.extra_classes:
        extras = [c.strip() for c in args.extra_classes.split(',') if c.strip()]
        # Insert before 'Unknown'
        idx = major_order.index('Unknown')
        for e in extras:
            if e not in major_order:
                major_order.insert(idx, e)
                idx += 1
    remap = dict(BASE_REMAP)

    print(f'Repeat major classes: {major_order}')

    # ── Scaffold list ──────────────────────────────────────────────────────────
    if args.scaffolds:
        scaffolds = [s.strip() for s in args.scaffolds.split(',')]
    else:
        # Infer from gene file
        scaffolds_seen = set()
        with open(args.genes) as fh:
            for line in fh:
                if line.startswith('#') or not line.strip():
                    continue
                p = line.split('\t')
                if len(p) >= 1:
                    scaffolds_seen.add(p[0])
        scaffolds = sorted(scaffolds_seen)
        print(f'  Scaffolds inferred from gene file: {scaffolds}')

    scaffolds_set = set(scaffolds)

    # ── Scaffold sizes ─────────────────────────────────────────────────────────
    sizes = {}
    if args.fai:
        with open(args.fai) as fh:
            for line in fh:
                p = line.strip().split('\t')
                if p[0] in scaffolds_set:
                    sizes[p[0]] = int(p[1])
        missing_sizes = [s for s in scaffolds if s not in sizes]
        if missing_sizes:
            print(f'  WARNING: no size in .fai for: {missing_sizes}')
    else:
        print('  No .fai provided — scaffold sizes will be inferred from GFF max '
              'coordinates (may underestimate true scaffold length).')

    # ── Load data ──────────────────────────────────────────────────────────────
    print('\nLoading gene annotations...')
    genes, exon_intervals, gene_feature = load_genes(args.genes, scaffolds_set)

    print('\nLoading repeat annotations...')
    repeats = load_repeats(args.repeats, scaffolds_set, major_order, remap)

    # Fill missing sizes from GFF max coordinates
    for scaf in scaffolds:
        if scaf not in sizes:
            g_max = max((e for _, e in genes.get(scaf, [])), default=0)
            r_max = max((e for ivs in repeats.get(scaf, {}).values()
                         for _, e in ivs), default=0)
            sizes[scaf] = max(g_max, r_max)
            if sizes[scaf] == 0:
                print(f'  WARNING: cannot determine size for {scaf} '
                      f'(no genes or repeats found) — set to 0')

    # ── Compute stats ──────────────────────────────────────────────────────────
    print('\nComputing per-scaffold statistics...')
    rows = compute_stats(scaffolds, sizes, genes, exon_intervals, repeats, major_order)

    # Print summary to stdout
    print(f'\n{"Scaffold":<20} {"Size (Mb)":>10} {"Genes":>7} '
          f'{"Gene%":>7} {"Repeat%":>8} {"Unannot%":>9}')
    print('-' * 65)
    for r in rows:
        print(f'{r["scaffold"]:<20} {r["size_bp"]/1e6:>10.2f} '
              f'{r["n_genes"]:>7} {r["gene_pct"]:>7.1f} '
              f'{r["repeat_pct"]:>8.1f} {r["unanno_pct"]:>9.1f}')

    # ── Write outputs ──────────────────────────────────────────────────────────
    print()
    tsv_path = os.path.join(args.outdir, f'{args.prefix}_composition_stats.tsv')
    write_stats(rows, tsv_path, major_order)

    if not args.no_plot:
        fig_base = os.path.join(args.outdir, f'{args.prefix}_composition')
        plot_composition(rows, major_order, fig_base)

    print('\nDone.')


if __name__ == '__main__':
    main()
