"""
Coverage and genomic context analysis of per-site 5mC methylation categories
in PaproF and PaproM autosomes.

Methylation categories:
  unmethylated : pct_modified <  20%
  intermediate : 20% <= pct_modified < 80%
  methylated   : pct_modified >= 80%

Analyses:
  1. Coverage (Nvalid_cov) distribution by category
  2. Genomic context (gene body, repeat class, intergenic) by category
"""

import bisect
from collections import defaultdict

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.patches as mpatches

# ── Paths ───────────────────────────────────────────────────────────────────
BASE   = '/media/data/projects/hornwort_sex_chromosomes/analysis/Paraphymatoceros_proskaueri'
BASE_F = f'{BASE}/Papro252-9.2.2'
BASE_M = f'{BASE}/Papro252-3'
OUTDIR = f'{BASE}/sex_chromosome_analyses'

METHYL_F = f'{BASE_F}/methylation/PaproF_genome.5mCG_5hmCG.bed'
METHYL_M = f'{BASE_M}/methylation/PaproM_genome.5mCG_5hmCG.bed'
GFF_F    = f'{BASE_F}/final_genome_prep/pilon.3_renamed.fasta.mod.EDTA.TEanno.gff3'
GTF_F    = f'{BASE_F}/final_genome_prep/braker_renamed.gtf'
GFF_M    = f'{BASE_M}/final_genome_prep/PaproM.renamed_renamed.fasta.mod.EDTA.TEanno.gff3'
GTF_M    = f'{BASE_M}/final_genome_prep/braker_renamed.gtf'

F_AUTO = ['PaproF.S1', 'PaproF.S2', 'PaproF.S3', 'PaproF.S4']
M_AUTO = ['PaproM.S1', 'PaproM.S2', 'PaproM.S3', 'PaproM.S4']
MIN_COV = 5

MAJOR_ORDER = ['LTR', 'DNA', 'LINE', 'MITE', 'pararetrovirus', 'Unknown']
CATS = ['unmethylated', 'intermediate', 'methylated']

# ── Figure style ─────────────────────────────────────────────────────────────
plt.rcParams.update({
    'font.family': 'Liberation Sans', 'font.size': 8,
    'axes.titlesize': 9, 'axes.labelsize': 8,
    'xtick.labelsize': 7.5, 'ytick.labelsize': 7.5,
    'axes.spines.top': False, 'axes.spines.right': False,
    'axes.linewidth': 0.6, 'xtick.major.width': 0.6, 'ytick.major.width': 0.6,
    'xtick.major.size': 3, 'ytick.major.size': 3,
    'pdf.fonttype': 42, 'ps.fonttype': 42,
    'figure.facecolor': 'white', 'axes.facecolor': 'white',
})

# Colours for methylation categories
CAT_COL = {
    'unmethylated': '#6a9ac4',
    'intermediate': '#c4a86a',
    'methylated':   '#c46a6a',
}

# Context colours
CTX_COLORS = {
    'intergenic':   '#d8d8d8',
    'gene_only':    '#7aabcf',
    'gene+LTR':     '#9a7a7a',
    'gene+DNA':     '#7aab8a',
    'gene+LINE':    '#8a8ab8',
    'gene+MITE':    '#c8a86e',
    'gene+para':    '#9aacb8',
    'gene+Unknown': '#7aaabf',
    'LTR':          '#c47a7a',
    'DNA':          '#5a9070',
    'LINE':         '#6a6a98',
    'MITE':         '#b08848',
    'pararetrovirus': '#7a9ca8',
    'Unknown':      '#5a8a9f',
}

# ── Helpers ──────────────────────────────────────────────────────────────────

def get_major(cls):
    major = cls.split('/')[0] if '/' in cls else cls
    if major == 'TIR':     return 'DNA'
    if major == 'Penelope': return 'LINE'
    return major


def load_meth(bed_file, scaffolds, min_cov=MIN_COV):
    """Return dict scaffold -> (pos, cov, pct) all sorted by pos."""
    raw = {s: ([], [], []) for s in scaffolds}
    scset = set(scaffolds)
    with open(bed_file) as fh:
        for line in fh:
            p = line.split('\t')
            if len(p) < 11 or p[3] != 'm':
                continue
            if p[0] not in scset:
                continue
            cov = int(p[9])
            if cov < min_cov:
                continue
            raw[p[0]][0].append(int(p[1]))
            raw[p[0]][1].append(cov)
            raw[p[0]][2].append(float(p[10]))
    out = {}
    for s in scaffolds:
        pos = np.array(raw[s][0], dtype=np.int64)
        cov = np.array(raw[s][1], dtype=np.int32)
        pct = np.array(raw[s][2], dtype=np.float32)
        idx = np.argsort(pos)
        out[s] = (pos[idx], cov[idx], pct[idx])
    return out


def load_genes(gtf_file, scaffolds):
    """Return dict scaffold -> sorted list of (start, end) gene-body intervals."""
    scaf_set = set(scaffolds)
    genes = defaultdict(list)
    with open(gtf_file) as fh:
        for line in fh:
            if line.startswith('#'):
                continue
            p = line.rstrip().split('\t')
            if len(p) < 5 or p[2] != 'gene' or p[0] not in scaf_set:
                continue
            genes[p[0]].append((int(p[3]), int(p[4])))
    for s in genes:
        genes[s].sort()
    return dict(genes)


def load_repeats(gff_file, scaffolds):
    """Return dict scaffold -> {major_class: sorted [(start, end)]}."""
    scaf_set = set(scaffolds)
    reps = defaultdict(lambda: defaultdict(list))
    with open(gff_file) as fh:
        for line in fh:
            if line.startswith('#'):
                continue
            p = line.rstrip().split('\t')
            if len(p) < 9 or p[0] not in scaf_set:
                continue
            attrs = {k: v for a in p[8].split(';') if '=' in a
                     for k, v in [a.split('=', 1)]}
            cls   = attrs.get('Classification', 'Unknown')
            major = get_major(cls)
            reps[p[0]][major].append((int(p[3]), int(p[4])))
    for s in reps:
        for maj in reps[s]:
            reps[s][maj].sort()
    return dict(reps)


def label_sites(pos_arr, gene_ivs, rep_ivs_by_major):
    """
    For each site in pos_arr (sorted), determine:
      is_gene : bool array  — site falls within a gene body
      rep_cls : str array   — major repeat class, or '' if not in a repeat
    Uses bisect for O(n_intervals * avg_width) total work.
    """
    n = len(pos_arr)
    is_gene = np.zeros(n, dtype=bool)
    rep_cls = np.full(n, '', dtype=object)

    for major, ivs in rep_ivs_by_major.items():
        for s, e in ivs:
            lo = bisect.bisect_left(pos_arr, s)
            hi = bisect.bisect_right(pos_arr, e)
            if lo < hi:
                rep_cls[lo:hi] = major   # last-write wins for overlapping intervals

    for s, e in gene_ivs:
        lo = bisect.bisect_left(pos_arr, s)
        hi = bisect.bisect_right(pos_arr, e)
        if lo < hi:
            is_gene[lo:hi] = True

    return is_gene, rep_cls


def categorize(pct_arr):
    cats = np.empty(len(pct_arr), dtype=object)
    cats[pct_arr < 20]                           = 'unmethylated'
    cats[(pct_arr >= 20) & (pct_arr < 80)]       = 'intermediate'
    cats[pct_arr >= 80]                           = 'methylated'
    return cats


# ── Context breakdown helper ─────────────────────────────────────────────────

def context_fracs(mask, is_gene, rep_cls):
    """
    Given a boolean mask selecting sites, return dict of context -> fraction.
    Contexts (mutually exclusive):
      intergenic | gene_only | gene+<major> | <major>   (repeat-only)
    """
    in_rep = rep_cls != ''
    sub_gene = is_gene[mask]
    sub_rep  = in_rep[mask]
    sub_cls  = rep_cls[mask]
    n = mask.sum()
    if n == 0:
        return {}
    out = {}
    # intergenic
    out['intergenic'] = float((~sub_gene & ~sub_rep).sum()) / n
    # gene only (no repeat)
    out['gene_only']  = float((sub_gene & ~sub_rep).sum()) / n
    # gene + repeat (by major class)
    for maj in MAJOR_ORDER:
        k = (sub_gene & (sub_cls == maj)).sum()
        label = 'gene+' + ('para' if maj == 'pararetrovirus' else maj)
        out[label] = float(k) / n
    # repeat only (by major class)
    for maj in MAJOR_ORDER:
        k = (~sub_gene & (sub_cls == maj)).sum()
        out[maj] = float(k) / n
    return out


# ══════════════════════════════════════════════════════════════════════════════
# Main loop
# ══════════════════════════════════════════════════════════════════════════════

all_results = {}

for label, scaffolds, methyl_file, gff_file, gtf_file in [
    ('Female (PaproF)', F_AUTO, METHYL_F, GFF_F, GTF_F),
    ('Male (PaproM)',   M_AUTO, METHYL_M, GFF_M, GTF_M),
]:
    print(f'\n{"="*65}')
    print(f'  {label}')
    print(f'{"="*65}')

    print('  Loading methylation...')
    meth = load_meth(methyl_file, scaffolds)
    print('  Loading genes...')
    genes = load_genes(gtf_file, scaffolds)
    print('  Loading repeats...')
    reps  = load_repeats(gff_file, scaffolds)

    # Concatenate across scaffolds
    pos_list, cov_list, pct_list, gene_list, rep_list = [], [], [], [], []
    for s in scaffolds:
        pos, cov, pct = meth[s]
        is_gene, rep_cls = label_sites(pos, genes.get(s, []),
                                        reps.get(s, {}))
        pos_list.append(pos);   cov_list.append(cov)
        pct_list.append(pct);   gene_list.append(is_gene)
        rep_list.append(rep_cls)

    all_pos  = np.concatenate(pos_list)
    all_cov  = np.concatenate(cov_list)
    all_pct  = np.concatenate(pct_list)
    all_gene = np.concatenate(gene_list)
    all_rep  = np.concatenate(rep_list)
    all_cats = categorize(all_pct)

    # ── 1. Coverage by category ───────────────────────────────────────────
    print(f'\n  Coverage (Nvalid_cov) by methylation category:')
    print(f'  {"Category":16s}  {"N":>10}  {"Median":>7}  {"Mean":>7}  '
          f'{"P10":>6}  {"P25":>6}  {"P75":>6}  {"P90":>6}  {"P99":>6}')
    cov_data = {}
    for cat in CATS:
        mask = all_cats == cat
        c = all_cov[mask]
        if len(c) == 0:
            print(f'  {cat:16s}  {"0":>10}')
            continue
        p10, p25, p75, p90, p99 = np.percentile(c, [10, 25, 75, 90, 99])
        print(f'  {cat:16s}  {len(c):>10,}  {np.median(c):>7.1f}  {np.mean(c):>7.1f}  '
              f'{p10:>6.1f}  {p25:>6.1f}  {p75:>6.1f}  {p90:>6.1f}  {p99:>6.1f}')
        cov_data[cat] = c

    # ── 2. Context by category ───────────────────────────────────────────
    in_rep = all_rep != ''
    print(f'\n  Context breakdown by methylation category:')
    ctx_cols = (['intergenic', 'gene_only'] +
                ['gene+' + ('para' if m == 'pararetrovirus' else m) for m in MAJOR_ORDER] +
                list(MAJOR_ORDER))
    hdr = f'  {"Category":16s}  {"N":>10}  ' + '  '.join(f'{c[:12]:>13}' for c in ctx_cols)
    print(hdr)
    ctx_data = {}
    for cat in CATS:
        mask = all_cats == cat
        cf = context_fracs(mask, all_gene, all_rep)
        ctx_data[cat] = cf
        vals = [cf.get(c, 0.0) * 100 for c in ctx_cols]
        row = (f'  {cat:16s}  {mask.sum():>10,}  ' +
               '  '.join(f'{v:>12.1f}%' for v in vals))
        print(row)

    # ── Condensed summary ────────────────────────────────────────────────
    print(f'\n  Condensed summary:')
    print(f'  {"Category":16s}  {"intergenic":>11}  {"gene body":>10}  '
          f'{"in repeat":>10}  {"gene+repeat":>12}')
    for cat in CATS:
        mask = all_cats == cat
        n = mask.sum()
        sub_g = all_gene[mask]
        sub_r = in_rep[mask]
        frac_interg = (~sub_g & ~sub_r).sum() / n * 100
        frac_gene   = sub_g.sum() / n * 100
        frac_rep    = sub_r.sum() / n * 100
        frac_both   = (sub_g & sub_r).sum() / n * 100
        print(f'  {cat:16s}  {frac_interg:>10.1f}%  {frac_gene:>9.1f}%  '
              f'{frac_rep:>9.1f}%  {frac_both:>11.1f}%')

    all_results[label] = {
        'cov_data': cov_data,
        'ctx_data': ctx_data,
        'all_cov': all_cov,
        'all_cats': all_cats,
    }

# ══════════════════════════════════════════════════════════════════════════════
# FIGURE — Coverage distributions + Context stacked bars
# ══════════════════════════════════════════════════════════════════════════════
print('\n\nPlotting figure...')

fig = plt.figure(figsize=(14, 10))
gs_outer = gridspec.GridSpec(2, 1, figure=fig, hspace=0.42,
                              top=0.93, bottom=0.07, left=0.07, right=0.97)

# Row 0: Coverage
# Row 1: Context

gs_cov = gridspec.GridSpecFromSubplotSpec(1, 2, subplot_spec=gs_outer[0], wspace=0.30)
gs_ctx = gridspec.GridSpecFromSubplotSpec(1, 2, subplot_spec=gs_outer[1], wspace=0.30)

genome_labels = ['Female (PaproF)', 'Male (PaproM)']

# ── Coverage violin plots ─────────────────────────────────────────────────
for col, glabel in enumerate(genome_labels):
    ax = fig.add_subplot(gs_cov[col])
    res = all_results[glabel]
    parts_data = [np.log10(res['cov_data'][cat]) for cat in CATS if len(res['cov_data'].get(cat, [])) > 0]
    cats_present = [cat for cat in CATS if len(res['cov_data'].get(cat, [])) > 0]

    vp = ax.violinplot(parts_data, positions=range(len(cats_present)),
                       showmedians=True, showextrema=False, widths=0.6)
    for body, cat in zip(vp['bodies'], cats_present):
        body.set_facecolor(CAT_COL[cat])
        body.set_edgecolor('none')
        body.set_alpha(0.75)
    vp['cmedians'].set_color('#222222')
    vp['cmedians'].set_linewidth(1.5)

    # Add median text labels
    for i, cat in enumerate(cats_present):
        c = res['cov_data'][cat]
        med = np.median(c)
        ax.text(i, np.log10(med) + 0.08, f'{med:.0f}×',
                ha='center', va='bottom', fontsize=6.5, color='#333333')

    ax.set_xticks(range(len(cats_present)))
    ax.set_xticklabels([c.capitalize() for c in cats_present], fontsize=8)
    ax.set_ylabel('Coverage (log₁₀)', fontsize=8)
    ax.set_title(glabel, fontsize=9, pad=6)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    # y ticks in original scale
    yticks_log = [1, 1.5, 2, 2.5, 3]
    ax.set_yticks(yticks_log)
    ax.set_yticklabels([f'{10**y:.0f}' for y in yticks_log], fontsize=7)
    ax.axhline(np.log10(MIN_COV), color='#aaaaaa', lw=0.8, ls='--')
    ax.text(len(cats_present) - 0.5, np.log10(MIN_COV) + 0.04, f'min cov={MIN_COV}',
            ha='right', va='bottom', fontsize=6, color='#aaaaaa')

fig.text(0.5, 0.975, 'A  Coverage by methylation category', ha='center', va='top',
         fontsize=10, fontweight='bold')

# ── Context stacked bar charts ────────────────────────────────────────────
# Stack order (bottom to top):
stack_order = (
    ['intergenic', 'gene_only'] +
    ['gene+' + ('para' if m == 'pararetrovirus' else m) for m in MAJOR_ORDER] +
    list(MAJOR_ORDER)
)
stack_labels_display = (
    ['Intergenic', 'Gene body only'] +
    [f'Gene + {m}' for m in ['LTR', 'DNA', 'LINE', 'MITE', 'Pararetrovirus', 'Unknown']] +
    [f'{m} (repeat only)' for m in ['LTR', 'DNA', 'LINE', 'MITE', 'Pararetrovirus', 'Unknown']]
)

for col, glabel in enumerate(genome_labels):
    ax = fig.add_subplot(gs_ctx[col])
    res = all_results[glabel]
    x = np.arange(len(CATS))
    bottoms = np.zeros(len(CATS))

    handles = []
    for key, disp in zip(stack_order, stack_labels_display):
        heights = np.array([res['ctx_data'][cat].get(key, 0.0) * 100 for cat in CATS])
        color   = CTX_COLORS.get(key, '#cccccc')
        bars = ax.bar(x, heights, bottom=bottoms, color=color, width=0.55,
                      linewidth=0, label=disp)
        bottoms += heights
        if any(h > 0.5 for h in heights):
            handles.append(mpatches.Patch(color=color, label=disp, linewidth=0))

    ax.set_xticks(x)
    ax.set_xticklabels([c.capitalize() for c in CATS], fontsize=8)
    ax.set_ylabel('% of sites', fontsize=8)
    ax.set_ylim(0, 100)
    ax.set_title(glabel, fontsize=9, pad=6)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    if col == 1:
        ax.legend(handles=handles[::-1], fontsize=6, loc='upper left',
                  bbox_to_anchor=(1.02, 1.0), frameon=False,
                  handlelength=0.9, handletextpad=0.4, labelspacing=0.35)

fig.text(0.5, 0.495, 'B  Genomic context by methylation category', ha='center', va='top',
         fontsize=10, fontweight='bold')

outpath = f'{OUTDIR}/PaproF_PaproM_methylation_context'
fig.savefig(f'{outpath}.pdf', dpi=300, bbox_inches='tight', facecolor='white')
fig.savefig(f'{outpath}.png', dpi=150, bbox_inches='tight', facecolor='white')
print(f'  Saved {outpath}.pdf/png')
plt.close(fig)

print('\nDone.')
