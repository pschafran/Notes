"""
Genome landscape: protein-coding genes + repeats for LedusM
Nature/Science style -- muted pastel palette, Liberation Sans font

Figure 1 (tracks): per-scaffold stacked area showing genes | unannotated | repeats
Figure 2 (summary): six-panel comparison using the same colour scheme
"""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.patches as mpatches
import matplotlib.ticker as ticker
from collections import defaultdict

# Global style
plt.rcParams.update({
    'font.family':        'Liberation Sans',
    'font.size':          8,
    'axes.titlesize':     9,
    'axes.labelsize':     8,
    'xtick.labelsize':    7,
    'ytick.labelsize':    7,
    'legend.fontsize':    7,
    'axes.spines.top':    False,
    'axes.spines.right':  False,
    'axes.linewidth':     0.6,
    'xtick.major.width':  0.6,
    'ytick.major.width':  0.6,
    'xtick.minor.width':  0.4,
    'ytick.minor.width':  0.4,
    'xtick.major.size':   3,
    'ytick.major.size':   3,
    'xtick.direction':    'out',
    'ytick.direction':    'out',
    'grid.color':         '#d4d4d4',
    'grid.linewidth':     0.4,
    'figure.facecolor':   'white',
    'axes.facecolor':     'white',
    'pdf.fonttype':       42,
    'ps.fonttype':        42,
})

# Paths
GFF   = 'LedusM_repeat_annotations.gff'
GTF   = 'LedusM_gene_annotations.gtf'
FASTA = 'LedusM_genome.fasta'
OUT   = 'LedusM_genome_landscape'

SCAFFOLDS = ['LedusM.S1', 'LedusM.S2', 'LedusM.S3', 'LedusM.S4', 'LedusM.S5']
SLABELS   = ['S1', 'S2', 'S3', 'S4', 'S5']
WINDOW, STEP = 200_000, 50_000

# Shared colour palette
# Gene / structural colours
GENE_COLOR   = '#5a87a5'   # gene bodies (tracks) / exons (summary)
INTRON_COLOR = '#96b8cc'   # introns in summary panel F only
UNANN_COLOR  = '#dcdcdc'   # unannotated intergenic space

# Repeat major classes
MAJOR_ORDER  = ['LTR', 'DNA', 'LINE', 'MITE', 'pararetrovirus', 'Unknown']
MAJOR_COLORS = {
    'LTR':           '#c47a7a',
    'DNA':           '#7aab8a',
    'LINE':          '#8a8ab8',
    'MITE':          '#c8a86e',
    'pararetrovirus':'#9aacb8',
    'Unknown':       '#7aaabf',
    'Unmasked':      '#e8e8e8',
}

# Detailed repeat classes (for Panel B stacked bar)
ORDERED_CLASSES = [
    'LTR/Gypsy', 'LTR/Copia', 'LTR/unknown',
    'DNA/DTM', 'DNA/DTC', 'DNA/Helitron', 'DNA/DTT', 'DNA/DTH', 'DNA/DTA',
    'TIR/Tc1_Mariner',
    'LINE/unknown', 'MITE/DTM', 'MITE/DTH', 'MITE/DTC', 'MITE/DTA',
    'pararetrovirus', 'Unknown',
]
CLASS_COLORS = {
    'LTR/Gypsy':       '#c47a7a',
    'LTR/Copia':       '#e8a87c',
    'LTR/unknown':     '#f0c9a8',
    'DNA/DTM':         '#7aab8a',
    'DNA/DTC':         '#a8c9b0',
    'DNA/Helitron':    '#7aafc0',
    'DNA/DTT':         '#a8c8d8',
    'DNA/DTH':         '#c9c47a',
    'DNA/DTA':         '#ddd8a0',
    'TIR/Tc1_Mariner': '#b8a888',
    'LINE/unknown':    '#8a8ab8',
    'MITE/DTM':        '#c4b8d8',
    'MITE/DTH':        '#d4a8a0',
    'MITE/DTC':        '#c8b8a8',
    'MITE/DTA':        '#e8c4a0',
    'pararetrovirus':  '#9aacb8',
    'Unknown':         '#7aaabf',
    'Unmasked':        '#e8e8e8',
}

# Scaffold highlight colours (S5 = sex chromosome, index 4 highlighted)
SCAF_PALETTE = ['#7396b8', '#7396b8', '#7396b8', '#7396b8', '#b87070']

def get_major(cls):
    major = cls.split('/')[0] if '/' in cls else cls
    if major == 'TIR':
        return 'DNA'
    return major

# 1. Genome sizes
print('Reading genome...')
genome = {}
curr = None
with open(FASTA) as f:
    for line in f:
        line = line.rstrip()
        if line.startswith('>'):
            curr = line.split()[0][1:]
            genome[curr] = 0
        elif curr:
            genome[curr] += len(line)
total_genome = sum(genome.values())

# 2. Repeats
print('Parsing repeats...')
rep_by_seq_major = defaultdict(lambda: defaultdict(list))
rep_by_seq_cls   = defaultdict(lambda: defaultdict(list))
all_by_class     = defaultdict(list)

with open(GFF) as f:
    for line in f:
        if line.startswith('#'):
            continue
        p = line.rstrip().split('\t')
        if len(p) < 9:
            continue
        seqid, _, _, start, end, _, _, _, attrs = p
        start, end = int(start), int(end)
        attr_dict = {}
        for a in attrs.split(';'):
            if '=' in a:
                k, v = a.split('=', 1)
                attr_dict[k] = v
        cls = attr_dict.get('Classification', 'Unknown')
        rep_by_seq_major[seqid][get_major(cls)].append((start, end))
        rep_by_seq_cls[seqid][cls].append((start, end))
        all_by_class[cls].append((start, end, seqid))

# 3. Genes
print('Parsing genes...')
gene_by_seq = defaultdict(list)
exon_by_seq = defaultdict(list)

with open(GTF) as f:
    for line in f:
        if line.startswith('#'):
            continue
        p = line.rstrip().split('\t')
        if len(p) < 9:
            continue
        seqid, _, ftype, start, end, _, strand, _, _ = p
        start, end = int(start), int(end)
        if ftype == 'gene':
            gene_by_seq[seqid].append((start, end, strand))
        elif ftype == 'exon':
            exon_by_seq[seqid].append((start, end))

# 4. Interval helpers
def merge_sorted(intervals):
    if not intervals:
        return 0, []
    iv = sorted(intervals)
    m = [list(iv[0])]
    for s, e in iv[1:]:
        if s <= m[-1][1]:
            m[-1][1] = max(m[-1][1], e)
        else:
            m.append([s, e])
    return sum(e - s + 1 for s, e in m), m

def merge_bp(intervals):
    return merge_sorted(intervals)[0]

def cov_in_win(merged, ws, we):
    """Fraction of window [ws, we] covered by pre-sorted merged intervals."""
    if not merged:
        return 0.0
    total = 0
    for s, e in merged:
        if e < ws:
            continue
        if s > we:
            break
        total += min(e, we) - max(s, ws) + 1
    return total / (we - ws + 1)

# 5. Build merged interval lists per scaffold
print('Building merged intervals...')
merged_rep    = {}   # seqid -> {major: [[s,e],...]}
merged_allrep = {}   # seqid -> merged list of ALL repeats (union)
merged_gene   = {}   # seqid -> merged list of gene bodies
merged_exon   = {}   # seqid -> merged list of exons
merged_union  = {}   # seqid -> merged list of gene union all repeats

for seqid in SCAFFOLDS:
    # Per-major-class repeats
    merged_rep[seqid] = {}
    all_rep_ivs = []
    for maj in MAJOR_ORDER:
        ivs = rep_by_seq_major[seqid][maj]
        merged_rep[seqid][maj] = merge_sorted(ivs)[1] if ivs else []
        all_rep_ivs.extend(ivs)

    merged_allrep[seqid] = merge_sorted(all_rep_ivs)[1]

    gene_ivs = [(s, e) for s, e, _ in gene_by_seq[seqid]]
    merged_gene[seqid] = merge_sorted(gene_ivs)[1] if gene_ivs else []

    exon_ivs = exon_by_seq[seqid]
    merged_exon[seqid] = merge_sorted(exon_ivs)[1] if exon_ivs else []

    union_ivs = all_rep_ivs + gene_ivs
    merged_union[seqid] = merge_sorted(union_ivs)[1] if union_ivs else []

# 6. Sliding-window computation
print('Computing windows...')

win_data = {}
for seqid in SCAFFOLDS:
    size = genome[seqid]
    positions = np.arange(0, size - WINDOW + 1, STEP, dtype=int)
    if len(positions) == 0:
        positions = np.array([0])
    n = len(positions)

    gene_only_arr   = np.zeros(n)   # gene body not covered by any repeat
    unannotated_arr = np.zeros(n)   # neither gene nor repeat
    rep_arr = {maj: np.zeros(n) for maj in MAJOR_ORDER}

    for i, ws in enumerate(positions):
        we = min(int(ws) + WINDOW - 1, size - 1)

        rep_total  = cov_in_win(merged_allrep[seqid], ws, we)
        union_frac = cov_in_win(merged_union[seqid], ws, we)

        gene_only_arr[i]   = max(0.0, union_frac - rep_total)
        unannotated_arr[i] = max(0.0, 1.0 - union_frac)

        # Per-class repeat fractions, scaled so they sum to rep_total exactly
        raw = {maj: cov_in_win(merged_rep[seqid][maj], ws, we)
               for maj in MAJOR_ORDER}
        raw_sum = sum(raw.values())
        if raw_sum > 1e-9:
            scale = rep_total / raw_sum
            for maj in MAJOR_ORDER:
                rep_arr[maj][i] = raw[maj] * scale
        # else all stay zero

    win_data[seqid] = {
        'pos_mb':      positions / 1e6,
        'size_mb':     size / 1e6,
        'gene_only':   gene_only_arr,
        'unannotated': unannotated_arr,
        'rep':         rep_arr,
    }

# 7. Summary statistics
print('Computing summary...')
summary = {}
for seqid in SCAFFOLDS:
    size = genome[seqid]
    d = {}
    for maj in MAJOR_ORDER:
        d[maj] = merge_bp(merged_rep[seqid][maj]) / size * 100
    d['rep_total']   = merge_bp(merged_allrep[seqid]) / size * 100
    d['union_pct']   = merge_bp(merged_union[seqid]) / size * 100
    d['unmasked']    = 100 - d['rep_total']
    d['gene_bp']     = merge_bp([(s, e) for s, e, _ in gene_by_seq[seqid]])
    d['exon_bp']     = merge_bp(exon_by_seq[seqid])
    d['gene_pct']    = d['gene_bp'] / size * 100
    d['exon_pct']    = d['exon_bp'] / size * 100
    d['intron_pct']  = d['gene_pct'] - d['exon_pct']
    d['gene_only_pct'] = max(0.0, d['union_pct'] - d['rep_total'])
    d['unannotated_pct'] = max(0.0, 100.0 - d['union_pct'])
    d['n_genes']     = len(gene_by_seq[seqid])
    d['gene_dens']   = d['n_genes'] / (size / 1e6)
    summary[seqid] = d

# =========================================================================
# FIGURE 1 -- Genome landscape tracks (unified stacked area)
# =========================================================================
print('Plotting Figure 1: landscape tracks...')

fig1 = plt.figure(figsize=(18, 14))
fig1.patch.set_facecolor('white')

gs1 = gridspec.GridSpec(
    len(SCAFFOLDS), 1, figure=fig1,
    top=0.84, bottom=0.07,
    left=0.07, right=0.97,
    hspace=0.28,
)
axs = [fig1.add_subplot(gs1[i]) for i in range(len(SCAFFOLDS))]

# Stack layer order (bottom to top): gene_only | unannotated | LTR | DNA | LINE | MITE | para | Unknown
STACK_LABELS = ['Protein-coding genes', 'Unannotated'] + \
               ['LTR retrotransposons', 'DNA transposons', 'LINEs', 'MITEs',
                'Pararetrovirus', 'Unknown repeat']
STACK_COLORS = [GENE_COLOR, UNANN_COLOR] + [MAJOR_COLORS[m] for m in MAJOR_ORDER]

MAJOR_NICE = {
    'LTR': 'LTR retrotransposons', 'DNA': 'DNA transposons',
    'LINE': 'LINEs', 'MITE': 'MITEs',
    'pararetrovirus': 'Pararetrovirus', 'Unknown': 'Unknown repeat',
}

for i, seqid in enumerate(SCAFFOLDS):
    ax    = axs[i]
    wd    = win_data[seqid]
    pos   = wd['pos_mb']
    size  = wd['size_mb']
    is_sex = seqid == 'LedusM.S5'

    layers = (
        [wd['gene_only'] * 100, wd['unannotated'] * 100] +
        [wd['rep'][m] * 100 for m in MAJOR_ORDER]
    )

    ax.stackplot(pos, layers, colors=STACK_COLORS, linewidth=0, zorder=2)

    ax.set_xlim(0, size)
    ax.set_ylim(0, 100)
    ax.set_yticks([0, 50, 100])
    ax.set_yticklabels(['0', '50', '100'], fontsize=7)
    ax.yaxis.set_minor_locator(ticker.MultipleLocator(25))
    ax.tick_params(axis='y', length=3, width=0.6)
    ax.set_ylabel('% of\nwindow', fontsize=7, labelpad=3)

    for sp in ['top', 'right']:
        ax.spines[sp].set_visible(False)
    ax.spines['left'].set_linewidth(0.5)
    ax.spines['bottom'].set_linewidth(0.5)
    ax.grid(axis='y', linewidth=0.3, alpha=0.4, zorder=0)
    ax.set_axisbelow(True)

    if i < len(SCAFFOLDS) - 1:
        ax.tick_params(axis='x', bottom=False, labelbottom=False)
        ax.spines['bottom'].set_visible(False)
    else:
        ax.set_xlabel('Position (Mb)', fontsize=8)
        ax.tick_params(axis='x', labelsize=7, length=3, width=0.6)

    # Scaffold title -- placed above the track in the inter-track gap
    ngenes   = summary[seqid]['n_genes']
    rep_pct  = summary[seqid]['rep_total']
    sex_note = '   [putative sex chromosome]' if is_sex else ''
    lc = '#7a2020' if is_sex else '#1a1a1a'
    ax.text(0.0, 1.04,
            f'S{i+1}  .  {size:.1f} Mb  .  {ngenes:,} genes  .  {rep_pct:.0f}% repetitive{sex_note}',
            transform=ax.transAxes, fontsize=8,
            fontweight='bold', va='bottom', color=lc,
            clip_on=False)

# Legend -- centered above the tracks
legend_patches = [mpatches.Patch(color=c, label=l, linewidth=0)
                  for c, l in zip(STACK_COLORS, STACK_LABELS)]
axs[0].legend(handles=legend_patches, fontsize=7.5, ncol=4,
              loc='lower center', frameon=False,
              handlelength=1.0, handletextpad=0.5, columnspacing=1.0,
              bbox_to_anchor=(0.5, 1.52))

fig1.text(0.5, 0.960,
          'Leiosporoceros dussii male (LedusM)  --  Genome landscape',
          ha='center', va='top', fontsize=11, fontweight='bold', color='#1a1a1a')
fig1.text(0.5, 0.945,
          '200 kb sliding window, 50 kb step',
          ha='center', va='top', fontsize=8, color='#555555')

fig1.savefig(f'{OUT}_tracks.pdf', dpi=300, bbox_inches='tight', facecolor='white')
fig1.savefig(f'{OUT}_tracks.png', dpi=180, bbox_inches='tight', facecolor='white')
print('  Saved tracks.')
plt.close(fig1)

# =========================================================================
# FIGURE 2 -- Summary panels (2 x 2)
# =========================================================================
print('Plotting Figure 2: summary panels...')

fig2 = plt.figure(figsize=(12, 9))
fig2.patch.set_facecolor('white')

gs2 = gridspec.GridSpec(2, 2, figure=fig2,
                         left=0.08, right=0.97, top=0.91, bottom=0.14,
                         hspace=0.55, wspace=0.36)

ax_A = fig2.add_subplot(gs2[0, 0])   # gene vs repeat scatter
ax_B = fig2.add_subplot(gs2[0, 1])   # gene size violin
ax_C = fig2.add_subplot(gs2[1, 0])   # gene density bar
ax_D = fig2.add_subplot(gs2[1, 1])   # genome compartment composition

x  = np.arange(len(SCAFFOLDS))
bw = 0.52

def style_ax(ax, grid_axis='y'):
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_linewidth(0.6)
    ax.spines['bottom'].set_linewidth(0.6)
    if grid_axis:
        ax.grid(axis=grid_axis, linewidth=0.35, alpha=0.55, zorder=0)
    ax.set_axisbelow(True)

def bar_labels(ax, bars, vals, fmt='{:,}', dy_frac=0.02):
    ymax = ax.get_ylim()[1]
    for bar, val in zip(bars, vals):
        ax.text(bar.get_x() + bar.get_width() / 2,
                bar.get_height() + ymax * dy_frac,
                fmt.format(val),
                ha='center', va='bottom', fontsize=7, color='#333333')

# A: Scatter -- gene coverage vs repeat coverage
WIN_SCAF_COL = {
    'LedusM.S1': '#a8c4d8',
    'LedusM.S2': '#7396b8',
    'LedusM.S3': '#8aaa8a',
    'LedusM.S4': '#a8b878',
    'LedusM.S5': '#c47a7a',
}

for seqid in SCAFFOLDS:
    wd = win_data[seqid]
    rep_pct_win  = (1.0 - wd['gene_only'] - wd['unannotated']) * 100
    gene_pct_win = wd['gene_only'] * 100
    ax_A.scatter(rep_pct_win, gene_pct_win,
                 c=WIN_SCAF_COL[seqid], s=7, alpha=0.40,
                 linewidths=0, zorder=3)

ax_A.set_xlabel('Repeat coverage per window (%)', fontsize=8)
ax_A.set_ylabel('Gene coverage per window (%)', fontsize=8)
ax_A.set_xlim(0, 100)
ax_A.set_ylim(0, None)
style_ax(ax_A, grid_axis=None)
ax_A.grid(linewidth=0.35, alpha=0.45, zorder=0)
ax_A.set_title('Gene vs. repeat coverage\n(per 200 kb window)', fontsize=8.5, pad=5)
ax_A.text(-0.14, 1.11, 'A', transform=ax_A.transAxes,
          fontsize=11, fontweight='bold', va='top')

leg_a = [mpatches.Patch(color=WIN_SCAF_COL[s],
                         label=f'S{i+1}  ({summary[s]["n_genes"]:,} genes)',
                         linewidth=0)
         for i, s in enumerate(SCAFFOLDS)]
ax_A.legend(handles=leg_a, fontsize=7, loc='upper right',
            frameon=False, handlelength=0.9, handletextpad=0.4)

# B: Gene size violin
gene_sizes = {s: [] for s in SCAFFOLDS}
with open(GTF) as f:
    for line in f:
        if line.startswith('#'):
            continue
        p = line.rstrip().split('\t')
        if len(p) < 9 or p[2] != 'gene':
            continue
        if p[0] in gene_sizes:
            gene_sizes[p[0]].append(int(p[4]) - int(p[3]) + 1)

vp = ax_B.violinplot(
    [gene_sizes[s] for s in SCAFFOLDS],
    positions=x, widths=0.55,
    showmedians=True, showextrema=False,
)
for i, pc in enumerate(vp['bodies']):
    pc.set_facecolor(SCAF_PALETTE[i])
    pc.set_alpha(0.65)
    pc.set_linewidth(0)
vp['cmedians'].set_color('#333333')
vp['cmedians'].set_linewidth(1.5)

ax_B.set_xticks(x)
ax_B.set_xticklabels(SLABELS, fontsize=8)
ax_B.set_ylabel('Gene span (bp)', fontsize=8)
ax_B.set_yscale('log')
ax_B.set_ylim(80, 35000)
ax_B.yaxis.set_major_formatter(ticker.FuncFormatter(lambda v, _: f'{int(v):,}'))
style_ax(ax_B)
ax_B.set_title('Gene size distribution', fontsize=8.5, pad=5)
ax_B.text(-0.14, 1.08, 'B', transform=ax_B.transAxes,
          fontsize=11, fontweight='bold', va='top')

# C: Gene density per Mb
gdens = [summary[s]['gene_dens'] for s in SCAFFOLDS]
bars_c = ax_C.bar(x, gdens, bw, color=SCAF_PALETTE, linewidth=0, zorder=3)
ax_C.set_xticks(x)
ax_C.set_xticklabels(SLABELS, fontsize=8)
ax_C.set_ylabel('Genes per Mb', fontsize=8)
style_ax(ax_C)
bar_labels(ax_C, bars_c, [round(v, 1) for v in gdens], fmt='{:.1f}')
ax_C.set_title('Gene density per scaffold', fontsize=8.5, pad=5)
ax_C.text(-0.14, 1.08, 'C', transform=ax_C.transAxes,
          fontsize=11, fontweight='bold', va='top')

# D: Genome compartment composition
# Uses same colour scheme as the tracks figure
comp_defs = [
    ('Exons',          GENE_COLOR),
    ('Introns',        INTRON_COLOR),
    ('Unannotated',    UNANN_COLOR),
    ('LTR TEs',        MAJOR_COLORS['LTR']),
    ('DNA TEs',        MAJOR_COLORS['DNA']),
    ('LINE + MITE',    MAJOR_COLORS['LINE']),
    ('Unknown repeat', MAJOR_COLORS['Unknown']),
]

def get_comp_d(seqid):
    s = summary[seqid]
    exon_p   = s['exon_pct']
    intron_p = s['intron_pct']
    unann_p  = s['unannotated_pct']
    ltr_p    = s['LTR']
    dna_p    = s['DNA']
    other_te = s['LINE'] + s['MITE'] + s.get('pararetrovirus', 0)
    unk_p    = s['Unknown']
    return [exon_p, intron_p, unann_p, ltr_p, dna_p, other_te, unk_p]

bottoms_d = np.zeros(len(SCAFFOLDS))
for ci, (label, color) in enumerate(comp_defs):
    vals = np.array([get_comp_d(s)[ci] for s in SCAFFOLDS])
    ax_D.bar(x, vals, bw, bottom=bottoms_d,
             color=color, label=label, linewidth=0, zorder=3)
    bottoms_d += vals

ax_D.set_xticks(x)
ax_D.set_xticklabels(SLABELS, fontsize=8)
ax_D.set_ylabel('Scaffold composition (%)', fontsize=8)
ax_D.set_ylim(0, 100)
ax_D.set_xlim(-0.55, len(SCAFFOLDS) - 0.45)
style_ax(ax_D)
ax_D.set_title('Genome compartment composition', fontsize=8.5, pad=5)
ax_D.text(-0.14, 1.08, 'D', transform=ax_D.transAxes,
          fontsize=11, fontweight='bold', va='top')
ax_D.legend(fontsize=7, ncol=4,
            loc='upper center', frameon=False,
            handlelength=0.9, handletextpad=0.4, columnspacing=0.8,
            bbox_to_anchor=(0.5, -0.14))

fig2.text(0.5, 0.975,
          'Leiosporoceros dussii male (LedusM)  --  Gene and repeat summary',
          ha='center', va='top', fontsize=10.5, fontweight='bold', color='#1a1a1a')

fig2.savefig(f'{OUT}_summary.pdf', dpi=300, bbox_inches='tight', facecolor='white')
fig2.savefig(f'{OUT}_summary.png', dpi=180, bbox_inches='tight', facecolor='white')
print('  Saved summary.')
plt.close(fig2)

print('Done.')
