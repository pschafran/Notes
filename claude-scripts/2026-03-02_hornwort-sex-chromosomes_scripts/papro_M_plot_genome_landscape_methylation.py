"""
Genome landscape: protein-coding genes + repeats + CpG methylation for PaproM
Nature/Science style — muted pastel palette, Liberation Sans font

Figure 1 (tracks): per-scaffold stacked area (genes | unannotated | repeats)
                   + 5mC methylation % line on twin y-axis
Figure 2 (summary): five-panel comparison using the same colour scheme
"""
import bisect
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.patches as mpatches
import matplotlib.ticker as ticker
from matplotlib.lines import Line2D
from collections import defaultdict

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

BASE       = '/media/data/projects/hornwort_sex_chromosomes/analysis/Paraphymatoceros_proskaueri/Papro252-3/final_genome_prep'
GFF        = f'{BASE}/PaproM.renamed_renamed.fasta.mod.EDTA.TEanno.gff3'
GTF        = f'{BASE}/braker_renamed.gtf'
FASTA      = f'{BASE}/PaproM.renamed_renamed.fasta'
METHYL_BED = f'{BASE}/../methylation/PaproM_genome.5mCG_5hmCG.bed'
OUT        = f'{BASE}/PaproM_genome_landscape'

SCAFFOLDS = ['PaproM.S1', 'PaproM.S2', 'PaproM.S3', 'PaproM.S4', 'PaproM.S5']
SLABELS   = ['S1', 'S2', 'S3', 'S4', 'S5']
SEX_CHR   = 'PaproM.S5'   # V chromosome
WINDOW, STEP = 200_000, 50_000
METHYL_MIN_COV = 5

GENE_COLOR    = '#5a87a5'
INTRON_COLOR  = '#96b8cc'
UNANN_COLOR   = '#dcdcdc'
METHYL_COLOR  = '#6a3d7a'   # muted violet for 5mC line
HMETH_COLOR   = '#b090c0'   # soft mauve for 5hmC dots

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

# S5 highlighted as sex chr (index 4)
SCAF_PALETTE = ['#6aab98', '#6aab98', '#6aab98', '#6aab98', '#c47a4a']

WIN_SCAF_COL = {
    'PaproM.S1': '#a8c8c0',
    'PaproM.S2': '#6aab98',
    'PaproM.S3': '#8ac0b0',
    'PaproM.S4': '#88a8b8',
    'PaproM.S5': '#c47a4a',   # V chromosome
}

def get_major(cls):
    major = cls.split('/')[0] if '/' in cls else cls
    if major == 'TIR':
        return 'DNA'
    return major

def load_methylation(bed_file, scaffolds, min_cov=5):
    """Load 5mC and 5hmC site data for windowed methylation and per-scaffold means.

    Returns
    -------
    data  : {seqid: (pos_m, met_m, pos_h, met_h)}  — position-sorted numpy arrays
    means : {seqid: (mC_mean_pct, hmC_mean_pct)}
    """
    raw_m = {s: ([], []) for s in scaffolds}
    raw_h = {s: ([], []) for s in scaffolds}
    scaf_set = set(scaffolds)
    with open(bed_file) as fh:
        for line in fh:
            p = line.split('\t')
            if len(p) < 11:
                continue
            seqid = p[0]
            if seqid not in scaf_set:
                continue
            mod = p[3]
            if mod not in ('m', 'h'):
                continue
            if int(p[9]) < min_cov:
                continue
            pos = int(p[1]); met = float(p[10])
            if mod == 'm':
                raw_m[seqid][0].append(pos)
                raw_m[seqid][1].append(met)
            else:
                raw_h[seqid][0].append(pos)
                raw_h[seqid][1].append(met)
    data  = {}
    means = {}
    for seqid in scaffolds:
        # 5mC
        if raw_m[seqid][0]:
            pos_m = np.array(raw_m[seqid][0], dtype=np.int64)
            met_m = np.array(raw_m[seqid][1])
            idx   = np.argsort(pos_m)
            pos_m, met_m = pos_m[idx], met_m[idx]
        else:
            pos_m, met_m = np.array([], dtype=np.int64), np.array([])
        # 5hmC
        if raw_h[seqid][0]:
            pos_h = np.array(raw_h[seqid][0], dtype=np.int64)
            met_h = np.array(raw_h[seqid][1])
            idx   = np.argsort(pos_h)
            pos_h, met_h = pos_h[idx], met_h[idx]
        else:
            pos_h, met_h = np.array([], dtype=np.int64), np.array([])
        data[seqid]  = (pos_m, met_m, pos_h, met_h)
        means[seqid] = (
            float(np.mean(met_m)) if len(met_m) else float('nan'),
            float(np.mean(met_h)) if len(met_h) else float('nan'),
        )
    return data, means

def win_meth(pos_arr, met_arr, ws, we):
    """Mean methylation % of sites with position in [ws, we]."""
    if len(pos_arr) == 0:
        return float('nan')
    lo = bisect.bisect_left(pos_arr, ws)
    hi = bisect.bisect_right(pos_arr, we)
    if lo >= hi:
        return float('nan')
    return float(np.mean(met_arr[lo:hi]))

# ── 1. Genome sizes ────────────────────────────────────────────────────────
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

# ── 2. Repeats ─────────────────────────────────────────────────────────────
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

# ── 3. Genes ────────────────────────────────────────────────────────────────
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

# ── 4. Methylation ──────────────────────────────────────────────────────────
print('Loading methylation...')
meth_data, meth_means = load_methylation(METHYL_BED, SCAFFOLDS, METHYL_MIN_COV)
for s in SCAFFOLDS:
    mC, hmC = meth_means[s]
    print(f'  {s}: 5mC={mC:.1f}%  5hmC={hmC:.1f}%')

# ── 5. Interval helpers ─────────────────────────────────────────────────────
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
    if not merged:
        return 0.0
    total = 0
    for s, e in merged:
        if e < ws: continue
        if s > we: break
        total += min(e, we) - max(s, ws) + 1
    return total / (we - ws + 1)

# ── 6. Build merged interval lists ────────────────────────────────────────
print('Building merged intervals...')
merged_rep    = {}
merged_allrep = {}
merged_gene   = {}
merged_exon   = {}
merged_union  = {}

for seqid in SCAFFOLDS:
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

# ── 7. Sliding-window computation ──────────────────────────────────────────
print('Computing windows...')
win_data = {}
for seqid in SCAFFOLDS:
    size = genome[seqid]
    positions = np.arange(0, size - WINDOW + 1, STEP, dtype=int)
    if len(positions) == 0:
        positions = np.array([0])
    n = len(positions)
    gene_only_arr   = np.zeros(n)
    unannotated_arr = np.zeros(n)
    rep_arr  = {maj: np.zeros(n) for maj in MAJOR_ORDER}
    mC_arr   = np.full(n, float('nan'))
    hmC_arr  = np.full(n, float('nan'))
    pos_m, met_m, pos_h, met_h = meth_data[seqid]
    for i, ws in enumerate(positions):
        we = min(int(ws) + WINDOW - 1, size - 1)
        rep_total  = cov_in_win(merged_allrep[seqid], ws, we)
        union_frac = cov_in_win(merged_union[seqid], ws, we)
        gene_only_arr[i]   = max(0.0, union_frac - rep_total)
        unannotated_arr[i] = max(0.0, 1.0 - union_frac)
        raw = {maj: cov_in_win(merged_rep[seqid][maj], ws, we) for maj in MAJOR_ORDER}
        raw_sum = sum(raw.values())
        if raw_sum > 1e-9:
            scale = rep_total / raw_sum
            for maj in MAJOR_ORDER:
                rep_arr[maj][i] = raw[maj] * scale
        mC_arr[i]  = win_meth(pos_m, met_m, ws, we)
        hmC_arr[i] = win_meth(pos_h, met_h, ws, we)
    win_data[seqid] = {
        'pos_mb':      positions / 1e6,
        'size_mb':     size / 1e6,
        'gene_only':   gene_only_arr,
        'unannotated': unannotated_arr,
        'rep':         rep_arr,
        'mC':          mC_arr,
        'hmC':         hmC_arr,
    }

# ── 8. Summary statistics ───────────────────────────────────────────────────
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
    d['mC_mean']     = meth_means[seqid][0]
    d['hmC_mean']    = meth_means[seqid][1]
    summary[seqid] = d

# ═══════════════════════════════════════════════════════════════════════════
# FIGURE 1 — Landscape tracks
# ═══════════════════════════════════════════════════════════════════════════
print('Plotting Figure 1: landscape tracks...')

fig1 = plt.figure(figsize=(18, 14))
fig1.patch.set_facecolor('white')

gs1 = gridspec.GridSpec(
    len(SCAFFOLDS), 1, figure=fig1,
    top=0.84, bottom=0.07,
    left=0.07, right=0.93,
    hspace=0.28,
)
axs = [fig1.add_subplot(gs1[i]) for i in range(len(SCAFFOLDS))]

STACK_LABELS = ['Protein-coding genes', 'Unannotated'] + \
               ['LTR retrotransposons', 'DNA transposons', 'LINEs', 'MITEs',
                'Pararetrovirus', 'Unknown repeat']
STACK_COLORS = [GENE_COLOR, UNANN_COLOR] + [MAJOR_COLORS[m] for m in MAJOR_ORDER]

for i, seqid in enumerate(SCAFFOLDS):
    ax    = axs[i]
    wd    = win_data[seqid]
    pos   = wd['pos_mb']
    size  = wd['size_mb']
    is_sex = seqid == SEX_CHR

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

    ngenes  = summary[seqid]['n_genes']
    rep_pct = summary[seqid]['rep_total']
    mC_pct  = summary[seqid]['mC_mean']
    sex_note = '   [V chromosome]' if is_sex else ''
    lc = '#7a2020' if is_sex else '#1a1a1a'
    ax.text(0.0, 1.04,
            f'S{i+1}  ·  {size:.1f} Mb  ·  {ngenes:,} genes  ·  '
            f'{rep_pct:.0f}% repetitive  ·  {mC_pct:.1f}% 5mC{sex_note}',
            transform=ax.transAxes, fontsize=8,
            fontweight='bold', va='bottom', color=lc,
            clip_on=False)

    # ── twin axis: 5mC methylation line ─────────────────────────────────
    ax2 = ax.twinx()
    mC_win = wd['mC']
    valid  = ~np.isnan(mC_win)
    if valid.any():
        ax2.plot(pos[valid], mC_win[valid], color=METHYL_COLOR,
                 lw=1.2, zorder=5, solid_capstyle='round')
        ax2.fill_between(pos[valid], 0, mC_win[valid],
                         color=METHYL_COLOR, alpha=0.08, zorder=4)
    ax2.set_xlim(0, size)
    ax2.set_ylim(0, 100)
    ax2.set_yticks([0, 50, 100])
    ax2.tick_params(axis='y', labelsize=6, colors=METHYL_COLOR,
                    length=2.5, width=0.5, pad=1)
    ax2.set_ylabel('5mC\n(%)', fontsize=6, color=METHYL_COLOR, labelpad=2)
    ax2.spines['right'].set_color(METHYL_COLOR)
    ax2.spines['right'].set_linewidth(0.5)
    ax2.spines['top'].set_visible(False)

meth_line = Line2D([0], [0], color=METHYL_COLOR, lw=1.6,
                   label='5mC methylation (%)')
legend_patches = [mpatches.Patch(color=c, label=l, linewidth=0)
                  for c, l in zip(STACK_COLORS, STACK_LABELS)]
axs[0].legend(handles=legend_patches + [meth_line], fontsize=7.5, ncol=5,
              loc='lower center', frameon=False,
              handlelength=1.0, handletextpad=0.5, columnspacing=1.0,
              bbox_to_anchor=(0.5, 1.52))

fig1.text(0.5, 0.960,
          'Paraphymatoceros proskaueri male (PaproM)  —  Genome landscape',
          ha='center', va='top', fontsize=11, fontweight='bold', color='#1a1a1a')
fig1.text(0.5, 0.945,
          '200 kb sliding window, 50 kb step',
          ha='center', va='top', fontsize=8, color='#555555')

fig1.savefig(f'{OUT}_tracks.pdf', dpi=300, bbox_inches='tight', facecolor='white')
fig1.savefig(f'{OUT}_tracks.png', dpi=180, bbox_inches='tight', facecolor='white')
print('  Saved tracks.')
plt.close(fig1)

# ═══════════════════════════════════════════════════════════════════════════
# FIGURE 2 — Summary panels (3 × 2, with methylation panel spanning bottom row)
# ═══════════════════════════════════════════════════════════════════════════
print('Plotting Figure 2: summary panels...')

fig2 = plt.figure(figsize=(12, 11))
fig2.patch.set_facecolor('white')

gs2 = gridspec.GridSpec(3, 2, figure=fig2,
                         left=0.08, right=0.97, top=0.91, bottom=0.07,
                         hspace=0.68, wspace=0.36,
                         height_ratios=[1, 1, 0.8])

ax_A = fig2.add_subplot(gs2[0, 0])
ax_B = fig2.add_subplot(gs2[0, 1])
ax_C = fig2.add_subplot(gs2[1, 0])
ax_D = fig2.add_subplot(gs2[1, 1])
ax_E = fig2.add_subplot(gs2[2, :])   # spans both columns

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

# ── A: Scatter ─────────────────────────────────────────────────────────────
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

# ── B: Gene size violin ─────────────────────────────────────────────────────
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

# ── C: Gene density ─────────────────────────────────────────────────────────
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

# ── D: Genome compartment composition ──────────────────────────────────────
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
    return [s['exon_pct'], s['intron_pct'], s['unannotated_pct'],
            s['LTR'], s['DNA'],
            s['LINE'] + s['MITE'] + s.get('pararetrovirus', 0),
            s['Unknown']]

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
ax_D.legend(fontsize=7, ncol=2, loc='upper right', frameon=False,
            handlelength=0.9, handletextpad=0.4, columnspacing=0.8)

# ── E: CpG methylation per scaffold ────────────────────────────────────────
mC_vals  = [summary[s]['mC_mean']  for s in SCAFFOLDS]
hmC_vals = [summary[s]['hmC_mean'] for s in SCAFFOLDS]

bars_e = ax_E.bar(x, mC_vals, bw, color=SCAF_PALETTE, linewidth=0, zorder=3,
                  label='5mC methylation')
bar_labels(ax_E, bars_e, mC_vals, fmt='{:.1f}')

ax_E.set_xticks(x)
ax_E.set_xticklabels(SLABELS, fontsize=8)
ax_E.set_ylabel('Mean 5mC (%)', fontsize=8)
ax_E.set_ylim(0, 100)
ax_E.set_xlim(-0.55, len(SCAFFOLDS) - 0.45)
style_ax(ax_E)
ax_E.set_title(f'CpG methylation per scaffold  (5mC and 5hmC, Nvalid_cov ≥ {METHYL_MIN_COV})',
               fontsize=8.5, pad=5)
ax_E.text(-0.06, 1.08, 'E', transform=ax_E.transAxes,
          fontsize=11, fontweight='bold', va='top')

# 5hmC on secondary y-axis
ax_E2 = ax_E.twinx()
ax_E2.plot(x, hmC_vals, 'o-', color=HMETH_COLOR, lw=1.2, ms=5,
           zorder=5, label='5hmC methylation')
ax_E2.set_ylim(0, 10)
ax_E2.set_yticks([0, 5, 10])
ax_E2.tick_params(axis='y', labelsize=6.5, colors=HMETH_COLOR,
                  length=2.5, width=0.5)
ax_E2.set_ylabel('Mean 5hmC (%)', fontsize=7, color=HMETH_COLOR, labelpad=2)
ax_E2.spines['right'].set_color(HMETH_COLOR)
ax_E2.spines['right'].set_linewidth(0.5)
ax_E2.spines['top'].set_visible(False)

# Combined legend for panel E
h1 = mpatches.Patch(color=SCAF_PALETTE[0], label='5mC methylation', linewidth=0)
h2 = Line2D([0], [0], color=HMETH_COLOR, marker='o', lw=1.2, ms=5,
            label='5hmC methylation (right axis)')
ax_E.legend(handles=[h1, h2], fontsize=7, loc='upper left',
            frameon=False, handlelength=1.0, handletextpad=0.5)

fig2.text(0.5, 0.975,
          'Paraphymatoceros proskaueri male (PaproM)  —  Gene, repeat and methylation summary',
          ha='center', va='top', fontsize=10.5, fontweight='bold', color='#1a1a1a')

fig2.savefig(f'{OUT}_summary.pdf', dpi=300, bbox_inches='tight', facecolor='white')
fig2.savefig(f'{OUT}_summary.png', dpi=180, bbox_inches='tight', facecolor='white')
print('  Saved summary.')
plt.close(fig2)

print('Done.')
