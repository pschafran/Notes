"""
Repeat composition visualization for PhchiM (Phaeomegaceros chiloensis male, 14765-4)
Nature/Science style — muted pastel palette, Liberation Sans font
"""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.gridspec as gridspec
import matplotlib.ticker as ticker
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

BASE  = '/media/data/projects/hornwort_sex_chromosomes/analysis/Phaeomegaceros_fimbriatus/Phaeomegaceros_14765-4/final_genome_prep'
GFF   = '/media/data/projects/hornwort_sex_chromosomes/analysis/Phaeomegaceros_fimbriatus/Phaeomegaceros_14765-4/annotation/EDTA/Phchi4_bacteriaRemoved.fasta.mod.EDTA.anno/Phchi4_bacteriaRemoved.fasta.mod.EDTA.TEanno.gff3'
FASTA = f'{BASE}/PhchiM_genome.fasta'
OUT   = f'{BASE}/PhchiM_repeat_composition'

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
    'LINE/unknown':    '#8a8ab8',
    'Penelope':        '#b0a0cc',
    'MITE/DTM':        '#c4b8d8',
    'MITE/DTH':        '#d4a8a0',
    'MITE/DTC':        '#c8b8a8',
    'MITE/DTA':        '#c0c8a8',
    'MITE/DTT':        '#b8c0a0',
    'TIR/Tc1_Mariner': '#9abcb8',
    'TIR/Sola2':       '#80b0a8',
    'polinton':        '#88b088',
    'Unknown':         '#aaaaaa',
    'Unmasked':        '#e8e8e8',
}

MAJOR_COLORS = {
    'LTR':     '#c47a7a',
    'DNA':     '#7aab8a',
    'LINE':    '#8a8ab8',
    'MITE':    '#c8a86e',
    'Unknown': '#aaaaaa',
    'Unmasked':'#e8e8e8',
}

ORDERED_CLASSES = [
    'LTR/Gypsy', 'LTR/Copia', 'LTR/unknown',
    'LINE/unknown', 'Penelope',
    'DNA/DTM', 'DNA/DTC', 'DNA/DTH', 'DNA/DTA', 'DNA/DTT', 'DNA/Helitron',
    'MITE/DTM', 'MITE/DTC', 'MITE/DTH', 'MITE/DTA', 'MITE/DTT',
    'TIR/Tc1_Mariner', 'TIR/Sola2', 'polinton',
    'Unknown',
]

MAJOR_ORDER = ['LTR', 'DNA', 'LINE', 'MITE', 'Unknown']

SCAFFOLDS = ['PhchiM.S1', 'PhchiM.S2', 'PhchiM.S3', 'PhchiM.S4', 'PhchiM.S5', 'PhchiM.S6']
SLABELS   = ['S1\n(37.3 Mb)', 'S2\n(34.0 Mb)', 'S3\n(26.7 Mb)', 'S4\n(18.5 Mb)', 'S5\n(9.9 Mb)', 'S6\n(6.2 Mb)']
SEX_CHR   = 'PhchiM.S6'   # V chromosome

WINDOW, STEP = 200_000, 50_000


def get_major(cls):
    major = cls.split('/')[0] if '/' in cls else cls
    if major in ('TIR', 'polinton'):
        return 'DNA'
    if major == 'Penelope':
        return 'LINE'
    return major


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

print('Parsing GFF...')
by_seq_class = defaultdict(lambda: defaultdict(list))
all_by_class = defaultdict(list)

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
        by_seq_class[seqid][cls].append((start, end))
        all_by_class[cls].append((start, end, seqid))


def merge_intervals(intervals):
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


def window_coverage(merged, ws, we):
    total = 0
    for s, e in merged:
        if e < ws: continue
        if s > we: break
        total += min(e, we) - max(s, ws) + 1
    return total / (we - ws + 1)


print('Computing sliding windows...')
scaffold_windows = {}
for seqid in SCAFFOLDS:
    size = genome[seqid]
    positions = list(range(0, size - WINDOW + 1, STEP)) or [0]
    major_merged = {}
    for cls, ivs in by_seq_class[seqid].items():
        maj = get_major(cls)
        major_merged.setdefault(maj, []).extend(ivs)
    for maj in major_merged:
        _, major_merged[maj] = merge_intervals(major_merged[maj])
    win_data = {maj: np.zeros(len(positions)) for maj in MAJOR_ORDER}
    for i, ws in enumerate(positions):
        we = min(ws + WINDOW - 1, size - 1)
        for maj in MAJOR_ORDER:
            if major_merged.get(maj):
                win_data[maj][i] = window_coverage(major_merged[maj], ws, we)
    scaffold_windows[seqid] = {
        'positions': np.array(positions) / 1e6,
        'data': win_data,
        'size': size / 1e6,
    }

print('Computing coverage...')
class_coverage = {}
for cls in ORDERED_CLASSES:
    ivs_by_seq = defaultdict(list)
    for s, e, seqid in all_by_class.get(cls, []):
        ivs_by_seq[seqid].append((s, e))
    class_coverage[cls] = sum(merge_intervals(ivs)[0] for ivs in ivs_by_seq.values())

scaffold_class_cov = {}
for seqid in SCAFFOLDS:
    scaffold_class_cov[seqid] = {}
    for cls in ORDERED_CLASSES:
        bp, _ = merge_intervals(by_seq_class[seqid].get(cls, []))
        scaffold_class_cov[seqid][cls] = bp

major_cov = defaultdict(int)
for cls, bp in class_coverage.items():
    major_cov[get_major(cls)] += bp
major_cov['Unmasked'] = max(0, total_genome - sum(v for k, v in major_cov.items() if k != 'Unmasked'))

print('Plotting...')
fig = plt.figure(figsize=(16, 20))
fig.patch.set_facecolor('white')

n_scaf = len(SCAFFOLDS)
gs_top = gridspec.GridSpec(1, 2, figure=fig,
                           left=0.07, right=0.97, top=0.97, bottom=0.74,
                           wspace=0.32)
gs_mid = gridspec.GridSpec(1, 1, figure=fig,
                           left=0.07, right=0.97, top=0.70, bottom=0.54)
gs_bot = gridspec.GridSpec(n_scaf, 1, figure=fig,
                           left=0.07, right=0.97, top=0.49, bottom=0.04,
                           hspace=0.07)

ax_pie  = fig.add_subplot(gs_top[0, 0])
ax_bar  = fig.add_subplot(gs_top[0, 1])
ax_heat = fig.add_subplot(gs_mid[0, 0])
axs_win = [fig.add_subplot(gs_bot[i, 0]) for i in range(n_scaf)]

# A: Pie
pie_order  = ['LTR', 'Unknown', 'DNA', 'MITE', 'LINE', 'Unmasked']
pie_sizes  = [major_cov[m] for m in pie_order]
pie_colors = [MAJOR_COLORS[m] for m in pie_order]
pie_pcts   = [v / total_genome * 100 for v in pie_sizes]

wedges, _, autotexts = ax_pie.pie(
    pie_sizes, colors=pie_colors,
    autopct=lambda p: f'{p:.1f}%' if p > 2 else '',
    pctdistance=0.72, startangle=90,
    wedgeprops=dict(linewidth=0.5, edgecolor='white'),
    textprops=dict(fontsize=7),
)
for at in autotexts:
    at.set_fontsize(7); at.set_fontweight('bold'); at.set_color('#333333')

legend_patches = [
    mpatches.Patch(color=MAJOR_COLORS[m], label=f'{m}  {pie_pcts[i]:.1f}%', linewidth=0)
    for i, m in enumerate(pie_order)
]
ax_pie.legend(handles=legend_patches, loc='lower left',
              bbox_to_anchor=(-0.08, -0.08),
              frameon=False, fontsize=7.5, handlelength=1.2, handletextpad=0.5)
ax_pie.set_title('Genome-wide composition\n(non-overlapping)', fontsize=8.5, pad=6)
ax_pie.text(-0.10, 1.08, 'A', transform=ax_pie.transAxes,
            fontsize=11, fontweight='bold', va='top')

# B: Stacked bar per scaffold
x = np.arange(n_scaf)
bw = 0.52
bottoms = np.zeros(n_scaf)

for cls in ORDERED_CLASSES:
    vals = np.array([scaffold_class_cov[s][cls] / genome[s] * 100 for s in SCAFFOLDS])
    ax_bar.bar(x, vals, bw, bottom=bottoms, color=CLASS_COLORS[cls], linewidth=0, zorder=3)
    bottoms += vals

for i, s in enumerate(SCAFFOLDS):
    used = sum(scaffold_class_cov[s][c] for c in ORDERED_CLASSES) / genome[s] * 100
    ax_bar.bar(i, 100 - used, bw, bottom=used,
               color=CLASS_COLORS['Unmasked'], linewidth=0, zorder=3)

ax_bar.set_xticks(x)
ax_bar.set_xticklabels(SLABELS, fontsize=7.5)
ax_bar.set_ylabel('Scaffold composition (%)', fontsize=8)
ax_bar.set_ylim(0, 100)
ax_bar.set_xlim(-0.55, n_scaf - 0.45)
ax_bar.yaxis.set_minor_locator(ticker.MultipleLocator(10))
ax_bar.grid(axis='y', zorder=0, alpha=0.6)
ax_bar.set_axisbelow(True)
ax_bar.set_title('Repeat class composition per scaffold', fontsize=8.5, pad=6)
ax_bar.text(-0.12, 1.08, 'B', transform=ax_bar.transAxes,
            fontsize=11, fontweight='bold', va='top')

# C: Heatmap
from matplotlib.colors import LinearSegmentedColormap
N_BINS = 500
heat_data = np.full((n_scaf, N_BINS), np.nan)

for row_i, seqid in enumerate(SCAFFOLDS):
    size = genome[seqid]
    all_ivs = []
    for ivs in by_seq_class[seqid].values():
        all_ivs.extend(ivs)
    _, merged = merge_intervals(all_ivs)
    bin_size = size / N_BINS
    for b in range(N_BINS):
        bs, be = int(b * bin_size), int((b + 1) * bin_size) - 1
        heat_data[row_i, b] = window_coverage(merged, bs, be) if merged else 0

pastel_cmap = LinearSegmentedColormap.from_list(
    'pastel_yor', ['#f5f0e8', '#f0d0a8', '#d4907a', '#a85858'], N=256)

im = ax_heat.imshow(heat_data, aspect='auto', cmap=pastel_cmap,
                    vmin=0, vmax=1, interpolation='nearest')
ax_heat.set_yticks(range(n_scaf))
ax_heat.set_yticklabels([f'S{i+1}' for i in range(n_scaf)], fontsize=8)
ax_heat.set_xticks([0, N_BINS // 4, N_BINS // 2, 3 * N_BINS // 4, N_BINS - 1])
ax_heat.set_xticklabels(['0%', '25%', '50%', '75%', '100%'], fontsize=7)
ax_heat.set_xlabel('Relative position along scaffold', fontsize=8)
ax_heat.spines['top'].set_visible(False); ax_heat.spines['right'].set_visible(False)
ax_heat.spines['left'].set_linewidth(0.6); ax_heat.spines['bottom'].set_linewidth(0.6)
ax_heat.tick_params(length=3, width=0.6)

cbar = plt.colorbar(im, ax=ax_heat, orientation='vertical', pad=0.015, shrink=0.85, aspect=18)
cbar.set_label('Repeat density', fontsize=7.5)
cbar.ax.tick_params(labelsize=6.5, length=2, width=0.5)
cbar.outline.set_linewidth(0.5)
cbar.set_ticks([0, 0.25, 0.5, 0.75, 1.0])
cbar.set_ticklabels(['0%', '25%', '50%', '75%', '100%'])
ax_heat.set_title('Repeat density across scaffolds', fontsize=8.5, pad=6)
ax_heat.text(-0.04, 1.12, 'C', transform=ax_heat.transAxes,
             fontsize=11, fontweight='bold', va='top')

# D: Stacked-area per scaffold
major_nice = {
    'LTR':     'LTR retrotransposons',
    'DNA':     'DNA transposons',
    'LINE':    'LINEs + Penelope',
    'MITE':    'MITEs',
    'Unknown': 'Unknown',
}

for ax_i, (ax, seqid) in enumerate(zip(axs_win, SCAFFOLDS)):
    wd   = scaffold_windows[seqid]
    pos  = wd['positions']
    size = wd['size']

    stacks = [wd['data'][m] * 100 for m in MAJOR_ORDER]
    ax.stackplot(pos, stacks, colors=[MAJOR_COLORS[m] for m in MAJOR_ORDER], linewidth=0)

    ax.set_xlim(0, size)
    ax.set_ylim(0, 100)
    ax.set_yticks([0, 50, 100])
    ax.set_yticklabels(['0', '50', '100'], fontsize=6.5)
    ax.yaxis.set_minor_locator(ticker.MultipleLocator(25))
    ax.tick_params(axis='y', length=3, width=0.6)
    ax.set_ylabel('% covered', fontsize=7, labelpad=3)

    for sp in ['top', 'right']:
        ax.spines[sp].set_visible(False)
    for sp in ['left', 'bottom']:
        ax.spines[sp].set_linewidth(0.5)
    ax.grid(axis='y', linewidth=0.3, alpha=0.5, zorder=0)
    ax.set_axisbelow(True)

    if ax_i < n_scaf - 1:
        ax.tick_params(axis='x', bottom=False, labelbottom=False)
        ax.spines['bottom'].set_visible(False)
    else:
        ax.set_xlabel('Position (Mb)', fontsize=8)
        ax.tick_params(axis='x', labelsize=7, length=3, width=0.6)

    is_sex = seqid == SEX_CHR
    label_color = '#8a3030' if is_sex else '#333333'
    suffix = '  [V chromosome]' if is_sex else ''
    ax.text(0.005, 0.90,
            f'S{ax_i + 1}  ({size:.1f} Mb){suffix}',
            transform=ax.transAxes, fontsize=7.5,
            fontweight='bold', va='top', color=label_color)

legend_patches_d = [
    mpatches.Patch(color=MAJOR_COLORS[m], label=major_nice[m], linewidth=0)
    for m in MAJOR_ORDER
]
axs_win[0].legend(handles=legend_patches_d, fontsize=7, ncol=3,
                  loc='upper right', frameon=False,
                  handlelength=1.0, handletextpad=0.4,
                  bbox_to_anchor=(1.0, 1.55))
axs_win[0].text(-0.04, 1.6, 'D', transform=axs_win[0].transAxes,
                fontsize=11, fontweight='bold', va='top')
axs_win[0].set_title('Repeat density per scaffold  (200 kb sliding window, 50 kb step)',
                     fontsize=8.5, loc='left', pad=26)

fig.text(0.5, 0.995,
         'Phaeomegaceros chiloensis male (PhchiM, 14765-4)  —  Repetitive element composition',
         ha='center', va='top', fontsize=10, fontweight='bold', color='#222222')

fig.savefig(f'{OUT}.pdf', dpi=300, bbox_inches='tight', facecolor='white')
fig.savefig(f'{OUT}.png', dpi=180, bbox_inches='tight', facecolor='white')
print(f'Saved {OUT}.pdf / .png')
