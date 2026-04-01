#!/usr/bin/env python3
"""
Plot BUSCO completeness vs. minimum contig size filter for LucruF_cmLunCruc20 genome.
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
from matplotlib.lines import Line2D

# ── Data ────────────────────────────────────────────────────────────────────
# min_size = 0 corresponds to LucruF_cmLunCruc20_genome.fasta (no filter applied)
data = [
    # min_size_bp, n_scaffolds, total_bp,      N50_bp,   C,   S,   D,  F,  M
    (         0,       3611,  778154215,  16000000,  809, 757,  52,  2, 11),
    (     50000,       1211,  702731932,  20000000,  807, 757,  50,  3, 12),
    (    100000,        463,  652442653,  23000000,  806, 765,  41,  3, 13),
    (    150000,        294,  631933632,  23000000,  806, 774,  32,  3, 13),
    (    200000,        201,  615682577,  23000000,  806, 780,  26,  3, 13),
    (    250000,        132,  600308963,  23000000,  806, 784,  22,  3, 13),
    (    300000,        102,  592189431,  23000000,  806, 785,  21,  3, 13),
    (    350000,         88,  587596532,  23000000,  806, 786,  20,  3, 13),
]

TOTAL_BUSCOS = 822

min_sizes   = np.array([d[0] for d in data])
n_scaffolds = np.array([d[1] for d in data])
total_bp    = np.array([d[2] for d in data])
n50_bp      = np.array([d[3] for d in data])
C           = np.array([d[4] for d in data])
S           = np.array([d[5] for d in data])
D           = np.array([d[6] for d in data])
F           = np.array([d[7] for d in data])
M           = np.array([d[8] for d in data])

# percentages
pS = S / TOTAL_BUSCOS * 100
pD = D / TOTAL_BUSCOS * 100
pF = F / TOTAL_BUSCOS * 100
pM = M / TOTAL_BUSCOS * 100

# ── Colours (muted pastels, Nature-style) ───────────────────────────────────
col_S   = '#4C9BE8'   # complete single-copy  – blue
col_D   = '#A8D1F5'   # complete duplicated   – light blue
col_F   = '#F5C842'   # fragmented            – amber
col_M   = '#E87070'   # missing               – coral/red
col_N50 = '#7A7A7A'   # N50 line              – grey

# ── x-axis tick labels ───────────────────────────────────────────────────────
x_labels = ['0\n(unfiltered)', '50 kb', '100 kb', '150 kb',
            '200 kb', '250 kb', '300 kb', '350 kb']
x = np.arange(len(min_sizes))
bar_w = 0.65

# ── Figure layout ───────────────────────────────────────────────────────────
fig = plt.figure(figsize=(9, 7))
gs  = fig.add_gridspec(2, 1, height_ratios=[3, 1.4], hspace=0.08)

ax1 = fig.add_subplot(gs[0])   # BUSCO stacked bars
ax2 = fig.add_subplot(gs[1], sharex=ax1)  # scaffold count + N50

# ── Panel 1 – stacked bar chart ──────────────────────────────────────────────
bars_S = ax1.bar(x, pS,       width=bar_w, color=col_S,   label='Complete (single-copy)', zorder=3)
bars_D = ax1.bar(x, pD,       width=bar_w, color=col_D,   label='Complete (duplicated)',  bottom=pS, zorder=3)
bars_F = ax1.bar(x, pF,       width=bar_w, color=col_F,   label='Fragmented',             bottom=pS+pD, zorder=3)
bars_M = ax1.bar(x, pM,       width=bar_w, color=col_M,   label='Missing',                bottom=pS+pD+pF, zorder=3)

# Annotate total complete % above each bar
for i, (s, d) in enumerate(zip(pS, pD)):
    total_c = s + d
    ax1.text(i, total_c + pF[i] + pM[i] + 0.8, f'{total_c:.1f}%',
             ha='center', va='bottom', fontsize=7.5, color='#333333')

ax1.set_ylabel('BUSCOs (%)', fontsize=10)
ax1.set_ylim(0, 108)
ax1.yaxis.set_major_locator(mticker.MultipleLocator(20))
ax1.set_xlim(-0.55, len(x) - 0.45)
ax1.spines[['top', 'right']].set_visible(False)
ax1.grid(axis='y', color='#e0e0e0', linewidth=0.7, zorder=0)
ax1.tick_params(axis='x', which='both', bottom=False, labelbottom=False)

# Legend inside panel 1
legend_handles = [
    Line2D([0], [0], marker='s', color='w', markerfacecolor=col_S, markersize=9, label='Complete (single-copy)'),
    Line2D([0], [0], marker='s', color='w', markerfacecolor=col_D, markersize=9, label='Complete (duplicated)'),
    Line2D([0], [0], marker='s', color='w', markerfacecolor=col_F, markersize=9, label='Fragmented'),
    Line2D([0], [0], marker='s', color='w', markerfacecolor=col_M, markersize=9, label='Missing'),
]
ax1.legend(handles=legend_handles, loc='lower left', fontsize=8,
           frameon=True, framealpha=0.9, edgecolor='#cccccc', ncol=2)

ax1.set_title('BUSCO completeness vs. minimum contig size filter\n'
              r'$\it{Lunularia\ cruciata}$ female assembly (LucruF\_cmLunCruc20)',
              fontsize=10.5, pad=8)

# ── Panel 2 – scaffold count (bars) + N50 (line) ────────────────────────────
col_scaf = '#B0C4DE'   # steel-blue-ish for scaffold bars

ax2b = ax2.twinx()

ax2.bar(x, n_scaffolds / 1000, width=bar_w, color=col_scaf,
        label='No. scaffolds', zorder=3, alpha=0.85)
ax2b.plot(x, n50_bp / 1e6, color=col_N50, marker='o', markersize=5,
          linewidth=1.6, label='Scaffold N50', zorder=4)

ax2.set_ylabel('Scaffolds (×10³)', fontsize=9)
ax2b.set_ylabel('N50 (Mb)', fontsize=9, color=col_N50)
ax2b.tick_params(axis='y', colors=col_N50)
ax2b.spines['right'].set_color(col_N50)

ax2.set_xticks(x)
ax2.set_xticklabels(x_labels, fontsize=8.5)
ax2.set_xlabel('Minimum contig length', fontsize=10)
ax2.set_ylim(0, n_scaffolds.max() / 1000 * 1.25)
ax2b.set_ylim(0, n50_bp.max() / 1e6 * 1.6)
ax2.spines[['top']].set_visible(False)
ax2b.spines[['top']].set_visible(False)
ax2.grid(axis='y', color='#e0e0e0', linewidth=0.7, zorder=0)
ax2.tick_params(axis='x', bottom=False)

# Combined legend for panel 2
handles2 = [
    Line2D([0], [0], marker='s', color='w', markerfacecolor=col_scaf,
           markersize=9, label='No. scaffolds', alpha=0.85),
    Line2D([0], [0], marker='o', color=col_N50, markersize=5,
           linewidth=1.6, label='Scaffold N50'),
]
ax2.legend(handles=handles2, loc='upper right', fontsize=8,
           frameon=True, framealpha=0.9, edgecolor='#cccccc')

# ── Font ────────────────────────────────────────────────────────────────────
for ax in [ax1, ax2, ax2b]:
    for label in ax.get_xticklabels() + ax.get_yticklabels():
        label.set_fontfamily('Liberation Sans')

# ── Save ────────────────────────────────────────────────────────────────────
out_base = 'LucruF_cmLunCruc20_busco_size_filter'
fig.savefig(f'{out_base}.pdf', dpi=300, bbox_inches='tight')
fig.savefig(f'{out_base}.png', dpi=300, bbox_inches='tight')
print(f'Saved {out_base}.pdf and {out_base}.png')
