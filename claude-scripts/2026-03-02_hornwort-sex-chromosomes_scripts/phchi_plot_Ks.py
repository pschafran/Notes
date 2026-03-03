#!/usr/bin/env python3
"""
Ks distribution figure — PhchiF vs PhchiM reciprocal-best-hit pairs
4 panels:
  A: Ks histogram (ok pairs) with median line
  B: Ka vs Ks scatter, sex-chromosome vs autosome pairs highlighted
  C: Ks per-scaffold violin plots (scaffold_F), sex chr highlighted
  D: Ka/Ks histogram with Ka/Ks = 1 reference line

Nature/Science style — Liberation Sans, muted pastel palette.
"""
import math
import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.patches as mpatches
import matplotlib.ticker as ticker
import statistics

# ── Paths ──────────────────────────────────────────────────────────────────────
WORKDIR = '/media/data/projects/hornwort_sex_chromosomes/analysis/Phaeomegaceros_fimbriatus/sex_chromosome_analyses'
TSV_IN  = os.path.join(WORKDIR, 'rbh_ks', 'PhchiF_PhchiM_RBH_Ks.tsv')
OUT     = os.path.join(WORKDIR, 'rbh_ks', 'PhchiF_PhchiM_Ks_distribution')

# Sex chromosomes
F_SEX_SCAF = 'PhchiF.S4'   # U chromosome (female)
M_SEX_SCAF = 'PhchiM.S6'   # V chromosome (male)

# ── Style ──────────────────────────────────────────────────────────────────────
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
    'grid.color':         '#d0d0d0',
    'grid.linewidth':     0.4,
    'figure.facecolor':   'white',
    'axes.facecolor':     'white',
    'pdf.fonttype':       42,
    'ps.fonttype':        42,
})

# Colours
C_AUTO    = '#7396b8'   # autosome pairs — muted steel blue
C_SEX_F   = '#b87070'   # U chromosome (female S4) — muted terracotta
C_SEX_M   = '#c47a4a'   # V chromosome (male S6) — muted burnt orange
C_SEX_BOTH= '#9a5080'   # pairs where BOTH genes are on sex chromosomes — muted purple
C_MEDIAN  = '#444444'   # median line
C_REF     = '#c04040'   # Ka/Ks = 1 reference line
C_HIST    = '#7396b8'   # main histogram bars

# ── Load data ──────────────────────────────────────────────────────────────────
print('Loading data...')
rows = []
with open(TSV_IN) as fh:
    header = fh.readline().rstrip('\n').split('\t')
    col = {name: i for i, name in enumerate(header)}
    for line in fh:
        parts = line.rstrip('\n').split('\t')
        flag = parts[col['flag']]
        if flag != 'ok':
            continue
        def safe_float(s):
            try:
                v = float(s)
                return v if math.isfinite(v) else float('nan')
            except ValueError:
                return float('nan')
        Ks     = safe_float(parts[col['Ks']])
        Ka     = safe_float(parts[col['Ka']])
        Ka_Ks  = safe_float(parts[col['Ka_Ks']])
        scaf_F = parts[col['scaffold_F']]
        scaf_M = parts[col['scaffold_M']]
        if math.isnan(Ks):
            continue
        rows.append({
            'Ks': Ks, 'Ka': Ka, 'Ka_Ks': Ka_Ks,
            'scaffold_F': scaf_F, 'scaffold_M': scaf_M,
        })

print(f'Loaded {len(rows):,} ok pairs with finite Ks')

# Classify each pair
def classify(scaf_F, scaf_M):
    on_f = scaf_F == F_SEX_SCAF
    on_m = scaf_M == M_SEX_SCAF
    if on_f and on_m: return 'sex_both'
    if on_f:          return 'sex_F'
    if on_m:          return 'sex_M'
    return 'auto'

for r in rows:
    r['class'] = classify(r['scaffold_F'], r['scaffold_M'])

# Subset arrays
ks_all    = np.array([r['Ks']    for r in rows])
ka_all    = np.array([r['Ka']    for r in rows if math.isfinite(r['Ka'])])
omega_all = np.array([r['Ka_Ks'] for r in rows if math.isfinite(r['Ka_Ks'])])

auto_rows  = [r for r in rows if r['class'] == 'auto']
sexF_rows  = [r for r in rows if r['class'] in ('sex_F', 'sex_both')]
sexM_rows  = [r for r in rows if r['class'] in ('sex_M', 'sex_both')]
sexB_rows  = [r for r in rows if r['class'] == 'sex_both']

ks_median = float(np.median(ks_all))
print(f'Ks median: {ks_median:.4f}  n={len(ks_all):,}')

# Per-scaffold Ks for violin panel (scaffold_F, S1–S7)
F_SCAFFOLDS = ['PhchiF.S1', 'PhchiF.S2', 'PhchiF.S3', 'PhchiF.S4',
               'PhchiF.S5', 'PhchiF.S6', 'PhchiF.S7']
F_LABELS    = ['S1', 'S2', 'S3', 'S4\n(U chr)', 'S5', 'S6', 'S7']
scaf_ks = {s: [] for s in F_SCAFFOLDS}
for r in rows:
    s = r['scaffold_F']
    if s in scaf_ks:
        scaf_ks[s].append(r['Ks'])

print('Per-scaffold counts:', {s: len(v) for s, v in scaf_ks.items()})

# ── Figure layout ──────────────────────────────────────────────────────────────
fig = plt.figure(figsize=(12, 10))
gs  = gridspec.GridSpec(2, 2, figure=fig,
                        left=0.08, right=0.97, top=0.92, bottom=0.09,
                        wspace=0.34, hspace=0.45)

ax_A = fig.add_subplot(gs[0, 0])
ax_B = fig.add_subplot(gs[0, 1])
ax_C = fig.add_subplot(gs[1, 0])
ax_D = fig.add_subplot(gs[1, 1])

def style_ax(ax, grid_axis='y'):
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_linewidth(0.6)
    ax.spines['bottom'].set_linewidth(0.6)
    ax.set_axisbelow(True)
    if grid_axis:
        ax.grid(axis=grid_axis, linewidth=0.35, alpha=0.55, zorder=0)

# ── A: Ks histogram ────────────────────────────────────────────────────────────
ks_plot = ks_all[ks_all <= 3.0]
bins_a  = np.linspace(0, min(3.0, ks_plot.max() + 0.05), 60)

ax_A.hist(ks_plot, bins=bins_a, color=C_HIST, edgecolor='white',
          linewidth=0.3, alpha=0.85, zorder=3)
ax_A.axvline(ks_median, color=C_MEDIAN, linewidth=1.2,
             linestyle='--', zorder=5,
             label=f'Median Ks = {ks_median:.3f}')
for r_list, color, label in [
    (sexF_rows, C_SEX_F,  f'U chr (F-S4, n={len(sexF_rows):,})'),
    (sexM_rows, C_SEX_M,  f'V chr (M-S6, n={len(sexM_rows):,})'),
]:
    ks_sub = np.array([r['Ks'] for r in r_list if r['Ks'] <= 3.0])
    if len(ks_sub):
        ax_A.hist(ks_sub, bins=bins_a, color=color, edgecolor='white',
                  linewidth=0.3, alpha=0.75, zorder=4, label=label)

ax_A.set_xlabel('Ks (synonymous substitutions per synonymous site)', fontsize=8)
ax_A.set_ylabel('Number of gene pairs', fontsize=8)
ax_A.set_title(f'Ks distribution  (n = {len(ks_plot):,} pairs)', fontsize=8.5, pad=5)
style_ax(ax_A)
ax_A.legend(fontsize=6.5, frameon=False, loc='upper right')
ax_A.text(-0.13, 1.08, 'A', transform=ax_A.transAxes,
          fontsize=11, fontweight='bold', va='top')

# ── B: Ka vs Ks scatter ────────────────────────────────────────────────────────
scatter_rows = [r for r in rows
                if math.isfinite(r['Ka']) and math.isfinite(r['Ks'])
                and r['Ks'] <= 3.0 and r['Ka'] <= 3.0]

layers = [
    ('auto',     C_AUTO,   2,  0.20, f"Autosomal pairs  (n = {len(auto_rows):,})"),
    ('sex_F',    C_SEX_F,  7,  0.55, f"U chr (F-S4)  (n = {len([r for r in rows if r['class']=='sex_F']):,})"),
    ('sex_M',    C_SEX_M,  7,  0.55, f"V chr (M-S6)  (n = {len([r for r in rows if r['class']=='sex_M']):,})"),
    ('sex_both', C_SEX_BOTH, 10, 0.75, f"Both sex chr  (n = {len(sexB_rows):,})"),
]

for cls, color, sz, alpha, lbl in layers:
    sub = [r for r in scatter_rows if r['class'] == cls]
    if sub:
        ax_B.scatter([r['Ks'] for r in sub], [r['Ka'] for r in sub],
                     c=color, s=sz, alpha=alpha, linewidths=0, zorder=3,
                     label=lbl)

lim = min(3.0, max(r['Ks'] for r in scatter_rows) * 1.05) if scatter_rows else 1.0
ax_B.plot([0, lim], [0, lim], color=C_REF, linewidth=0.8,
          linestyle=':', zorder=2, label='Ka = Ks (neutral)')

ax_B.set_xlabel('Ks', fontsize=8)
ax_B.set_ylabel('Ka', fontsize=8)
ax_B.set_xlim(0, lim)
ax_B.set_ylim(0, lim)
ax_B.set_title('Ka vs Ks per RBH pair', fontsize=8.5, pad=5)
style_ax(ax_B, grid_axis=None)
ax_B.grid(linewidth=0.35, alpha=0.45, zorder=0)
ax_B.legend(fontsize=6.5, frameon=False, loc='upper left',
            handlelength=0.8, markerscale=2.0)
ax_B.text(-0.13, 1.08, 'B', transform=ax_B.transAxes,
          fontsize=11, fontweight='bold', va='top')

# ── C: Ks violin plots per scaffold_F ─────────────────────────────────────────
positions = list(range(1, len(F_SCAFFOLDS) + 1))
violin_data = [scaf_ks[s] for s in F_SCAFFOLDS]
scaf_colors = [C_SEX_F if s == F_SEX_SCAF else C_AUTO for s in F_SCAFFOLDS]

for pos, data, color in zip(positions, violin_data, scaf_colors):
    if len(data) >= 5:
        vp = ax_C.violinplot(data, positions=[pos], widths=0.7,
                             showmedians=True, showextrema=False)
        for pc in vp['bodies']:
            pc.set_facecolor(color)
            pc.set_alpha(0.65)
            pc.set_edgecolor('none')
        vp['cmedians'].set_color(C_MEDIAN)
        vp['cmedians'].set_linewidth(1.2)
    elif data:
        jitter = np.random.uniform(-0.15, 0.15, size=len(data))
        ax_C.scatter(np.array([pos] * len(data)) + jitter, data,
                     c=color, s=6, alpha=0.6, linewidths=0, zorder=3)

for pos, s in zip(positions, F_SCAFFOLDS):
    n = len(scaf_ks[s])
    ax_C.text(pos, -0.08, f'n={n:,}', ha='center', va='top',
              fontsize=6, color='#555555',
              transform=ax_C.get_xaxis_transform())

ax_C.set_xticks(positions)
ax_C.set_xticklabels(F_LABELS, fontsize=7.5)
ax_C.set_ylabel('Ks', fontsize=8)
ax_C.set_title('Ks distribution per female scaffold', fontsize=8.5, pad=5)
ax_C.set_xlim(0.3, len(F_SCAFFOLDS) + 0.7)
style_ax(ax_C)
ax_C.text(-0.13, 1.08, 'C', transform=ax_C.transAxes,
          fontsize=11, fontweight='bold', va='top')

auto_patch = mpatches.Patch(color=C_AUTO,  label='Autosome', linewidth=0)
sexF_patch = mpatches.Patch(color=C_SEX_F, label='U chromosome (S4)', linewidth=0)
ax_C.legend(handles=[auto_patch, sexF_patch], fontsize=6.5,
            frameon=False, loc='upper right')

# ── D: Ka/Ks histogram ────────────────────────────────────────────────────────
omega_plot = omega_all[(omega_all >= 0) & (omega_all <= 5)]
bins_d     = np.linspace(0, min(5.0, omega_plot.max() + 0.05), 60) if len(omega_plot) else [0, 1]
omega_median = float(np.median(omega_plot)) if len(omega_plot) else float('nan')

ax_D.hist(omega_plot, bins=bins_d, color=C_HIST, edgecolor='white',
          linewidth=0.3, alpha=0.85, zorder=3)
ax_D.axvline(1.0, color=C_REF, linewidth=1.0,
             linestyle='--', zorder=5, label='Ka/Ks = 1 (neutral)')
if math.isfinite(omega_median):
    ax_D.axvline(omega_median, color=C_MEDIAN, linewidth=1.2,
                 linestyle='--', zorder=5,
                 label=f'Median Ka/Ks = {omega_median:.3f}')

ax_D.set_xlabel('Ka/Ks (ω)', fontsize=8)
ax_D.set_ylabel('Number of gene pairs', fontsize=8)
ax_D.set_title(f'Ka/Ks distribution  (n = {len(omega_plot):,} pairs)', fontsize=8.5, pad=5)
style_ax(ax_D)
ax_D.legend(fontsize=6.5, frameon=False, loc='upper right')
ax_D.text(-0.13, 1.08, 'D', transform=ax_D.transAxes,
          fontsize=11, fontweight='bold', va='top')

# ── Figure title ───────────────────────────────────────────────────────────────
fig.text(0.5, 0.985,
         'Phaeomegaceros chiloensis — Synonymous substitution rates (Ks)',
         ha='center', va='top', fontsize=10.5, fontweight='bold', color='#1a1a1a')
fig.text(0.5, 0.970,
         f'Reciprocal-best-BLAST orthologs · YN00 method · n = {len(rows):,} ok pairs',
         ha='center', va='top', fontsize=8, color='#555555')

# ── Save ───────────────────────────────────────────────────────────────────────
fig.savefig(f'{OUT}.pdf', dpi=300, bbox_inches='tight', facecolor='white')
fig.savefig(f'{OUT}.png', dpi=180, bbox_inches='tight', facecolor='white')
print(f'Saved {OUT}.pdf / .png')

# ── Caption ────────────────────────────────────────────────────────────────────
caption = f"""Figure. Synonymous substitution rates (Ks) between reciprocal-best-BLAST ortholog pairs in Phaeomegaceros chiloensis female (PhchiF) and male (PhchiM) genomes.

(A) Histogram of Ks values for all {len(rows):,} ortholog pairs with valid YN00 estimates. The dashed vertical line indicates the median Ks ({ks_median:.3f}). Pairs involving the U chromosome (female S4, terracotta) or V chromosome (male S6, burnt orange) are overlaid as separate histograms.

(B) Scatter plot of Ka (dN) versus Ks (dS) for all pairs with finite Ka and Ks ≤ 3. Autosomal pairs (blue) and sex-chromosome pairs (coloured) are plotted separately. The dotted diagonal line marks Ka = Ks (neutral evolution, ω = 1).

(C) Violin plots of Ks distributions for gene pairs grouped by female scaffold (S1–S7). The U chromosome (S4, terracotta) contains {len(scaf_ks.get('PhchiF.S4', []))} pairs with valid Ks estimates. Horizontal bars within violins show medians.

(D) Histogram of Ka/Ks (ω) values for pairs with finite omega. The dashed red line marks Ka/Ks = 1 (strict neutrality); the black dashed line marks the median Ka/Ks ({f'{omega_median:.3f}' if math.isfinite(omega_median) else 'N/A'}).
"""

caption_path = f'{OUT}.caption.txt'
with open(caption_path, 'w') as fh:
    fh.write(caption)
print(f'Saved {caption_path}')
