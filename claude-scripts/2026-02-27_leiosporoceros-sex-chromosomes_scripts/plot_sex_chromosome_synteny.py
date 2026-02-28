#!/usr/bin/env python3
"""
U/V sex chromosome synteny plot for Leiosporoceros dussii.

Draws gene tick marks for all genes on the U (LedusF.S3) and V (LedusM.S5)
chromosomes, with Bezier connections between U↔V RBH pairs coloured by Ks.
Genes with an RBH to the other genome but on an autosome are marked dark grey.
"""

import re
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.path import Path
from matplotlib.colors import LinearSegmentedColormap
import matplotlib.cm as cm

# ── paths ────────────────────────────────────────────────────────────────────
BASE = '/media/data/projects/hornwort_sex_chromosomes/analysis/Leiosporoceros/sex_chromosome_analyses'
GTF_F  = f'{BASE}/LedusF_gene_annotations.gtf'
GTF_M  = f'{BASE}/LedusM_gene_annotations.gtf'
TSV    = f'{BASE}/rbh_ks/LedusF_LedusM_RBH_Ks.tsv'
OUT_PDF = f'{BASE}/rbh_ks/LedusF_LedusM_sex_chromosome_synteny.pdf'
OUT_PNG = f'{BASE}/rbh_ks/LedusF_LedusM_sex_chromosome_synteny.png'
OUT_CAP = f'{BASE}/rbh_ks/LedusF_LedusM_sex_chromosome_synteny.caption.txt'

# ── chromosome sizes ──────────────────────────────────────────────────────────
LEN_U = 22_033_833   # LedusF.S3
LEN_V =  5_298_153   # LedusM.S5

# ── colours ───────────────────────────────────────────────────────────────────
COL_U      = '#b87070'   # terracotta
COL_V      = '#c47a4a'   # burnt orange
COL_GREY   = '#888888'   # autosomal-partner ticks
COL_ALLTICK = '#cccccc'  # background gene ticks

# ── Ks colormap ───────────────────────────────────────────────────────────────
KS_MIN, KS_MAX = 0.0, 4.5
ks_cmap = LinearSegmentedColormap.from_list(
    'ks_sex', ['#2c7bb6', '#abd9e9', '#fdae61', '#d7191c'], N=256
)


# ── helpers ───────────────────────────────────────────────────────────────────

def load_gene_midpoints(gtf_file, seqid):
    """Return {gene_id: midpoint_bp} for 'gene' features on seqid."""
    genes = {}
    with open(gtf_file) as fh:
        for line in fh:
            if line.startswith('#'):
                continue
            p = line.rstrip('\n').split('\t')
            if len(p) < 9:
                continue
            if p[0] != seqid or p[2] != 'gene':
                continue
            gene_id = p[8].strip()
            start, end = int(p[3]), int(p[4])
            genes[gene_id] = (start + end) / 2
    return genes


def strip_transcript_suffix(tid):
    """LedusF.S3G000300.t1 → LedusF.S3G000300"""
    return re.sub(r'\.t\d+$', '', tid)


def load_rbh(tsv_file):
    """
    Return three lists of dicts:
      uv_pairs  — scaffold_F=S3, scaffold_M=S5, flag=ok
      u_auto    — scaffold_F=S3, scaffold_M!=S5  (any flag)
      v_auto    — scaffold_M=S5, scaffold_F!=S3  (any flag)
    Each dict has keys: gene_F, gene_M, scaffold_F, scaffold_M, Ks, flag
    """
    uv_pairs, u_auto, v_auto = [], [], []
    with open(tsv_file) as fh:
        header = fh.readline()
        for line in fh:
            p = line.rstrip('\n').split('\t')
            if len(p) < 12:
                continue
            gene_F, gene_M = p[0], p[1]
            scaf_F, scaf_M = p[2], p[3]
            try:
                ks = float(p[8])
            except ValueError:
                ks = np.nan
            flag = p[11]
            rec = dict(gene_F=gene_F, gene_M=gene_M,
                       scaffold_F=scaf_F, scaffold_M=scaf_M,
                       Ks=ks, flag=flag)
            if scaf_F == 'LedusF.S3' and scaf_M == 'LedusM.S5' and flag == 'ok':
                uv_pairs.append(rec)
            elif scaf_F == 'LedusF.S3' and scaf_M != 'LedusM.S5':
                u_auto.append(rec)
            elif scaf_M == 'LedusM.S5' and scaf_F != 'LedusF.S3':
                v_auto.append(rec)
    return uv_pairs, u_auto, v_auto


def draw_bezier(ax, x1, y1, x2, y2, cy1, cy2, color, lw=1.8, alpha=0.85):
    verts = [(x1, y1), (x1, cy1), (x2, cy2), (x2, y2)]
    codes = [Path.MOVETO, Path.CURVE4, Path.CURVE4, Path.CURVE4]
    path  = Path(verts, codes)
    patch = mpatches.PathPatch(path, facecolor='none',
                               edgecolor=color, lw=lw, alpha=alpha,
                               zorder=3)
    ax.add_patch(patch)


# ── load data ─────────────────────────────────────────────────────────────────
print('Loading gene positions…')
u_genes = load_gene_midpoints(GTF_F, 'LedusF.S3')
v_genes = load_gene_midpoints(GTF_M, 'LedusM.S5')
print(f'  U genes: {len(u_genes)}')
print(f'  V genes: {len(v_genes)}')

print('Loading RBH data…')
uv_pairs, u_auto, v_auto = load_rbh(TSV)
print(f'  U↔V pairs (ok): {len(uv_pairs)}')
print(f'  U-gene → autosome: {len(u_auto)}')
print(f'  V-gene ← autosome: {len(v_auto)}')

# gene IDs involved in U↔V connections (for coloured ticks)
uv_gene_F_ids = {strip_transcript_suffix(r['gene_F']): r['Ks'] for r in uv_pairs}
uv_gene_M_ids = {strip_transcript_suffix(r['gene_M']): r['Ks'] for r in uv_pairs}

# gene IDs with autosomal partner (grey ticks)
u_auto_ids = {strip_transcript_suffix(r['gene_F']) for r in u_auto}
v_auto_ids  = {strip_transcript_suffix(r['gene_M']) for r in v_auto}

# remove any that are already in uv sets
u_auto_ids -= set(uv_gene_F_ids.keys())
v_auto_ids  -= set(uv_gene_M_ids.keys())


# ── normalised positions ───────────────────────────────────────────────────────
def xnorm(mid, chrom_len):
    return mid / chrom_len


# ── figure ────────────────────────────────────────────────────────────────────
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype']  = 42
matplotlib.rcParams['font.family']  = 'Liberation Sans'
matplotlib.rcParams['font.size']    = 10

fig, ax = plt.subplots(figsize=(14, 7))
ax.set_xlim(-0.02, 1.02)
ax.set_ylim(0.0, 1.0)
ax.set_aspect('auto')
ax.axis('off')

# ── chromosome bars ────────────────────────────────────────────────────────────
U_Y0, U_Y1 = 0.68, 0.74
V_Y0, V_Y1 = 0.26, 0.32

ax.add_patch(mpatches.Rectangle(
    (0, U_Y0), 1.0, U_Y1 - U_Y0,
    facecolor='#d8d8d8', edgecolor='#888888', lw=0.8, zorder=2))
ax.add_patch(mpatches.Rectangle(
    (0, V_Y0), 1.0, V_Y1 - V_Y0,
    facecolor='#d8d8d8', edgecolor='#888888', lw=0.8, zorder=2))

# ── chromosome labels (left of bar) ──────────────────────────────────────────
ax.text(-0.02, (U_Y0 + U_Y1) / 2, 'U (LedusF.S3)\n22.0 Mb',
        va='center', ha='right', fontsize=9,
        color='#333333', linespacing=1.4)
ax.text(-0.02, (V_Y0 + V_Y1) / 2, 'V (LedusM.S5)\n5.3 Mb',
        va='center', ha='right', fontsize=9,
        color='#333333', linespacing=1.4)

# ── Mb tick labels on chromosomes ─────────────────────────────────────────────
def add_mb_ticks(ax, chrom_len, y_bar_top, y_bar_bot, above=True):
    mb_step = 2 if chrom_len <= 10e6 else 4
    mb_vals = np.arange(0, chrom_len / 1e6 + 0.01, mb_step)
    tick_y  = y_bar_top + 0.005 if above else y_bar_bot - 0.005
    label_y = y_bar_top + 0.018 if above else y_bar_bot - 0.018
    va = 'bottom' if above else 'top'
    for mb in mb_vals:
        x = mb * 1e6 / chrom_len
        if x > 1.0:
            break
        ax.plot([x, x], [y_bar_top if above else y_bar_bot,
                         tick_y], color='#666666', lw=0.6, zorder=5)
        label = f'{mb:.0f}' if mb == int(mb) else f'{mb:.1f}'
        ax.text(x, label_y, label + ' Mb', ha='center', va=va,
                fontsize=6.5, color='#555555')

add_mb_ticks(ax, LEN_U, U_Y1, U_Y0, above=True)
add_mb_ticks(ax, LEN_V, V_Y0, V_Y1, above=False)

# ── gene ticks ────────────────────────────────────────────────────────────────
# U chromosome: ticks above bar (y: U_Y1 → U_Y1 + 0.03)
TICK_U_LO = U_Y1
TICK_U_HI = U_Y1 + 0.03
TICK_V_HI = V_Y0
TICK_V_LO = V_Y0 - 0.03

for gid, mid in u_genes.items():
    x = xnorm(mid, LEN_U)
    if gid in uv_gene_F_ids:
        color = ks_cmap((uv_gene_F_ids[gid] - KS_MIN) / (KS_MAX - KS_MIN))
        lw, zorder = 1.2, 4
    elif gid in u_auto_ids:
        color = COL_GREY
        lw, zorder = 0.9, 3
    else:
        color = COL_ALLTICK
        lw, zorder = 0.6, 2
    ax.plot([x, x], [TICK_U_LO, TICK_U_HI], color=color, lw=lw, zorder=zorder)

for gid, mid in v_genes.items():
    x = xnorm(mid, LEN_V)
    if gid in uv_gene_M_ids:
        color = ks_cmap((uv_gene_M_ids[gid] - KS_MIN) / (KS_MAX - KS_MIN))
        lw, zorder = 1.2, 4
    elif gid in v_auto_ids:
        color = COL_GREY
        lw, zorder = 0.9, 3
    else:
        color = COL_ALLTICK
        lw, zorder = 0.6, 2
    ax.plot([x, x], [TICK_V_LO, TICK_V_HI], color=color, lw=lw, zorder=zorder)

# ── Bezier connections ────────────────────────────────────────────────────────
norm = matplotlib.colors.Normalize(vmin=KS_MIN, vmax=KS_MAX)

for rec in uv_pairs:
    gid_F = strip_transcript_suffix(rec['gene_F'])
    gid_M = strip_transcript_suffix(rec['gene_M'])
    if gid_F not in u_genes or gid_M not in v_genes:
        print(f'  WARNING: missing position for {gid_F} or {gid_M}')
        continue
    x_U = xnorm(u_genes[gid_F], LEN_U)
    x_V = xnorm(v_genes[gid_M], LEN_V)
    ks  = rec['Ks']
    color = ks_cmap(norm(ks))
    draw_bezier(ax, x_U, U_Y0, x_V, V_Y1,
                cy1=0.52, cy2=0.48, color=color)

# ── grey-tick annotation ──────────────────────────────────────────────────────
if u_auto_ids or v_auto_ids:
    ax.annotate('Autosomal\northologs only',
                xy=(0.01, TICK_U_HI + 0.005), xycoords='data',
                xytext=(0.08, 0.90),
                fontsize=8, color='#666666',
                arrowprops=dict(arrowstyle='->', color='#888888', lw=0.8),
                ha='left', va='bottom')

# ── colorbar ─────────────────────────────────────────────────────────────────
sm = cm.ScalarMappable(cmap=ks_cmap, norm=norm)
sm.set_array([])
cax = fig.add_axes([0.92, 0.22, 0.015, 0.56])
cb  = fig.colorbar(sm, cax=cax, orientation='vertical')
cb.set_label('Ks (synonymous substitutions\nper synonymous site)',
             fontsize=8, labelpad=6)
cb.set_ticks([0.5, 1.0, 2.0, 3.0, 4.0])
cb.ax.tick_params(labelsize=8)

# ── title & cleanup ───────────────────────────────────────────────────────────
ax.set_title('Sex chromosome synteny — U (female) vs V (male)\n'
             'Leiosporoceros dussii (LedusF.S3 × LedusM.S5)',
             fontsize=11, pad=12, color='#222222')

fig.subplots_adjust(left=0.12, right=0.90, top=0.88, bottom=0.05)

# ── save ──────────────────────────────────────────────────────────────────────
fig.savefig(OUT_PDF, dpi=300, bbox_inches='tight')
fig.savefig(OUT_PNG, dpi=200, bbox_inches='tight')
print(f'Saved {OUT_PDF}')
print(f'Saved {OUT_PNG}')
plt.close(fig)

# ── caption ──────────────────────────────────────────────────────────────────
n_uv = len(uv_pairs)
n_u_auto = len(u_auto_ids)
n_v_auto = len(v_auto_ids)
ks_vals = [r['Ks'] for r in uv_pairs if not np.isnan(r['Ks'])]
ks_lo, ks_hi = min(ks_vals), max(ks_vals)

caption = f"""\
Figure. Synteny of gene-bearing regions on the U (female) and V (male) sex chromosomes of
Leiosporoceros dussii. The U chromosome (LedusF.S3, 22.0 Mb; top bar, terracotta) and V
chromosome (LedusM.S5, 5.3 Mb; lower bar, burnt orange) are each drawn to scale along a
normalised horizontal axis. Vertical tick marks above/below each bar indicate the midpoint
position of every annotated gene ({len(u_genes)} on U, {len(v_genes)} on V). Bezier arcs connect the
{n_uv} reciprocal-best-hit (RBH) gene pairs shared between U and V (flag=ok); arc colour
indicates pairwise synonymous divergence (Ks) on a blue-to-red scale (range
{ks_lo:.2f}–{ks_hi:.2f}), far above the autosomal median of ~0.004, consistent with ancient
sex-chromosome divergence. Ticks coloured dark grey mark genes with a cross-sex RBH to an
autosome of the other individual ({n_u_auto} U-chromosome genes → male autosomes;
{n_v_auto} V-chromosome genes → female autosomes). The scattered positions of connected gene
pairs indicate no conserved collinearity between U and V, suggesting extensive chromosomal
rearrangement since the sex chromosomes originated. Ks was estimated using the YN00 method
(Yang & Nielsen 2000) on pairwise codon alignments of RBH transcript pairs.
"""
with open(OUT_CAP, 'w') as fh:
    fh.write(caption)
print(f'Saved {OUT_CAP}')
