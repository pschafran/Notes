"""
Male vs. female genome comparison: Phaeomegaceros chiloensis (PhchiF vs PhchiM)
Emphasis on sex chromosome differentiation (U chromosome in female S4, V chromosome in male S6)
Nature/Science style — muted pastel palette, Liberation Sans font
"""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.patches as mpatches
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

BASE       = '/media/data/projects/hornwort_sex_chromosomes/analysis/Phaeomegaceros_fimbriatus'
FEMALE_DIR = f'{BASE}/Phaeomegaceros_14765-5/final_genome_prep'
MALE_DIR   = f'{BASE}/Phaeomegaceros_14765-4/final_genome_prep'
OUT        = f'{BASE}/sex_chromosome_analyses/PhchiF_PhchiM_sex_chromosome_comparison'

F_SCAFFOLDS = ['PhchiF.S1', 'PhchiF.S2', 'PhchiF.S3', 'PhchiF.S4',
               'PhchiF.S5', 'PhchiF.S6', 'PhchiF.S7']
F_SEX       = 'PhchiF.S4'   # U chromosome (female)

M_SCAFFOLDS = ['PhchiM.S1', 'PhchiM.S2', 'PhchiM.S3', 'PhchiM.S4', 'PhchiM.S5', 'PhchiM.S6']
M_SEX       = 'PhchiM.S6'   # V chromosome (male)

F_AUTO  = '#7396b8'
F_SEX_C = '#b87070'   # U chromosome (terracotta)
M_AUTO  = '#6aab98'
M_SEX_C = '#c47a4a'   # V chromosome (burnt orange)

MAJOR_ORDER = ['LTR', 'DNA', 'LINE', 'MITE', 'Unknown']
MAJOR_COLORS = {
    'LTR':     '#c47a7a',
    'DNA':     '#7aab8a',
    'LINE':    '#8a8ab8',
    'MITE':    '#c8a86e',
    'Unknown': '#7aaabf',
    'Unmasked':'#e8e8e8',
}
GENE_COLOR   = '#5a87a5'
INTRON_COLOR = '#96b8cc'
UNANN_COLOR  = '#dcdcdc'

WINDOW, STEP = 200_000, 50_000

M_GFF = (f'{BASE}/Phaeomegaceros_14765-4/annotation/EDTA/'
         'Phchi4_bacteriaRemoved.fasta.mod.EDTA.anno/'
         'Phchi4_bacteriaRemoved.fasta.mod.EDTA.TEanno.gff3')
F_GFF = f'{FEMALE_DIR}/PhchiF_genome.fasta.mod.EDTA.TEanno.gff3'

def get_major(cls):
    major = cls.split('/')[0] if '/' in cls else cls
    if major in ('TIR', 'polinton'):
        return 'DNA'
    if major == 'Penelope':
        return 'LINE'
    return major

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

def load_genome(fasta_file):
    genome = {}; curr = None
    with open(fasta_file) as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith('>'):
                curr = line.split()[0][1:]
                genome[curr] = 0
            elif curr:
                genome[curr] += len(line)
    return genome

def load_repeats(gff_file):
    by_seq_major = defaultdict(lambda: defaultdict(list))
    all_rep      = defaultdict(list)
    with open(gff_file) as fh:
        for line in fh:
            if line.startswith('#'): continue
            p = line.rstrip().split('\t')
            if len(p) < 9: continue
            seqid = p[0]; start, end = int(p[3]), int(p[4]); attrs = p[8]
            attr_dict = {}
            for a in attrs.split(';'):
                if '=' in a:
                    k, v = a.split('=', 1)
                    attr_dict[k] = v
            cls = attr_dict.get('Classification', 'Unknown')
            maj = get_major(cls)
            by_seq_major[seqid][maj].append((start, end))
            all_rep[seqid].append((start, end))
    return by_seq_major, all_rep

def load_genes(gtf_file):
    gene_by_seq = defaultdict(list)
    exon_by_seq = defaultdict(list)
    with open(gtf_file) as fh:
        for line in fh:
            if line.startswith('#'): continue
            p = line.rstrip().split('\t')
            if len(p) < 9: continue
            seqid, ftype = p[0], p[2]
            start, end = int(p[3]), int(p[4])
            if ftype == 'gene': gene_by_seq[seqid].append((start, end))
            elif ftype == 'exon': exon_by_seq[seqid].append((start, end))
    return gene_by_seq, exon_by_seq

def compute_summary(scaffolds, sex_scaf, genome, by_seq_major, all_rep,
                    gene_by_seq, exon_by_seq):
    summary = {}
    for seqid in scaffolds:
        size = genome[seqid]
        d = {'size_mb': size / 1e6, 'is_sex': seqid == sex_scaf}
        d['rep_total'] = merge_bp(all_rep[seqid]) / size * 100
        for maj in MAJOR_ORDER:
            d[maj] = merge_bp(by_seq_major[seqid][maj]) / size * 100
        gene_ivs = gene_by_seq[seqid]
        exon_ivs = exon_by_seq[seqid]
        d['n_genes']    = len(gene_ivs)
        d['gene_dens']  = d['n_genes'] / (size / 1e6)
        d['gene_pct']   = merge_bp(gene_ivs) / size * 100
        d['exon_pct']   = merge_bp(exon_ivs) / size * 100
        d['intron_pct'] = max(0.0, d['gene_pct'] - d['exon_pct'])
        union_ivs = all_rep[seqid] + gene_ivs
        union_pct = merge_bp(union_ivs) / size * 100
        d['gene_only_pct']   = max(0.0, union_pct - d['rep_total'])
        d['unannotated_pct'] = max(0.0, 100.0 - union_pct)
        summary[seqid] = d
    return summary

def autosomal_avg(scaffolds, sex_scaf, genome, summary):
    autos      = [s for s in scaffolds if s != sex_scaf]
    total_size = sum(genome[s] for s in autos)
    keys = (['rep_total', 'gene_dens', 'gene_pct', 'exon_pct', 'intron_pct',
              'gene_only_pct', 'unannotated_pct'] + MAJOR_ORDER)
    return {k: sum(summary[s][k] * genome[s] for s in autos) / total_size
            for k in keys}

def compute_scatter(scaffolds, sex_scaf, genome, all_rep, gene_by_seq):
    auto_rep = []; auto_gene = []
    sex_rep  = []; sex_gene  = []
    for seqid in scaffolds:
        size = genome[seqid]
        positions = np.arange(0, size - WINDOW + 1, STEP, dtype=int)
        if len(positions) == 0: positions = np.array([0])
        gene_ivs      = gene_by_seq[seqid]
        merged_allrep = merge_sorted(all_rep[seqid])[1]
        merged_union  = merge_sorted(all_rep[seqid] + gene_ivs)[1]
        is_sex = seqid == sex_scaf
        trep = sex_rep if is_sex else auto_rep
        tgen = sex_gene if is_sex else auto_gene
        for ws in positions:
            we = min(int(ws) + WINDOW - 1, size - 1)
            rep_total  = cov_in_win(merged_allrep, ws, we)
            union_frac = cov_in_win(merged_union, ws, we)
            trep.append(rep_total * 100)
            tgen.append(max(0.0, union_frac - rep_total) * 100)
    return (np.array(auto_rep), np.array(auto_gene),
            np.array(sex_rep),  np.array(sex_gene))

# ── Load data ──────────────────────────────────────────────────────────────
print('Loading female genome...')
f_genome = load_genome(f'{FEMALE_DIR}/PhchiF_genome.fasta')
print('Loading female repeats...')
f_rep_major, f_all_rep = load_repeats(F_GFF)
print('Loading female genes...')
f_genes, f_exons = load_genes(f'{FEMALE_DIR}/PhchiF_gene_annotations.gtf')

print('Loading male genome...')
m_genome = load_genome(f'{MALE_DIR}/PhchiM_genome.fasta')
print('Loading male repeats...')
m_rep_major, m_all_rep = load_repeats(M_GFF)
print('Loading male genes...')
m_genes, m_exons = load_genes(f'{MALE_DIR}/PhchiM_gene_annotations.gtf')

print('Computing female summary stats...')
f_summary = compute_summary(F_SCAFFOLDS, F_SEX, f_genome,
                             f_rep_major, f_all_rep, f_genes, f_exons)
f_auto_avg = autosomal_avg(F_SCAFFOLDS, F_SEX, f_genome, f_summary)

print('Computing male summary stats...')
m_summary = compute_summary(M_SCAFFOLDS, M_SEX, m_genome,
                             m_rep_major, m_all_rep, m_genes, m_exons)
m_auto_avg = autosomal_avg(M_SCAFFOLDS, M_SEX, m_genome, m_summary)

print('Computing scatter windows...')
f_aut_rep, f_aut_gen, f_sex_rep, f_sex_gen = compute_scatter(
    F_SCAFFOLDS, F_SEX, f_genome, f_all_rep, f_genes)
m_aut_rep, m_aut_gen, m_sex_rep, m_sex_gen = compute_scatter(
    M_SCAFFOLDS, M_SEX, m_genome, m_all_rep, m_genes)

# ═══════════════════════════════════════════════════════════════════════════
# FIGURE
# ═══════════════════════════════════════════════════════════════════════════
print('Plotting...')

fig = plt.figure(figsize=(14, 15))
fig.patch.set_facecolor('white')

gs_top = gridspec.GridSpec(1, 2, figure=fig,
                           left=0.08, right=0.97, top=0.93, bottom=0.73,
                           wspace=0.32)
gs_mid = gridspec.GridSpec(1, 1, figure=fig,
                           left=0.08, right=0.97, top=0.67, bottom=0.46)
gs_bot = gridspec.GridSpec(1, 2, figure=fig,
                           left=0.08, right=0.97, top=0.40, bottom=0.06,
                           wspace=0.32)

ax_A = fig.add_subplot(gs_top[0, 0])
ax_B = fig.add_subplot(gs_top[0, 1])
ax_C = fig.add_subplot(gs_mid[0, 0])
ax_D = fig.add_subplot(gs_bot[0, 0])
ax_E = fig.add_subplot(gs_bot[0, 1])

def style_ax(ax, grid_axis='y'):
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_linewidth(0.6)
    ax.spines['bottom'].set_linewidth(0.6)
    if grid_axis:
        ax.grid(axis=grid_axis, linewidth=0.35, alpha=0.55, zorder=0)
    ax.set_axisbelow(True)

# x-positions: female at 0–6 (7 scaffolds), male at 8.5–13.5 (6 scaffolds)
f_x = np.arange(7, dtype=float)
m_x = np.arange(6, dtype=float) + 8.5
bw  = 0.65
GAP_X = 7.75

f_colors = [F_SEX_C if s == F_SEX else F_AUTO for s in F_SCAFFOLDS]
m_colors = [M_SEX_C if s == M_SEX else M_AUTO for s in M_SCAFFOLDS]

f_xlabels = ['S1', 'S2', 'S3', 'S4\n(U)', 'S5', 'S6', 'S7']
m_xlabels = ['S1', 'S2', 'S3', 'S4', 'S5', 'S6\n(V)']
all_x      = list(f_x) + list(m_x)
all_labels = f_xlabels + m_xlabels

# ── A: Gene density ──────────────────────────────────────────────────────────
f_gdens   = [f_summary[s]['gene_dens'] for s in F_SCAFFOLDS]
m_gdens   = [m_summary[s]['gene_dens'] for s in M_SCAFFOLDS]
all_gdens = f_gdens + m_gdens
ymax_a    = max(all_gdens) * 1.22

bars_af = ax_A.bar(f_x, f_gdens, bw, color=f_colors, linewidth=0, zorder=3)
bars_am = ax_A.bar(m_x, m_gdens, bw, color=m_colors, linewidth=0, zorder=3)

for bar, val in zip(list(bars_af) + list(bars_am), all_gdens):
    ax_A.text(bar.get_x() + bar.get_width() / 2,
              bar.get_height() + ymax_a * 0.02,
              f'{val:.0f}', ha='center', va='bottom', fontsize=6.5, color='#333333')

ax_A.axvline(GAP_X, color='#cccccc', linewidth=0.8, linestyle='--', zorder=0)
ax_A.set_xticks(all_x)
ax_A.set_xticklabels(all_labels, fontsize=7)
ax_A.set_xlim(-0.6, max(m_x) + 0.6)
ax_A.set_ylim(0, ymax_a)
ax_A.set_ylabel('Genes per Mb', fontsize=8)
style_ax(ax_A)
ax_A.set_title('Gene density per scaffold', fontsize=8.5, pad=5)
ax_A.text(np.mean(f_x), ymax_a * 0.99, 'Female (PhchiF)',
          ha='center', va='top', fontsize=7.5, fontweight='bold', color=F_AUTO)
ax_A.text(np.mean(m_x), ymax_a * 0.99, 'Male (PhchiM)',
          ha='center', va='top', fontsize=7.5, fontweight='bold', color=M_AUTO)
ax_A.text(-0.12, 1.10, 'A', transform=ax_A.transAxes,
          fontsize=11, fontweight='bold', va='top')

# ── B: Total repeat fraction ─────────────────────────────────────────────────
f_rep     = [f_summary[s]['rep_total'] for s in F_SCAFFOLDS]
m_rep     = [m_summary[s]['rep_total'] for s in M_SCAFFOLDS]
all_rep_v = f_rep + m_rep

bars_rf = ax_B.bar(f_x, f_rep, bw, color=f_colors, linewidth=0, zorder=3)
bars_rm = ax_B.bar(m_x, m_rep, bw, color=m_colors, linewidth=0, zorder=3)

for bar, val in zip(list(bars_rf) + list(bars_rm), all_rep_v):
    ax_B.text(bar.get_x() + bar.get_width() / 2,
              bar.get_height() + 1.5,
              f'{val:.0f}%', ha='center', va='bottom', fontsize=6.5, color='#333333')

ax_B.axvline(GAP_X, color='#cccccc', linewidth=0.8, linestyle='--', zorder=0)
ax_B.set_xticks(all_x)
ax_B.set_xticklabels(all_labels, fontsize=7)
ax_B.set_xlim(-0.6, max(m_x) + 0.6)
ax_B.set_ylim(0, 110)
ax_B.set_ylabel('Total repeat fraction (%)', fontsize=8)
style_ax(ax_B)
ax_B.set_title('Repeat content per scaffold', fontsize=8.5, pad=5)
ax_B.text(np.mean(f_x), 108, 'Female (PhchiF)',
          ha='center', va='top', fontsize=7.5, fontweight='bold', color=F_AUTO)
ax_B.text(np.mean(m_x), 108, 'Male (PhchiM)',
          ha='center', va='top', fontsize=7.5, fontweight='bold', color=M_AUTO)
ax_B.text(-0.12, 1.10, 'B', transform=ax_B.transAxes,
          fontsize=11, fontweight='bold', va='top')

# ── C: Scatter ────────────────────────────────────────────────────────────────
ax_C.scatter(f_aut_rep, f_aut_gen, c=F_AUTO, s=5, alpha=0.30, linewidths=0, zorder=3,
             label=f'Female autosomes  (n = {len(f_aut_rep):,} windows)')
ax_C.scatter(m_aut_rep, m_aut_gen, c=M_AUTO, s=5, alpha=0.30, linewidths=0, zorder=3,
             label=f'Male autosomes  (n = {len(m_aut_rep):,} windows)')
ax_C.scatter(f_sex_rep, f_sex_gen, c=F_SEX_C, s=9, alpha=0.55, linewidths=0, zorder=4,
             label=f'Female S4 / U chromosome  (n = {len(f_sex_rep):,} windows)')
ax_C.scatter(m_sex_rep, m_sex_gen, c=M_SEX_C, s=9, alpha=0.60, linewidths=0, zorder=4,
             label=f'Male S6 / V chromosome  (n = {len(m_sex_rep):,} windows)')

ax_C.set_xlabel('Repeat coverage per 200 kb window (%)', fontsize=8)
ax_C.set_ylabel('Gene-only coverage per 200 kb window (%)', fontsize=8)
ax_C.set_xlim(0, 100)
ax_C.set_ylim(0, None)
style_ax(ax_C, grid_axis=None)
ax_C.grid(linewidth=0.35, alpha=0.45, zorder=0)
ax_C.set_title('Gene vs. repeat coverage  (200 kb sliding windows, both genomes)',
               fontsize=8.5, pad=5)
ax_C.legend(fontsize=7.5, loc='upper right', frameon=False,
            handlelength=1.0, handletextpad=0.5, markerscale=2.0)
ax_C.text(-0.06, 1.05, 'C', transform=ax_C.transAxes,
          fontsize=11, fontweight='bold', va='top')

# ── D & E: 4-group comparison ─────────────────────────────────────────────────
d_groups = [
    ('Autosome\nmean (F)', f_auto_avg, F_AUTO),
    ('U chr\n(F-S4)',      f_summary[F_SEX], F_SEX_C),
    ('Autosome\nmean (M)', m_auto_avg, M_AUTO),
    ('V chr\n(M-S6)',      m_summary[M_SEX], M_SEX_C),
]
d_x  = np.arange(4)
bw_d = 0.65

# ── D: Repeat class composition ──────────────────────────────────────────────
bottoms_d = np.zeros(4)
for maj in MAJOR_ORDER:
    vals = np.array([g[maj] for _, g, _ in d_groups])
    ax_D.bar(d_x, vals, bw_d, bottom=bottoms_d,
             color=MAJOR_COLORS[maj], label=maj, linewidth=0, zorder=3)
    bottoms_d += vals

for i in range(4):
    unmasked = max(0.0, 100.0 - bottoms_d[i])
    ax_D.bar(i, unmasked, bw_d, bottom=bottoms_d[i],
             color=MAJOR_COLORS['Unmasked'], linewidth=0, zorder=3)

ax_D.axvline(1.5, color='#cccccc', linewidth=0.8, linestyle='--', zorder=0)
ax_D.set_xticks(d_x)
ax_D.set_xticklabels([l for l, _, _ in d_groups], fontsize=8)
ax_D.set_ylim(0, 100)
ax_D.set_ylabel('Genome fraction (%)', fontsize=8)
style_ax(ax_D)
ax_D.set_title('Repeat class composition', fontsize=8.5, pad=5)
ax_D.text(0.5, 103, 'Female', ha='center', va='bottom', fontsize=7.5,
          fontweight='bold', color=F_AUTO, clip_on=False)
ax_D.text(2.5, 103, 'Male', ha='center', va='bottom', fontsize=7.5,
          fontweight='bold', color=M_AUTO, clip_on=False)

for tick, (_, _, col) in zip(ax_D.get_xticklabels(), d_groups):
    tick.set_color(col)

d_legend = [mpatches.Patch(color=MAJOR_COLORS[m], label=m, linewidth=0)
            for m in MAJOR_ORDER]
d_legend.append(mpatches.Patch(color=MAJOR_COLORS['Unmasked'], label='Unmasked', linewidth=0))
ax_D.legend(handles=d_legend, fontsize=7, ncol=2, loc='upper left',
            frameon=False, handlelength=0.9, handletextpad=0.4,
            bbox_to_anchor=(0.01, 0.99))
ax_D.text(-0.14, 1.08, 'D', transform=ax_D.transAxes,
          fontsize=11, fontweight='bold', va='top')

# ── E: Genome compartment composition ────────────────────────────────────────
comp_defs = [
    ('Exons',          GENE_COLOR),
    ('Introns',        INTRON_COLOR),
    ('Unannotated',    UNANN_COLOR),
    ('LTR TEs',        MAJOR_COLORS['LTR']),
    ('DNA TEs',        MAJOR_COLORS['DNA']),
    ('LINE + MITE',    MAJOR_COLORS['LINE']),
    ('Unknown repeat', MAJOR_COLORS['Unknown']),
]

def get_comp(g):
    return [
        g['exon_pct'],
        g['intron_pct'],
        g['unannotated_pct'],
        g['LTR'],
        g['DNA'],
        g['LINE'] + g['MITE'],
        g['Unknown'],
    ]

bottoms_e = np.zeros(4)
for ci, (label, color) in enumerate(comp_defs):
    vals = np.array([get_comp(g)[ci] for _, g, _ in d_groups])
    ax_E.bar(d_x, vals, bw_d, bottom=bottoms_e,
             color=color, label=label, linewidth=0, zorder=3)
    bottoms_e += vals

ax_E.axvline(1.5, color='#cccccc', linewidth=0.8, linestyle='--', zorder=0)
ax_E.set_xticks(d_x)
ax_E.set_xticklabels([l for l, _, _ in d_groups], fontsize=8)
ax_E.set_ylim(0, 100)
ax_E.set_ylabel('Scaffold composition (%)', fontsize=8)
style_ax(ax_E)
ax_E.set_title('Genome compartment composition', fontsize=8.5, pad=5)
ax_E.text(0.5, 103, 'Female', ha='center', va='bottom', fontsize=7.5,
          fontweight='bold', color=F_AUTO, clip_on=False)
ax_E.text(2.5, 103, 'Male', ha='center', va='bottom', fontsize=7.5,
          fontweight='bold', color=M_AUTO, clip_on=False)

for tick, (_, _, col) in zip(ax_E.get_xticklabels(), d_groups):
    tick.set_color(col)

ax_E.legend(fontsize=7, ncol=2, loc='upper center', frameon=False,
            handlelength=0.9, handletextpad=0.4, columnspacing=0.8,
            bbox_to_anchor=(0.5, -0.14))
ax_E.text(-0.14, 1.08, 'E', transform=ax_E.transAxes,
          fontsize=11, fontweight='bold', va='top')

fig.text(0.5, 0.990,
         'Phaeomegaceros chiloensis  —  Female vs. male genome comparison',
         ha='center', va='top', fontsize=11, fontweight='bold', color='#1a1a1a')
fig.text(0.5, 0.974,
         'U chromosome (female S4)  ·  V chromosome (male S6)',
         ha='center', va='top', fontsize=8.5, color='#555555')

fig.savefig(f'{OUT}.pdf', dpi=300, bbox_inches='tight', facecolor='white')
fig.savefig(f'{OUT}.png', dpi=180, bbox_inches='tight', facecolor='white')
print(f'Saved {OUT}.pdf / .png')
