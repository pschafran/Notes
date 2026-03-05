"""
Whole-genome methylation comparison: PhphyF (female, Phymato_6) vs PhphyM (male, ref_genome)
Phymatoceros phymatodes

NOTE on methylation data:
  Female (PhphyF): PacBio HiFi + jasmine — 5mC only (18-column bedMethyl)
  Male   (PhphyM): Oxford Nanopore + modkit — 5mC + 5hmC (standard bedMethyl)
  Both formats have: chrom, start, end, mod_type, score, strand, ..., Nvalid_cov, pct_modified
  Only 5mC rows (col 3 == 'm') are used; 5hmC rows are skipped for both.

Steps:
  1. Extract autosomes → run minimap2 asm5 WGA (cached in wga/)
  2. Load 5mC methylation for all scaffolds
  3. Compute 200 kb windowed 5mC at syntenic positions (PAF-guided)
  4. Load repeat/gene landscape
  5. Per-site coverage and genomic context analysis
  6. Figure 1: Miami plot (female top, landscape middle, male bottom)
  7. Figure 2: Correlation scatter
  8. Figure 3: Coverage + context by methylation category

Sex chromosomes:
  U chromosome = PhphyF.S5 (4.9 Mb, female)
  V chromosome = PhphyM.S5 (4.0 Mb, male)
  Autosomes: S1-S4 in both genomes

Outputs (sex_chromosome_analyses/):
  PhphyF_PhphyM_methylation_manhattan.pdf/png
  PhphyF_PhphyM_methylation_manhattan_correlation.pdf/png
  PhphyF_PhphyM_methylation_context.pdf/png
  wga/PhphyF_PhphyM_autosomes.paf
"""

import bisect
import math
import os
import subprocess
import sys
from collections import defaultdict

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.patches as mpatches

# ── Global style ─────────────────────────────────────────────────────────────
plt.rcParams.update({
    'font.family':        'Liberation Sans',
    'font.size':          8,
    'axes.titlesize':     9,
    'axes.labelsize':     8,
    'xtick.labelsize':    8,
    'ytick.labelsize':    7,
    'legend.fontsize':    7.5,
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

# ── Paths ─────────────────────────────────────────────────────────────────────
BASE   = '/media/data/projects/hornwort_sex_chromosomes/analysis/Phymatoceros_phymatodes'
BASE_F = f'{BASE}/Phymato_6'
BASE_M = f'{BASE}/Phymato_ref_genome'

FASTA_F  = f'{BASE_F}/final_genome_prep/PhphyF_genome.chromosomes.fasta'
FASTA_M  = f'{BASE_M}/final_genome_prep/PhphyM_genome.chromosomes.fasta'
METHYL_F = f'{BASE_F}/methylation/PhphyF_genome.5mC.bed'
METHYL_M = f'{BASE_M}/methylation/PhphyM_genome.5mCG_5hmCG.bed'

GFF_F = f'{BASE_F}/final_genome_prep/PhphyF_repeat_annotations.gff'
GTF_F = f'{BASE_F}/final_genome_prep/braker_renamed.gtf'
GFF_M = f'{BASE_M}/final_genome_prep/pilon.3_renamed.fasta.mod.EDTA.TEanno.gff3'
GTF_M = f'{BASE_M}/final_genome_prep/braker_renamed.gtf'

MINIMAP2  = '/home/peter/bin/minimap2'
WGA_DIR   = 'wga'
PAF       = f'{WGA_DIR}/PhphyF_PhphyM_autosomes.paf'
AUTO_F_FA = f'{WGA_DIR}/PhphyF_autosomes.fasta'
AUTO_M_FA = f'{WGA_DIR}/PhphyM_autosomes.fasta'
OUT       = 'PhphyF_PhphyM_methylation_manhattan'

# ── Scaffold lists ────────────────────────────────────────────────────────────
F_AUTO   = ['PhphyF.S1', 'PhphyF.S2', 'PhphyF.S3', 'PhphyF.S4']
M_AUTO   = ['PhphyM.S1', 'PhphyM.S2', 'PhphyM.S3', 'PhphyM.S4']
F_SEXCHR = 'PhphyF.S5'   # U chromosome
M_SEXCHR = 'PhphyM.S5'   # V chromosome
F_ALL    = F_AUTO + [F_SEXCHR]
M_ALL    = M_AUTO + [M_SEXCHR]

SCAF_SIZES = {
    'PhphyF.S1': 48_759_091, 'PhphyF.S2': 44_430_696,
    'PhphyF.S3': 36_958_177, 'PhphyF.S4': 34_832_688,
    'PhphyF.S5':  4_914_197,
    'PhphyM.S1': 47_679_935, 'PhphyM.S2': 43_881_656,
    'PhphyM.S3': 35_704_750, 'PhphyM.S4': 34_574_324,
    'PhphyM.S5':  3_979_947,
}

WINDOW         = 200_000
METHYL_MIN_COV = 5
MIN_BLOCK_LEN  = 50_000
MIN_MAPQ       = 5

# ── Colours ───────────────────────────────────────────────────────────────────
F_COL   = '#3d6e9e'
M_COL   = '#3a7a52'
U_COLOR = '#b87070'
V_COLOR = '#c47a4a'

CHROM_COLS = {
    'PhphyF.S1': '#90b8d4', 'PhphyF.S2': '#6090b8',
    'PhphyF.S3': '#90c4a4', 'PhphyF.S4': '#5a9070',
}

MAJOR_ORDER  = ['LTR', 'DNA', 'LINE', 'MITE', 'Unknown']
MAJOR_COLORS = {
    'LTR':     '#c47a7a', 'DNA':     '#7aab8a',
    'LINE':    '#8a8ab8', 'MITE':    '#c8a86e',
    'Unknown': '#7aaabf',
}
GENE_COLOR   = '#5a87a5'
UNANN_COLOR  = '#dcdcdc'
STACK_COLORS = [GENE_COLOR, UNANN_COLOR] + [MAJOR_COLORS[m] for m in MAJOR_ORDER]
STACK_LABELS = ['Genes', 'Unannotated', 'LTR', 'DNA', 'LINE', 'MITE', 'Unknown repeat']

CATS    = ['unmethylated', 'intermediate', 'methylated']
CAT_COL = {'unmethylated': '#6a9ac4', 'intermediate': '#c4a86a', 'methylated': '#c46a6a'}
CTX_COLORS = {
    'intergenic': '#d8d8d8', 'gene_only': '#7aabcf',
    'gene+LTR':   '#9a7a7a', 'gene+DNA':   '#7aab8a',
    'gene+LINE':  '#8a8ab8', 'gene+MITE':  '#c8a86e',
    'gene+Unknown': '#7aaabf',
    'LTR':     '#c47a7a', 'DNA':     '#5a9070',
    'LINE':    '#6a6a98', 'MITE':    '#b08848',
    'Unknown': '#5a8a9f',
}

# ══════════════════════════════════════════════════════════════════════════════
# Helpers
# ══════════════════════════════════════════════════════════════════════════════

def pearsonr(x, y):
    xm, ym = x - np.mean(x), y - np.mean(y)
    r = float(np.dot(xm, ym) / np.sqrt(np.dot(xm, xm) * np.dot(ym, ym)))
    n = len(x)
    t = r * math.sqrt((n - 2) / max(1e-15, 1 - r ** 2))
    return r, float(math.erfc(abs(t) / math.sqrt(2)))

def spearmanr(x, y):
    return pearsonr(np.argsort(np.argsort(x)).astype(float),
                    np.argsort(np.argsort(y)).astype(float))

def linregress(x, y):
    xm, ym = np.mean(x), np.mean(y)
    sxy = float(np.sum((x - xm) * (y - ym)))
    sxx = float(np.sum((x - xm) ** 2))
    slope = sxy / sxx
    return slope, float(ym - slope * xm)

def merge_intervals(intervals):
    if not intervals: return []
    iv = sorted(intervals); m = [list(iv[0])]
    for s, e in iv[1:]:
        if s <= m[-1][1]: m[-1][1] = max(m[-1][1], e)
        else: m.append([s, e])
    return m

def cov_in_win(merged, ws, we):
    if not merged: return 0.0
    total = 0
    for s, e in merged:
        if e < ws: continue
        if s > we: break
        total += min(e, we) - max(s, ws) + 1
    return total / (we - ws + 1)

def get_major(cls):
    major = cls.split('/')[0] if '/' in cls else cls
    if major == 'TIR':      return 'DNA'
    if major == 'Penelope': return 'LINE'
    return major

def win_meth(pos_arr, met_arr, ws, we):
    if len(pos_arr) == 0: return float('nan')
    lo = bisect.bisect_left(pos_arr, ws)
    hi = bisect.bisect_right(pos_arr, we)
    if lo >= hi: return float('nan')
    return float(np.mean(met_arr[lo:hi]))

def load_genes(gtf_file, scaffolds):
    scaf_set = set(scaffolds); genes = defaultdict(list)
    with open(gtf_file) as fh:
        for line in fh:
            if line.startswith('#'): continue
            p = line.rstrip().split('\t')
            if len(p) < 5 or p[2] != 'gene' or p[0] not in scaf_set: continue
            genes[p[0]].append((int(p[3]), int(p[4])))
    for s in genes: genes[s].sort()
    return dict(genes)

def load_repeats_by_major(gff_file, scaffolds):
    scaf_set = set(scaffolds); reps = defaultdict(lambda: defaultdict(list))
    with open(gff_file) as fh:
        for line in fh:
            if line.startswith('#'): continue
            p = line.rstrip().split('\t')
            if len(p) < 9 or p[0] not in scaf_set: continue
            attrs = {k: v for a in p[8].split(';') if '=' in a
                     for k, v in [a.split('=', 1)]}
            reps[p[0]][get_major(attrs.get('Classification', 'Unknown'))].append(
                (int(p[3]), int(p[4])))
    for s in reps:
        for maj in reps[s]: reps[s][maj].sort()
    return dict(reps)

def label_sites(pos_arr, gene_ivs, rep_ivs_by_major):
    n = len(pos_arr)
    is_gene = np.zeros(n, dtype=bool)
    rep_cls = np.full(n, '', dtype=object)
    for major, ivs in rep_ivs_by_major.items():
        for s, e in ivs:
            lo = bisect.bisect_left(pos_arr, s)
            hi = bisect.bisect_right(pos_arr, e)
            if lo < hi: rep_cls[lo:hi] = major
    for s, e in gene_ivs:
        lo = bisect.bisect_left(pos_arr, s)
        hi = bisect.bisect_right(pos_arr, e)
        if lo < hi: is_gene[lo:hi] = True
    return is_gene, rep_cls

def categorize(pct_arr):
    cats = np.empty(len(pct_arr), dtype=object)
    cats[pct_arr < 20]                     = 'unmethylated'
    cats[(pct_arr >= 20) & (pct_arr < 80)] = 'intermediate'
    cats[pct_arr >= 80]                    = 'methylated'
    return cats

def context_fracs(mask, is_gene, rep_cls):
    in_rep = rep_cls != ''
    sub_g, sub_r, sub_cls = is_gene[mask], in_rep[mask], rep_cls[mask]
    n = mask.sum()
    if n == 0: return {}
    out = {'intergenic': float((~sub_g & ~sub_r).sum()) / n,
           'gene_only':  float((sub_g  & ~sub_r).sum()) / n}
    for maj in MAJOR_ORDER:
        out['gene+' + maj] = float((sub_g  & (sub_cls == maj)).sum()) / n
        out[maj]           = float((~sub_g & (sub_cls == maj)).sum()) / n
    return out

# ══════════════════════════════════════════════════════════════════════════════
# STEP 1 — Whole-genome alignment
# ══════════════════════════════════════════════════════════════════════════════
def write_fasta_subset(input_fa, output_fa, keep_seqs):
    keep, writing = set(keep_seqs), False
    with open(input_fa) as fi, open(output_fa, 'w') as fo:
        for line in fi:
            if line.startswith('>'):
                writing = line.split()[0][1:] in keep
            if writing: fo.write(line)

if not os.path.exists(PAF):
    print('Running whole-genome alignment...')
    os.makedirs(WGA_DIR, exist_ok=True)
    print('  Extracting female autosomes...')
    write_fasta_subset(FASTA_F, AUTO_F_FA, F_AUTO)
    print('  Extracting male autosomes...')
    write_fasta_subset(FASTA_M, AUTO_M_FA, M_AUTO)
    print('  Running minimap2 asm5...')
    result = subprocess.run(
        [MINIMAP2, '-x', 'asm5', '-t', '16', AUTO_F_FA, AUTO_M_FA, '-o', PAF],
        capture_output=True, text=True)
    if result.returncode != 0:
        print(f'minimap2 failed:\n{result.stderr}', file=sys.stderr); sys.exit(1)
    print(f'  Saved: {PAF}')
else:
    print(f'Using existing alignment: {PAF}')

# ══════════════════════════════════════════════════════════════════════════════
# STEP 2 — Parse PAF
# ══════════════════════════════════════════════════════════════════════════════
print('Parsing PAF...')
blocks = []
with open(PAF) as fh:
    for line in fh:
        p = line.split('\t')
        if len(p) < 12: continue
        q_name, q_start, q_end = p[0], int(p[2]), int(p[3])
        strand = p[4]
        r_name, r_start, r_end = p[5], int(p[7]), int(p[8])
        mapq = int(p[11])
        if mapq < MIN_MAPQ or (r_end - r_start) < MIN_BLOCK_LEN: continue
        if r_name not in F_AUTO or q_name not in M_AUTO: continue
        blocks.append((r_name, r_start, r_end, q_name, q_start, q_end, strand))
total_aln_bp = sum(r_end - r_start for _, r_start, r_end, *_ in blocks)
print(f'  {len(blocks)} blocks  |  {total_aln_bp/1e6:.1f} Mb aligned')

# ══════════════════════════════════════════════════════════════════════════════
# STEP 3 — Load methylation
# ══════════════════════════════════════════════════════════════════════════════
print('Loading methylation...')

def load_meth(bed_file, scaffolds, min_cov=METHYL_MIN_COV):
    raw = {s: ([], []) for s in scaffolds}; scset = set(scaffolds)
    with open(bed_file) as fh:
        for line in fh:
            p = line.split('\t')
            if len(p) < 11 or p[3] != 'm': continue
            if p[0] not in scset or int(p[9]) < min_cov: continue
            raw[p[0]][0].append(int(p[1]))
            raw[p[0]][1].append(float(p[10]))
    out = {}
    for s in scaffolds:
        pos = np.array(raw[s][0], dtype=np.int64)
        met = np.array(raw[s][1])
        idx = np.argsort(pos)
        out[s] = (pos[idx], met[idx])
        print(f'  {s}: {len(pos):,} sites  mean={np.mean(met):.1f}%' if len(met) else f'  {s}: 0 sites')
    return out

print('  Female (PacBio/jasmine 5mC)...')
meth_f = load_meth(METHYL_F, F_ALL)
print('  Male (ONT/modkit 5mC)...')
meth_m = load_meth(METHYL_M, M_ALL)

# ══════════════════════════════════════════════════════════════════════════════
# STEP 4 — Syntenic 200 kb windows
# ══════════════════════════════════════════════════════════════════════════════
print('Computing methylation at syntenic windows...')
syn_r_name, syn_f_center, syn_f_5mC, syn_m_5mC = [], [], [], []

for (r_name, r_start, r_end, q_name, q_start, q_end, strand) in blocks:
    r_len = r_end - r_start; q_len = q_end - q_start
    n_wins = max(1, r_len // WINDOW)
    for i in range(n_wins):
        fw_s = r_start + i * (r_len // n_wins)
        fw_e = min(r_start + (i + 1) * (r_len // n_wins), r_end)
        f_5mC = win_meth(*meth_f[r_name], fw_s, fw_e)
        frac_s, frac_e = i / n_wins, (i + 1) / n_wins
        if strand == '+':
            mw_s = q_start + int(frac_s * q_len); mw_e = q_start + int(frac_e * q_len)
        else:
            mw_s = q_end - int(frac_e * q_len);   mw_e = q_end - int(frac_s * q_len)
        m_5mC = win_meth(*meth_m[q_name], mw_s, mw_e)
        if not (np.isnan(f_5mC) or np.isnan(m_5mC)):
            syn_r_name.append(r_name); syn_f_center.append((fw_s + fw_e) / 2)
            syn_f_5mC.append(f_5mC);  syn_m_5mC.append(m_5mC)

syn_f_center = np.array(syn_f_center)
syn_f_5mC    = np.array(syn_f_5mC)
syn_m_5mC    = np.array(syn_m_5mC)
print(f'  {len(syn_f_5mC):,} syntenic windows with valid methylation')

# ══════════════════════════════════════════════════════════════════════════════
# STEP 5 — Sex chromosome + full autosomal windowed 5mC
# ══════════════════════════════════════════════════════════════════════════════
def scaffold_windows(meth_dict, scaf, size, window=WINDOW):
    pos_arr, met_arr = meth_dict[scaf]
    starts = np.arange(0, size, window, dtype=int)
    centers, vals = [], []
    for ws in starts:
        we = min(ws + window, size)
        centers.append((ws + we) / 2); vals.append(win_meth(pos_arr, met_arr, ws, we))
    return np.array(centers), np.array(vals)

u_centers, u_mC = scaffold_windows(meth_f, F_SEXCHR, SCAF_SIZES[F_SEXCHR])
v_centers, v_mC = scaffold_windows(meth_m, M_SEXCHR, SCAF_SIZES[M_SEXCHR])
f_auto_win = {s: scaffold_windows(meth_f, s, SCAF_SIZES[s]) for s in F_AUTO}

# ══════════════════════════════════════════════════════════════════════════════
# STEP 6 — x-axis layout
# ══════════════════════════════════════════════════════════════════════════════
GAP_BP = 3_000_000
f_auto_offset = {}; cum = 0
for s in F_AUTO: f_auto_offset[s] = cum; cum += SCAF_SIZES[s]
auto_total_bp = cum
u_offset_bp = auto_total_bp + GAP_BP
v_offset_bp = u_offset_bp + SCAF_SIZES[F_SEXCHR] + GAP_BP
total_bp    = v_offset_bp + SCAF_SIZES[M_SEXCHR]
total_mb    = total_bp / 1e6

syn_x_mb = np.array([(f_auto_offset[r] + c) / 1e6
                      for r, c in zip(syn_r_name, syn_f_center)])

f_auto_x_mb, f_auto_y = [], []
for s in F_AUTO:
    ctr, val = f_auto_win[s]
    for c, v in zip(ctr, val):
        if not np.isnan(v):
            f_auto_x_mb.append((f_auto_offset[s] + c) / 1e6); f_auto_y.append(v)
f_auto_x_mb = np.array(f_auto_x_mb); f_auto_y = np.array(f_auto_y)
u_x_mb = (u_offset_bp + u_centers) / 1e6
v_x_mb = (v_offset_bp + v_centers) / 1e6

# ══════════════════════════════════════════════════════════════════════════════
# STEP 7 — Summary statistics
# ══════════════════════════════════════════════════════════════════════════════
f_auto_mean = float(np.nanmean(syn_f_5mC))
m_auto_mean = float(np.nanmean(syn_m_5mC))
u_mean      = float(np.nanmean(u_mC[~np.isnan(u_mC)]))
v_mean      = float(np.nanmean(v_mC[~np.isnan(v_mC)]))
r_pearson,  p_pearson  = pearsonr(syn_f_5mC, syn_m_5mC)
r_spearman, p_spearman = spearmanr(syn_f_5mC, syn_m_5mC)

print(f'\nSummary statistics:')
print(f'  Autosomal 5mC: female={f_auto_mean:.1f}%  male={m_auto_mean:.1f}%')
print(f'  U chr 5mC: {u_mean:.1f}%   V chr 5mC: {v_mean:.1f}%')
print(f'  Pearson r={r_pearson:.4f}  Spearman ρ={r_spearman:.4f}  n={len(syn_f_5mC):,}')

# ══════════════════════════════════════════════════════════════════════════════
# STEP 8 — Landscape annotation data
# ══════════════════════════════════════════════════════════════════════════════
def load_landscape_windows(gff_file, gtf_file, scaffolds, scaf_sizes, window=WINDOW):
    scaf_set = set(scaffolds)
    print(f'  Parsing {os.path.basename(gff_file)}...')
    rep_by_seq = defaultdict(lambda: defaultdict(list))
    with open(gff_file) as fh:
        for line in fh:
            if line.startswith('#'): continue
            p = line.rstrip().split('\t')
            if len(p) < 9 or p[0] not in scaf_set: continue
            attrs = {k: v for a in p[8].split(';') if '=' in a
                     for k, v in [a.split('=', 1)]}
            rep_by_seq[p[0]][get_major(attrs.get('Classification', 'Unknown'))].append(
                (int(p[3]), int(p[4])))
    print(f'  Parsing {os.path.basename(gtf_file)}...')
    gene_by_seq = defaultdict(list)
    with open(gtf_file) as fh:
        for line in fh:
            if line.startswith('#'): continue
            p = line.rstrip().split('\t')
            if len(p) < 5 or p[2] != 'gene' or p[0] not in scaf_set: continue
            gene_by_seq[p[0]].append((int(p[3]), int(p[4])))

    win_data = {}
    for seqid in scaffolds:
        merged_rep   = {maj: merge_intervals(rep_by_seq[seqid][maj]) for maj in MAJOR_ORDER}
        all_rep_ivs  = [iv for maj in MAJOR_ORDER for iv in rep_by_seq[seqid][maj]]
        merged_all   = merge_intervals(all_rep_ivs)
        gene_ivs     = gene_by_seq[seqid]
        merged_union = merge_intervals(all_rep_ivs + list(gene_ivs))
        size         = scaf_sizes[seqid]
        positions    = np.arange(0, size, window, dtype=int)
        n            = len(positions)
        gene_only_arr   = np.zeros(n); unannotated_arr = np.zeros(n)
        rep_arr = {maj: np.zeros(n) for maj in MAJOR_ORDER}
        for i, ws in enumerate(positions):
            we         = min(ws + window - 1, size - 1)
            rep_total  = cov_in_win(merged_all,   ws, we)
            union_frac = cov_in_win(merged_union, ws, we)
            gene_only_arr[i]   = max(0.0, union_frac - rep_total)
            unannotated_arr[i] = max(0.0, 1.0 - union_frac)
            raw = {maj: cov_in_win(merged_rep[maj], ws, we) for maj in MAJOR_ORDER}
            raw_sum = sum(raw.values())
            if raw_sum > 1e-9:
                scale = rep_total / raw_sum
                for maj in MAJOR_ORDER: rep_arr[maj][i] = raw[maj] * scale
        win_data[seqid] = {
            'centers': (positions + window // 2).astype(float),
            'gene_only': gene_only_arr, 'unannotated': unannotated_arr, 'rep': rep_arr,
        }
    return win_data

print('\nLoading female landscape annotations...')
land_f = load_landscape_windows(GFF_F, GTF_F, F_ALL, SCAF_SIZES)
print('Loading male landscape annotations (V chr only)...')
land_m = load_landscape_windows(GFF_M, GTF_M, [M_SEXCHR], SCAF_SIZES)
print('Done loading annotations.')

# ══════════════════════════════════════════════════════════════════════════════
# STEP 9 — Per-site context analysis
# ══════════════════════════════════════════════════════════════════════════════
print('\nRunning per-site coverage and context analysis...')

def load_meth_full(bed_file, scaffolds, min_cov=METHYL_MIN_COV):
    raw = {s: ([], [], []) for s in scaffolds}; scset = set(scaffolds)
    with open(bed_file) as fh:
        for line in fh:
            p = line.split('\t')
            if len(p) < 11 or p[3] != 'm': continue
            if p[0] not in scset: continue
            cov = int(p[9])
            if cov < min_cov: continue
            raw[p[0]][0].append(int(p[1])); raw[p[0]][1].append(cov)
            raw[p[0]][2].append(float(p[10]))
    out = {}
    for s in scaffolds:
        pos = np.array(raw[s][0], dtype=np.int64)
        cov = np.array(raw[s][1], dtype=np.int32)
        pct = np.array(raw[s][2], dtype=np.float32)
        idx = np.argsort(pos); out[s] = (pos[idx], cov[idx], pct[idx])
    return out

ctx_results = {}
for glabel, scaffolds, methyl_file, gff_file, gtf_file in [
    ('Female (PhphyF)', F_AUTO, METHYL_F, GFF_F, GTF_F),
    ('Male (PhphyM)',   M_AUTO, METHYL_M, GFF_M, GTF_M),
]:
    print(f'\n  {glabel}')
    meth_full = load_meth_full(methyl_file, scaffolds)
    genes_ctx = load_genes(gtf_file, scaffolds)
    reps_ctx  = load_repeats_by_major(gff_file, scaffolds)

    pos_l, cov_l, pct_l, gene_l, rep_l = [], [], [], [], []
    for s in scaffolds:
        pos, cov, pct = meth_full[s]
        ig, rc = label_sites(pos, genes_ctx.get(s, []), reps_ctx.get(s, {}))
        pos_l.append(pos); cov_l.append(cov); pct_l.append(pct)
        gene_l.append(ig); rep_l.append(rc)

    all_cov  = np.concatenate(cov_l); all_pct  = np.concatenate(pct_l)
    all_gene = np.concatenate(gene_l); all_rep  = np.concatenate(rep_l)
    all_cats = categorize(all_pct); in_rep = all_rep != ''

    print(f'  Coverage by category (min_cov={METHYL_MIN_COV}):')
    print(f'  {"Category":16s}  {"N":>10}  {"Median":>7}  {"Mean":>7}  {"P10":>6}  {"P90":>6}')
    cov_data = {}
    for cat in CATS:
        mask = all_cats == cat; c = all_cov[mask]; cov_data[cat] = c
        if len(c) == 0: print(f'  {cat:16s}  {"0":>10}'); continue
        print(f'  {cat:16s}  {len(c):>10,}  {np.median(c):>7.1f}  {np.mean(c):>7.1f}  '
              f'{np.percentile(c,10):>6.1f}  {np.percentile(c,90):>6.1f}')

    print(f'  Condensed context summary:')
    print(f'  {"Category":16s}  {"N":>10}  {"intergenic":>11}  {"gene body":>10}  '
          f'{"in repeat":>10}  {"gene+repeat":>12}')
    ctx_data = {}
    for cat in CATS:
        mask = all_cats == cat; n = mask.sum()
        if n == 0: continue
        sub_g = all_gene[mask]; sub_r = in_rep[mask]
        print(f'  {cat:16s}  {n:>10,}  '
              f'{(~sub_g & ~sub_r).sum()/n*100:>10.1f}%  {sub_g.sum()/n*100:>9.1f}%  '
              f'{sub_r.sum()/n*100:>9.1f}%  {(sub_g & sub_r).sum()/n*100:>11.1f}%')
        ctx_data[cat] = context_fracs(mask, all_gene, all_rep)
    pct_data = {cat: all_pct[all_cats == cat] for cat in CATS}
    ctx_results[glabel] = {'cov_data': cov_data, 'ctx_data': ctx_data, 'pct_data': pct_data}

# ══════════════════════════════════════════════════════════════════════════════
# FIGURE 1 — Three-panel Miami + landscape
# ══════════════════════════════════════════════════════════════════════════════
print('\n\nPlotting Figure 1: three-panel Miami Manhattan...')

fig1 = plt.figure(figsize=(24, 12))
gs = gridspec.GridSpec(3, 1, figure=fig1, height_ratios=[2, 3, 2],
                       hspace=0.04, top=0.88, bottom=0.09, left=0.055, right=0.98)
ax_f = fig1.add_subplot(gs[0])
ax_land = fig1.add_subplot(gs[1], sharex=ax_f)
ax_m = fig1.add_subplot(gs[2], sharex=ax_f)

def shade_sections(ax):
    for i, s in enumerate(F_AUTO):
        col = '#f7f7f7' if i % 2 == 0 else '#eeeeee'
        ax.axvspan(f_auto_offset[s]/1e6, (f_auto_offset[s]+SCAF_SIZES[s])/1e6,
                   color=col, alpha=1.0, zorder=0)
    ax.axvspan(u_offset_bp/1e6, (u_offset_bp+SCAF_SIZES[F_SEXCHR])/1e6,
               color='#f9eded', alpha=1.0, zorder=0)
    ax.axvspan(v_offset_bp/1e6, (v_offset_bp+SCAF_SIZES[M_SEXCHR])/1e6,
               color='#fdf5ec', alpha=1.0, zorder=0)

shade_sections(ax_f); shade_sections(ax_land); shade_sections(ax_m)

sep_positions = [f_auto_offset[s]/1e6 for s in F_AUTO]
sep_positions += [(f_auto_offset[F_AUTO[-1]]+SCAF_SIZES[F_AUTO[-1]])/1e6,
                  u_offset_bp/1e6, (u_offset_bp+SCAF_SIZES[F_SEXCHR])/1e6,
                  v_offset_bp/1e6]
for ax in [ax_f, ax_land, ax_m]:
    for xp in sep_positions: ax.axvline(xp, color='#c0c0c0', lw=0.5, zorder=1)

# Top: female 5mC
ax_f.scatter(f_auto_x_mb, f_auto_y, c=F_COL, s=5, alpha=0.35, linewidths=0,
             zorder=3, rasterized=True)
ax_f.scatter(syn_x_mb, syn_f_5mC, c=F_COL, s=5, alpha=0.65, linewidths=0,
             zorder=4, rasterized=True)
valid_u = ~np.isnan(u_mC)
ax_f.scatter(u_x_mb[valid_u], u_mC[valid_u], c=U_COLOR, s=8, alpha=0.9,
             linewidths=0, zorder=5)
ax_f.plot([u_offset_bp/1e6, (u_offset_bp+SCAF_SIZES[F_SEXCHR])/1e6],
          [u_mean, u_mean], color=U_COLOR, lw=2.0, solid_capstyle='round', zorder=6)
xmax_auto_frac = auto_total_bp / total_bp
ax_f.axhline(f_auto_mean, color=F_COL, lw=1.0, ls='--', alpha=0.75, zorder=2,
             xmin=0, xmax=xmax_auto_frac)
ax_f.text(auto_total_bp*0.005/1e6, f_auto_mean+1.5,
          f'mean {f_auto_mean:.0f}%', color=F_COL, fontsize=7, va='bottom', style='italic')
ax_f.text((u_offset_bp+SCAF_SIZES[F_SEXCHR]/2)/1e6, u_mean+2, f'{u_mean:.0f}%',
          ha='center', va='bottom', fontsize=8, color=U_COLOR, fontweight='bold')
ax_f.set_ylim(0, 105); ax_f.set_yticks([0, 25, 50, 75, 100])
ax_f.set_yticklabels(['0','25','50','75','100'], fontsize=7)
ax_f.set_ylabel('5mC (%)', fontsize=9, labelpad=4)
ax_f.tick_params(axis='x', bottom=False, labelbottom=False)
for sp in ['bottom','top','right']: ax_f.spines[sp].set_visible(False)
ax_f.grid(axis='y', linewidth=0.3, alpha=0.5, zorder=0)
ax_f.text(-0.005, 0.85, '♀  female', transform=ax_f.transAxes,
          ha='right', va='center', fontsize=9.5, color=F_COL, style='italic')

# Middle: landscape
for seqid in F_ALL:
    wd  = land_f[seqid]
    off = f_auto_offset[seqid] if seqid in f_auto_offset else u_offset_bp
    x_mb = (off + wd['centers']) / 1e6
    layers = [wd['gene_only']*100, wd['unannotated']*100] + \
             [wd['rep'][m]*100 for m in MAJOR_ORDER]
    ax_land.stackplot(x_mb, layers, colors=STACK_COLORS, linewidth=0, zorder=2)
wd = land_m[M_SEXCHR]
x_mb = (v_offset_bp + wd['centers']) / 1e6
ax_land.stackplot(x_mb,
    [wd['gene_only']*100, wd['unannotated']*100] + [wd['rep'][m]*100 for m in MAJOR_ORDER],
    colors=STACK_COLORS, linewidth=0, zorder=2)
ax_land.set_ylim(0, 100); ax_land.set_yticks([0, 50, 100])
ax_land.set_yticklabels(['0','50','100'], fontsize=6)
ax_land.set_ylabel('% window', fontsize=7, labelpad=4)
ax_land.tick_params(axis='x', bottom=False, labelbottom=False)
for sp in ['top','bottom','right']: ax_land.spines[sp].set_visible(False)

# Bottom: male 5mC
ax_m.scatter(syn_x_mb, syn_m_5mC, c=M_COL, s=5, alpha=0.65, linewidths=0,
             zorder=4, rasterized=True)
valid_v = ~np.isnan(v_mC)
ax_m.scatter(v_x_mb[valid_v], v_mC[valid_v], c=V_COLOR, s=8, alpha=0.9,
             linewidths=0, zorder=5)
ax_m.plot([v_offset_bp/1e6, (v_offset_bp+SCAF_SIZES[M_SEXCHR])/1e6],
          [v_mean, v_mean], color=V_COLOR, lw=2.0, solid_capstyle='round', zorder=6)
ax_m.axhline(m_auto_mean, color=M_COL, lw=1.0, ls='--', alpha=0.75, zorder=2,
             xmin=0, xmax=xmax_auto_frac)
ax_m.text(auto_total_bp*0.005/1e6, m_auto_mean+1.5,
          f'mean {m_auto_mean:.0f}%', color=M_COL, fontsize=7, va='bottom', style='italic')
ax_m.text((v_offset_bp+SCAF_SIZES[M_SEXCHR]/2)/1e6, v_mean+2, f'{v_mean:.0f}%',
          ha='center', va='bottom', fontsize=8, color=V_COLOR, fontweight='bold')
ax_m.set_ylim(0, 105); ax_m.invert_yaxis()
ax_m.set_yticks([0, 25, 50, 75, 100])
ax_m.set_yticklabels(['0','25','50','75','100'], fontsize=7)
ax_m.set_ylabel('5mC (%)', fontsize=9, labelpad=4)
for sp in ['top','right']: ax_m.spines[sp].set_visible(False)
ax_m.grid(axis='y', linewidth=0.3, alpha=0.5, zorder=0)
ax_m.text(-0.005, 0.15, '♂  male', transform=ax_m.transAxes,
          ha='right', va='center', fontsize=9.5, color=M_COL, style='italic')

# x-axis labels
xtick_pos = [(f_auto_offset[s]+SCAF_SIZES[s]/2)/1e6 for s in F_AUTO]
xtick_lab  = [f'S{i+1}\n{SCAF_SIZES[s]/1e6:.0f} Mb' for i, s in enumerate(F_AUTO)]
xtick_pos += [(u_offset_bp+SCAF_SIZES[F_SEXCHR]/2)/1e6,
              (v_offset_bp+SCAF_SIZES[M_SEXCHR]/2)/1e6]
xtick_lab  += [f'U  (S5)\n{SCAF_SIZES[F_SEXCHR]/1e6:.1f} Mb',
               f'V  (S5)\n{SCAF_SIZES[M_SEXCHR]/1e6:.1f} Mb']
ax_m.set_xticks(xtick_pos); ax_m.set_xticklabels(xtick_lab, fontsize=8.5)
ax_m.tick_params(axis='x', length=0, pad=4)
ax_m.set_xlim(-1, total_mb + 1)

# Section labels
for txt, xpos, col in [
    ('Autosomes', auto_total_bp/2/1e6, '#333333'),
    ('U chr', (u_offset_bp+SCAF_SIZES[F_SEXCHR]/2)/1e6, U_COLOR),
    ('V chr', (v_offset_bp+SCAF_SIZES[M_SEXCHR]/2)/1e6, V_COLOR),
]:
    ax_f.text(xpos, 1.10, txt, transform=ax_f.get_xaxis_transform(),
              ha='center', va='bottom', fontsize=9.5, color=col)
ax_f.annotate('', xy=(0, 1.06), xytext=(auto_total_bp/1e6, 1.06),
              xycoords=ax_f.get_xaxis_transform(), textcoords=ax_f.get_xaxis_transform(),
              arrowprops=dict(arrowstyle='-', color='#888888', lw=0.8))

# Legends
land_patches = [mpatches.Patch(color=c, label=l, linewidth=0)
                for c, l in zip(STACK_COLORS, STACK_LABELS)]
meth_elems = [
    mpatches.Patch(color=F_COL,   label=f'♀ 5mC (autosomes, {f_auto_mean:.0f}% mean)  [PacBio]'),
    mpatches.Patch(color=M_COL,   label=f'♂ 5mC (autosomes, {m_auto_mean:.0f}% mean)  [ONT]'),
    mpatches.Patch(color=U_COLOR, label=f'♀ 5mC (U chr, {u_mean:.0f}% mean)'),
    mpatches.Patch(color=V_COLOR, label=f'♂ 5mC (V chr, {v_mean:.0f}% mean)'),
]
ax_land.legend(handles=land_patches, fontsize=6.5, ncol=4, loc='upper right',
               frameon=True, framealpha=0.9, edgecolor='none',
               handlelength=0.9, handletextpad=0.4, columnspacing=0.8, borderpad=0.4)
ax_f.legend(handles=meth_elems, fontsize=7.5, ncol=4, loc='lower center',
            frameon=False, handlelength=0.9, handletextpad=0.4,
            columnspacing=1.0, bbox_to_anchor=(0.5, 1.22))

fig1.text(0.5, 0.97,
          'Phymatoceros phymatodes — DNA methylation and genome composition (200 kb windows)',
          ha='center', va='top', fontsize=11, fontweight='bold', color='#1a1a1a')
fig1.text(0.5, 0.955,
          'Female (♀ PacBio/jasmine) above  ·  Male (♂ ONT/modkit) below  ·  '
          'Male autosomal windows at syntenic female position  ·  Nvalid_cov ≥ 5',
          ha='center', va='top', fontsize=8, color='#555555')

fig1.savefig(f'{OUT}.pdf', dpi=300, bbox_inches='tight', facecolor='white')
fig1.savefig(f'{OUT}.png', dpi=180, bbox_inches='tight', facecolor='white')
print(f'  Saved {OUT}.pdf/png')
plt.close(fig1)

# ══════════════════════════════════════════════════════════════════════════════
# FIGURE 2 — Correlation scatter
# ══════════════════════════════════════════════════════════════════════════════
print('Plotting Figure 2: correlation scatter...')
slope, intercept = linregress(syn_f_5mC, syn_m_5mC)
fig2, ax2 = plt.subplots(figsize=(7, 7)); fig2.patch.set_facecolor('white')
ax2.scatter(syn_f_5mC, syn_m_5mC,
            c=[CHROM_COLS[s] for s in syn_r_name], s=7, alpha=0.5,
            linewidths=0, zorder=3, rasterized=True)
ax2.plot([0,100],[0,100], color='#aaaaaa', ls='--', lw=1.0, zorder=1, label='1:1')
x_fit = np.array([0,100])
ax2.plot(x_fit, slope*x_fit+intercept, color='#333333', lw=1.5, zorder=4,
         label=f'Regression  (slope={slope:.2f}, intercept={intercept:.1f})')
ax2.text(0.97, 0.05,
         f'Pearson r = {r_pearson:.3f}  (p = {p_pearson:.1e})\n'
         f'Spearman ρ = {r_spearman:.3f}  (p = {p_spearman:.1e})\n'
         f'n = {len(syn_f_5mC):,} syntenic 200 kb windows',
         transform=ax2.transAxes, ha='right', va='bottom', fontsize=8,
         bbox=dict(boxstyle='round,pad=0.4', facecolor='white', edgecolor='#cccccc'))
ax2.set_xlabel('Female 5mC per 200 kb window (%) [PacBio]', fontsize=9)
ax2.set_ylabel('Male 5mC per 200 kb window (%) [ONT]', fontsize=9)
ax2.set_xlim(0,100); ax2.set_ylim(0,100); ax2.set_aspect('equal')
ax2.spines['top'].set_visible(False); ax2.spines['right'].set_visible(False)
ax2.legend(fontsize=7.5, loc='upper left', frameon=False)
chrom_patches = [mpatches.Patch(color=CHROM_COLS[s],
                                label=f'S{i+1}  ({SCAF_SIZES[s]/1e6:.0f} Mb)')
                 for i, s in enumerate(F_AUTO)]
ax2.add_artist(ax2.legend(handles=chrom_patches, fontsize=7.5, loc='lower right',
                           frameon=False, title='Chromosome', title_fontsize=7.5))
ax2.set_title('Methylation correlation at syntenic 200 kb windows\n'
              'Phymatoceros phymatodes  ·  autosomes only', fontsize=9, pad=8)
fig2.savefig(f'{OUT}_correlation.pdf', dpi=300, bbox_inches='tight', facecolor='white')
fig2.savefig(f'{OUT}_correlation.png', dpi=180, bbox_inches='tight', facecolor='white')
print(f'  Saved {OUT}_correlation.pdf/png')
plt.close(fig2)

# ══════════════════════════════════════════════════════════════════════════════
# FIGURE 3 — Coverage + context by methylation category
# ══════════════════════════════════════════════════════════════════════════════
print('Plotting Figure 3: context and coverage...')
ctx_stack_order  = (['intergenic','gene_only'] +
                    ['gene+'+m for m in MAJOR_ORDER] + list(MAJOR_ORDER))
ctx_stack_labels = (['Intergenic','Gene body only'] +
                    [f'Gene + {m}' for m in ['LTR','DNA','LINE','MITE','Unknown']] +
                    [f'{m} (repeat only)' for m in ['LTR','DNA','LINE','MITE','Unknown']])

fig3 = plt.figure(figsize=(14, 10))
gs_outer = gridspec.GridSpec(2, 1, figure=fig3, hspace=0.42,
                              top=0.93, bottom=0.07, left=0.07, right=0.97)
gs_cov = gridspec.GridSpecFromSubplotSpec(1, 2, subplot_spec=gs_outer[0], wspace=0.30)
gs_ctx = gridspec.GridSpecFromSubplotSpec(1, 2, subplot_spec=gs_outer[1], wspace=0.30)

for col, glabel in enumerate(['Female (PhphyF)', 'Male (PhphyM)']):
    ax = fig3.add_subplot(gs_cov[col])
    res = ctx_results[glabel]
    cats_present = [c for c in CATS if len(res['cov_data'].get(c, [])) > 0]
    vp = ax.violinplot([np.log10(res['cov_data'][c]) for c in cats_present],
                       positions=range(len(cats_present)), showmedians=True,
                       showextrema=False, widths=0.6)
    for body, cat in zip(vp['bodies'], cats_present):
        body.set_facecolor(CAT_COL[cat]); body.set_edgecolor('none'); body.set_alpha(0.75)
    vp['cmedians'].set_color('#222222'); vp['cmedians'].set_linewidth(1.5)
    for i, cat in enumerate(cats_present):
        med = np.median(res['cov_data'][cat])
        ax.text(i, np.log10(med)+0.08, f'{med:.0f}×',
                ha='center', va='bottom', fontsize=6.5, color='#333333')
    ax.set_xticks(range(len(cats_present)))
    ax.set_xticklabels([c.capitalize() for c in cats_present], fontsize=8)
    ax.set_ylabel('Coverage (log₁₀)', fontsize=8); ax.set_title(glabel, fontsize=9, pad=6)
    ax.spines['top'].set_visible(False); ax.spines['right'].set_visible(False)
    ax.set_yticks([1,1.5,2,2.5,3])
    ax.set_yticklabels([f'{10**y:.0f}' for y in [1,1.5,2,2.5,3]], fontsize=7)
    ax.axhline(np.log10(METHYL_MIN_COV), color='#aaaaaa', lw=0.8, ls='--')

fig3.text(0.5, 0.975, 'A  Coverage by methylation category',
          ha='center', va='top', fontsize=10, fontweight='bold')

all_handles_dict = {}  # label → patch; accumulated across both panels
for col, glabel in enumerate(['Female (PhphyF)', 'Male (PhphyM)']):
    ax = fig3.add_subplot(gs_ctx[col])
    res = ctx_results[glabel]
    x = np.arange(len(CATS)); bottoms = np.zeros(len(CATS))
    for key, disp in zip(ctx_stack_order, ctx_stack_labels):
        heights = np.array([res['ctx_data'].get(cat, {}).get(key, 0.0)*100 for cat in CATS])
        color   = CTX_COLORS.get(key, '#cccccc')
        ax.bar(x, heights, bottom=bottoms, color=color, width=0.55, linewidth=0)
        bottoms += heights
        if any(h > 0.5 for h in heights) and disp not in all_handles_dict:
            all_handles_dict[disp] = mpatches.Patch(color=color, label=disp, linewidth=0)
    ax.set_xticks(x); ax.set_xticklabels([c.capitalize() for c in CATS], fontsize=8)
    ax.set_ylabel('% of sites', fontsize=8); ax.set_ylim(0,100)
    ax.set_title(glabel, fontsize=9, pad=6)
    ax.spines['top'].set_visible(False); ax.spines['right'].set_visible(False)
# legend on male axes (ax = last iteration), covering both panels' categories
ax.legend(handles=list(all_handles_dict.values())[::-1], fontsize=6, loc='upper left',
          bbox_to_anchor=(1.02, 1.0), frameon=False,
          handlelength=0.9, handletextpad=0.4, labelspacing=0.35)

fig3.text(0.5, 0.495, 'B  Genomic context by methylation category',
          ha='center', va='top', fontsize=10, fontweight='bold')

ctx_out = 'PhphyF_PhphyM_methylation_context'
fig3.savefig(f'{ctx_out}.pdf', dpi=300, bbox_inches='tight', facecolor='white')
fig3.savefig(f'{ctx_out}.png', dpi=150, bbox_inches='tight', facecolor='white')
print(f'  Saved {ctx_out}.pdf/png')
plt.close(fig3)

# ══════════════════════════════════════════════════════════════════════════════
# STATS — Save context/coverage statistics to text file
# ══════════════════════════════════════════════════════════════════════════════
import datetime
print('\nSaving statistics...')

def _rep_key(maj):
    return 'gene+' + maj

_genome_labels = ['Female (PhphyF)', 'Male (PhphyM)']
stats_path = 'PhphyF_PhphyM_methylation_stats.txt'

with open(stats_path, 'w') as fout:
    fout.write('Phymatoceros phymatodes — per-site 5mC methylation statistics\n')
    fout.write(f'Generated: {datetime.date.today()}\n')
    fout.write(f'Min coverage: {METHYL_MIN_COV}×  |  Window: {WINDOW//1000} kb\n')
    fout.write('Note: Female = PacBio/jasmine (5mC only); Male = ONT/modkit (5mC + 5hmC, only 5mC used)\n')
    fout.write(f'\nSummary (200 kb autosomal windows):\n')
    fout.write(f'  Autosomal 5mC: female={f_auto_mean:.1f}%  male={m_auto_mean:.1f}%\n')
    fout.write(f'  U chr 5mC: {u_mean:.1f}%   V chr 5mC: {v_mean:.1f}%\n')
    fout.write(f'  Pearson r={r_pearson:.4f}  Spearman ρ={r_spearman:.4f}  n={len(syn_f_5mC):,} windows\n')
    for glabel in _genome_labels:
        res = ctx_results[glabel]
        total_n = sum(len(res['cov_data'].get(cat, [])) for cat in CATS)
        fout.write(f'\n{"="*60}\n{glabel}  (total autosomal sites: {total_n:,})\n{"="*60}\n')
        fout.write(f'\nSite counts by category:\n')
        fout.write(f'  {"Category":16s}  {"N":>10}  {"% total":>8}\n')
        for cat in CATS:
            n = len(res['cov_data'].get(cat, []))
            fout.write(f'  {cat:16s}  {n:>10,}  {n/total_n*100:>7.1f}%\n')
        fout.write(f'\nCoverage statistics:\n')
        fout.write(f'  {"Category":16s}  {"N":>10}  {"Median":>7}  {"Mean":>7}  {"P10":>6}  {"P90":>6}\n')
        for cat in CATS:
            c = res['cov_data'].get(cat, np.array([]))
            if len(c) == 0:
                fout.write(f'  {cat:16s}  {"0":>10}\n'); continue
            fout.write(f'  {cat:16s}  {len(c):>10,}  {np.median(c):>7.1f}  {np.mean(c):>7.1f}  '
                       f'{np.percentile(c,10):>6.1f}  {np.percentile(c,90):>6.1f}\n')
        fout.write(f'\nCondensed context (intergenic / gene body / in repeat / gene+repeat):\n')
        fout.write(f'  {"Category":16s}  {"N":>10}  {"intergenic":>11}  {"gene body":>10}  '
                   f'{"in repeat":>10}  {"gene+repeat":>12}\n')
        for cat in CATS:
            c = res['cov_data'].get(cat, np.array([]))
            n = len(c)
            if n == 0: continue
            ctx = res['ctx_data'].get(cat, {})
            fi = ctx.get('intergenic', 0) * 100
            fg = ctx.get('gene_only', 0) * 100
            fr = sum(ctx.get(maj, 0) for maj in MAJOR_ORDER) * 100
            fb = sum(ctx.get(_rep_key(maj), 0) for maj in MAJOR_ORDER) * 100
            fout.write(f'  {cat:16s}  {n:>10,}  {fi:>10.1f}%  {fg:>9.1f}%  {fr:>9.1f}%  {fb:>11.1f}%\n')
        fout.write(f'\nRepeat context by major class (repeat-only + gene+repeat combined):\n')
        fout.write(f'  {"Category":16s}' + ''.join(f'  {m:>14}' for m in MAJOR_ORDER) + '\n')
        for cat in CATS:
            c = res['cov_data'].get(cat, np.array([]))
            if len(c) == 0: continue
            ctx = res['ctx_data'].get(cat, {})
            row = f'  {cat:16s}'
            for maj in MAJOR_ORDER:
                tot = (ctx.get(maj, 0) + ctx.get(_rep_key(maj), 0)) * 100
                row += f'  {tot:>13.1f}%'
            fout.write(row + '\n')
print(f'  Saved {stats_path}')

# ══════════════════════════════════════════════════════════════════════════════
# FIGURE 4 — Per-site methylation distribution histogram
# ══════════════════════════════════════════════════════════════════════════════
print('Plotting Figure 4: methylation distribution...')

fig4, axes4 = plt.subplots(1, 2, figsize=(12, 5))
fig4.patch.set_facecolor('white')
bins = np.arange(0, 102, 2)
bin_centers = (bins[:-1] + bins[1:]) / 2

for col, glabel in enumerate(_genome_labels):
    ax = axes4[col]
    res = ctx_results[glabel]
    total_n = sum(len(res['pct_data'].get(cat, [])) for cat in CATS)
    for cat in CATS:
        pct = res['pct_data'].get(cat, np.array([]))
        if len(pct) == 0: continue
        counts, _ = np.histogram(pct, bins=bins)
        ax.bar(bin_centers, counts / total_n * 100, width=2.0,
               color=CAT_COL[cat], alpha=0.82, linewidth=0,
               label=f'{cat.capitalize()}  ({len(pct)/total_n*100:.1f}%,  n={len(pct):,})')
    ax.axvline(20, color='#666666', lw=0.8, ls='--', alpha=0.6, zorder=5)
    ax.axvline(80, color='#666666', lw=0.8, ls='--', alpha=0.6, zorder=5)
    ax.set_xlabel('% modified (5mC per CpG site)', fontsize=9)
    ax.set_ylabel('% of autosomal CpG sites', fontsize=9)
    ax.set_xlim(0, 100)
    ax.set_title(glabel, fontsize=10, pad=6)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.legend(fontsize=8, frameon=False, loc='upper center')

tech_note = '  Female: PacBio/jasmine [5mC only]    Male: ONT/modkit [5mC only used]'
fig4.suptitle('Phymatoceros phymatodes — autosomal CpG methylation distribution\n'
              f'per-site % modified  ·  min_cov ≥ 5×  ·  autosomes only\n{tech_note}',
              fontsize=10, fontweight='bold')
plt.tight_layout()
dist_out = 'PhphyF_PhphyM_methylation_distribution'
fig4.savefig(f'{dist_out}.pdf', dpi=300, bbox_inches='tight', facecolor='white')
fig4.savefig(f'{dist_out}.png', dpi=150, bbox_inches='tight', facecolor='white')
print(f'  Saved {dist_out}.pdf/png')
plt.close(fig4)

print('\nDone.')
