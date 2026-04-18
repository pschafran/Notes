#!/usr/bin/env python3
"""
Plot Ks vs. chromosomal position for U-V gametolog pairs, coloured by TE status.

For each species with a GTF (Papro, Phchi, Phphy, Ledus):
  - Two panels per species: U gene position (top) and V gene position (bottom)
  - x = gene midpoint on its chromosome; y = Ks
  - TE-flagged genes highlighted
"""

import os
import re
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# ── Config ─────────────────────────────────────────────────────────────────────

KS_MIN = 0.05
KS_MAX = 5.0

BASE  = '/media/data/projects/hornwort_sex_chromosomes/analysis/ks'
ANNO  = '/media/data/projects/hornwort_sex_chromosomes/analysis/functional_annotations'

SPECIES = {
    'Papro': {'u': 'Papro.U',  'v': 'Papro.V',  'label': 'P. proskaueri',
              'gtf': 'Papro_UV_genome.genes.gtf',
              'anno_f': 'PaproF_PROT.emapper.annotations',
              'anno_m': 'PaproM_PROT.emapper.annotations',
              'prefix_f': 'PaproF', 'prefix_m': 'PaproM',
              'scaf_f': ('S5', 'U'), 'scaf_m': ('S5', 'V')},
    'Phchi': {'u': 'Phchi.U',  'v': 'Phchi.V',  'label': 'P. chiloensis',
              'gtf': 'Phchi_UV_genome.genes.gtf',
              'anno_f': 'PhchiF_PROT.emapper.annotations',
              'anno_m': 'PhchiM_PROT.emapper.annotations',
              'prefix_f': 'PhchiF', 'prefix_m': 'PhchiM',
              'scaf_f': ('S4', 'U'), 'scaf_m': ('S6', 'V')},
    'Phphy': {'u': 'Phphy.U',  'v': 'Phphy.V',  'label': 'P. phymatodes',
              'gtf': 'Phphy_UV_genome.genes.gtf',
              'anno_f': 'PhphyF_PROT.emapper.annotations',
              'anno_m': 'PhphyM_PROT.emapper.annotations',
              'prefix_f': 'PhphyF', 'prefix_m': 'PhphyM',
              'scaf_f': ('S5', 'U'), 'scaf_m': ('S5', 'V')},
    'Ledus': {'u': 'Ledus.U',  'v': 'Ledus.V',  'label': 'L. dussii',
              'gtf': 'Ledus_UV_genome.genes.gtf',
              'anno_f': 'LedusF_PROT.emapper.annotations',
              'anno_m': 'LedusM_PROT.emapper.annotations',
              'prefix_f': 'LedusF', 'prefix_m': 'LedusM',
              'scaf_f': ('S3', 'U'), 'scaf_m': ('S5', 'V')},
}

TE_PATTERN = re.compile(
    r'transpos|retrotranspos|reverse transcriptase|integrase|'
    r'gypsy|copia|helitron|mule|mudr|piggybac|gag.{0,20}(ltr|retro|transpos)',
    re.IGNORECASE
)
TE_EXCLUDE = re.compile(r'telomerase', re.IGNORECASE)

COLOR_NORMAL = '#5B8DB8'
COLOR_TE     = '#E07B54'

# ── Loaders ────────────────────────────────────────────────────────────────────

def strip_transcript(name):
    return re.sub(r'\.t\d+$', '', name)


def load_positions(gtf_path):
    """Return dict gene_id -> midpoint (bp)."""
    pos = {}
    with open(gtf_path) as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.strip().split('\t')
            if len(parts) < 9 or parts[2] != 'gene':
                continue
            gene_id = parts[8].strip()
            mid = (int(parts[3]) + int(parts[4])) / 2
            pos[gene_id] = mid
    return pos


def load_te_genes(anno_path, genotype_prefix, merged_prefix, scaf_rename):
    """Return set of merged-genome gene IDs flagged as TE-related.

    scaf_rename: (old_scaf, new_scaf) e.g. ('S5', 'U') to map
    {prefix}.S5G... -> {merged}.UG...
    """
    old_scaf, new_scaf = scaf_rename
    # Build the exact string to replace in the gene ID, e.g. 'PaproF.S5G' -> 'Papro.UG'
    old_frag = f'{genotype_prefix}.{old_scaf}G'
    new_frag = f'{merged_prefix}.{new_scaf}G'
    te_genes = set()
    with open(anno_path) as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.strip().split('\t')
            if len(parts) < 8:
                continue
            desc = parts[7]
            if TE_PATTERN.search(desc) and not TE_EXCLUDE.search(desc):
                query = parts[0]  # e.g. PaproF.S5G038800.t1
                merged_id = query.replace(old_frag, new_frag, 1)
                # Fallback for autosomal genes: just replace species prefix
                if merged_id == query:
                    merged_id = query.replace(genotype_prefix, merged_prefix, 1)
                te_genes.add(strip_transcript(merged_id))
    return te_genes


def load_uv_ks(sp, cfg):
    """Return DataFrame of U-V gametolog pairs with dS values."""
    tsv = os.path.join(BASE, f'wgd_ksd_{sp}', f'gametologs_{sp}.tsv.ks.tsv')
    df = pd.read_csv(tsv, sep='\t')
    u, v = cfg['u'], cfg['v']
    uv_mask = (
        (df['gene1'].str.startswith(u) & df['gene2'].str.startswith(v)) |
        (df['gene1'].str.startswith(v) & df['gene2'].str.startswith(u))
    )
    df = df[uv_mask & df['dS'].notna()].copy()
    df = df[(df['dS'] >= KS_MIN) & (df['dS'] <= KS_MAX)]
    # Normalise so u_gene is always the U copy
    swap = df['gene1'].str.startswith(v)
    df.loc[swap, ['gene1', 'gene2']] = df.loc[swap, ['gene2', 'gene1']].values
    df = df.rename(columns={'gene1': 'u_gene', 'gene2': 'v_gene'})
    df['u_id'] = df['u_gene'].apply(strip_transcript)
    df['v_id'] = df['v_gene'].apply(strip_transcript)
    return df


# ── Plot ───────────────────────────────────────────────────────────────────────

fig, axes = plt.subplots(2, 4, figsize=(14, 6))
fig.subplots_adjust(hspace=0.5, wspace=0.35)

for col, (sp, cfg) in enumerate(SPECIES.items()):
    positions = load_positions(os.path.join(BASE, cfg['gtf']))
    te_f = load_te_genes(os.path.join(ANNO, cfg['anno_f']), cfg['prefix_f'], sp, cfg['scaf_f'])
    te_m = load_te_genes(os.path.join(ANNO, cfg['anno_m']), cfg['prefix_m'], sp, cfg['scaf_m'])
    te_all = te_f | te_m

    df = load_uv_ks(sp, cfg)
    df['u_pos'] = df['u_id'].map(positions)
    df['v_pos'] = df['v_id'].map(positions)
    df['is_te'] = df['u_id'].isin(te_all) | df['v_id'].isin(te_all)

    n_te = df['is_te'].sum()
    n_total = len(df)

    for row, (pos_col, chrom_label) in enumerate([('u_pos', 'U'), ('v_pos', 'V')]):
        ax = axes[row, col]
        sub = df.dropna(subset=[pos_col])
        normal = sub[~sub['is_te']]
        te     = sub[sub['is_te']]

        ax.scatter(normal[pos_col] / 1e6, normal['dS'],
                   c=COLOR_NORMAL, s=12, alpha=0.6, linewidths=0, label='Non-TE')
        ax.scatter(te[pos_col] / 1e6, te['dS'],
                   c=COLOR_TE, s=16, alpha=0.8, linewidths=0, label='TE-related')

        ax.set_ylim(0, KS_MAX)
        ax.set_xlabel('Position (Mb)', fontsize=8)
        ax.set_ylabel('Ks', fontsize=8)
        ax.tick_params(labelsize=7)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)

        title = f'{cfg["label"]}' if row == 0 else ''
        if title:
            ax.set_title(title, fontsize=9, style='italic')
        ax.text(0.02, 0.97, chrom_label, transform=ax.transAxes,
                ha='left', va='top', fontsize=9, fontweight='bold')

    # TE count annotation on top panel
    axes[0, col].text(0.98, 0.97, f'TE: {n_te}/{n_total}',
                      transform=axes[0, col].transAxes,
                      ha='right', va='top', fontsize=7, color='#555555')

# ── Legend ─────────────────────────────────────────────────────────────────────
from matplotlib.lines import Line2D
legend_elements = [
    Line2D([0], [0], marker='o', color='w', markerfacecolor=COLOR_NORMAL,
           markersize=6, label='Non-TE gametolog'),
    Line2D([0], [0], marker='o', color='w', markerfacecolor=COLOR_TE,
           markersize=6, label='TE-related gametolog'),
]
fig.legend(handles=legend_elements, loc='lower center', ncol=2,
           fontsize=8, frameon=False, bbox_to_anchor=(0.5, -0.03))

plt.rcParams['font.family'] = 'Liberation Sans'

for fmt in ('pdf', 'png'):
    fig.savefig(os.path.join(BASE, f'ks_vs_position.{fmt}'),
                dpi=300, bbox_inches='tight')

print("Saved ks_vs_position.pdf and ks_vs_position.png")
