#!/usr/bin/env python3
"""
Flatten CoGe chain-encoded GFF into standard GFF3.

CoGe encodes some genes (V1A-suffix IDs) as a linked-list of mRNA features
where each node contains one exon/CDS, and each node's Parent is the
previous mRNA (not the gene). This script converts them to standard GFF3:
  gene → one mRNA → multiple exon/CDS children

Plain genes (non-V1A) are passed through unchanged.
"""

import sys
import re

GFF_IN  = "Syntrichia_ruralis_annos1-cds0-id_typename-nu1-upa1-add_chr0.gid66749_renamed.gff"
GFF_OUT = "Syntrichia_ruralis_annos1-cds0-id_typename-nu1-upa1-add_chr0.gid66749_renamed.fixed.gff"

def parse_attrs(col9):
    attrs = {}
    for field in col9.split(';'):
        field = field.strip()
        if '=' in field:
            k, v = field.split('=', 1)
            attrs[k] = v
    return attrs

# ── Pass 1: load everything into memory ────────────────────────────────────
genes   = {}   # id → raw line + parsed fields
mrnas   = {}   # id → {line, fields, parent, children_cds, children_exon}
cdss    = []   # list of (parent_mrna_id, raw_line, fields)
exons   = []   # list of (parent_mrna_id, raw_line, fields)
others  = []   # comment lines, other feature types

with open(GFF_IN) as fh:
    for raw in fh:
        line = raw.rstrip('\n')
        if line.startswith('#') or line.strip() == '':
            others.append(('comment', line))
            continue
        fields = line.split('\t')
        if len(fields) < 9:
            others.append(('comment', line))
            continue
        ftype = fields[2]
        attrs = parse_attrs(fields[8])
        fid   = attrs.get('ID', '')
        par   = attrs.get('Parent', '')

        if ftype == 'gene':
            genes[fid] = {'line': line, 'fields': fields, 'attrs': attrs}
        elif ftype == 'mRNA':
            mrnas[fid] = {'line': line, 'fields': fields, 'attrs': attrs,
                          'parent': par, 'exons': [], 'cdss': []}
        elif ftype == 'CDS':
            cdss.append((par, line, fields, attrs))
        elif ftype == 'exon':
            exons.append((par, line, fields, attrs))
        else:
            others.append((ftype, line))

# Attach CDS and exon features to their immediate mRNA parent
for par, line, fields, attrs in cdss:
    if par in mrnas:
        mrnas[par]['cdss'].append((line, fields, attrs))

for par, line, fields, attrs in exons:
    if par in mrnas:
        mrnas[par]['exons'].append((line, fields, attrs))

# ── Pass 2: identify chain structure ───────────────────────────────────────
# mRNA → next mRNA in chain (its single mRNA child)
mrna_child = {}   # parent_mrna_id → child_mrna_id
for mid, m in mrnas.items():
    par = m['parent']
    if par in mrnas:           # parent is an mRNA, not a gene
        mrna_child[par] = mid

# For each gene, find its direct mRNA child
gene_mrna = {}   # gene_id → first mRNA ID (direct child)
for mid, m in mrnas.items():
    par = m['parent']
    if par in genes:
        gene_mrna[par] = mid

# ── Pass 3: write output ────────────────────────────────────────────────────
v1a_written = 0
plain_written = 0
chain_total_exons = 0

with open(GFF_OUT, 'w') as out:
    # Write comments first
    for ftype, line in others:
        if ftype == 'comment':
            out.write(line + '\n')

    for gid, gene in genes.items():
        gfields = gene['fields']
        out.write(gene['line'] + '\n')

        first_mrna_id = gene_mrna.get(gid)
        if first_mrna_id is None:
            continue   # gene with no mRNA — skip

        is_v1a = gid.endswith('_V1A')

        if not is_v1a:
            # ── Plain gene: write mRNA then its exons/CDS unchanged ──────
            m = mrnas[first_mrna_id]
            out.write(m['line'] + '\n')
            for line, fields, attrs in m['exons']:
                out.write(line + '\n')
            for line, fields, attrs in m['cdss']:
                out.write(line + '\n')
            plain_written += 1

        else:
            # ── V1A gene: follow chain, collect all exons/CDS ─────────────
            # Build synthetic mRNA spanning the gene coordinates
            gseq   = gfields[0]
            gsrc   = gfields[1]
            gstart = gfields[3]
            gend   = gfields[4]
            gstrand= gfields[6]
            gphase = gfields[7]

            # Use the gene ID minus _V1A suffix as base for mRNA ID
            base = re.sub(r'_V1A$', '', gid)
            mrna_id = base + '.mRNA1'

            mrna_attrs = (f"ID={mrna_id};Name={mrna_id};"
                          f"Parent={gid};gene={base}")
            mrna_line = '\t'.join([gseq, gsrc, 'mRNA',
                                   gstart, gend, '.', gstrand, '.', mrna_attrs])
            out.write(mrna_line + '\n')

            # Walk the chain
            chain_exons = []
            chain_cdss  = []
            cur = first_mrna_id
            n = 0
            while cur is not None:
                m = mrnas[cur]
                chain_exons.extend(m['exons'])
                chain_cdss.extend(m['cdss'])
                cur = mrna_child.get(cur)
                n += 1

            chain_total_exons += n

            # Re-parent all exons and CDS to the new mRNA ID
            exon_num = 1
            for line, fields, attrs in chain_exons:
                attrs['Parent'] = mrna_id
                attrs['ID']     = f"{mrna_id}.exon{exon_num}"
                new_attrs = ';'.join(f"{k}={v}" for k, v in attrs.items())
                out.write('\t'.join(fields[:8] + [new_attrs]) + '\n')
                exon_num += 1

            cds_num = 1
            for line, fields, attrs in chain_cdss:
                attrs['Parent'] = mrna_id
                attrs['ID']     = f"{mrna_id}.CDS{cds_num}"
                new_attrs = ';'.join(f"{k}={v}" for k, v in attrs.items())
                out.write('\t'.join(fields[:8] + [new_attrs]) + '\n')
                cds_num += 1

            v1a_written += 1

print(f"Plain genes written:  {plain_written}")
print(f"V1A genes flattened:  {v1a_written}")
print(f"Total chain nodes:    {chain_total_exons}  (avg {chain_total_exons/v1a_written:.1f} exons/gene)")
print(f"Output: {GFF_OUT}")
