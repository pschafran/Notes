#!/usr/bin/env python3
"""
split_mt_per_sample.py
Split phylogeny_paper/mt/all_mt_scaffolds.fasta into one file per sample.
Headers in source are already formatted as sample_id|scaffold_N.
"""

from pathlib import Path

MTDIR = Path('/media/data/projects/iva_phylogeny/analysis/phylogeny_paper/mt')
combined = MTDIR / 'all_mt_scaffolds.fasta'

samples = {}
current_header, current_seq = None, []

with open(combined) as f:
    for line in f:
        line = line.rstrip()
        if line.startswith('>'):
            if current_header is not None:
                sid = current_header.split('|')[0]
                samples.setdefault(sid, []).append((current_header, ''.join(current_seq)))
            current_header = line[1:]
            current_seq = []
        else:
            current_seq.append(line)
    if current_header is not None:
        sid = current_header.split('|')[0]
        samples.setdefault(sid, []).append((current_header, ''.join(current_seq)))

rows = [('sample_id', 'n_scaffolds', 'total_bp', 'largest_scaffold_bp')]
for sid in sorted(samples):
    recs = samples[sid]
    out_path = MTDIR / f'{sid}_mt.fasta'
    with open(out_path, 'w') as f:
        for header, seq in recs:
            f.write(f'>{header}\n')
            for i in range(0, len(seq), 80):
                f.write(seq[i:i+80] + '\n')
    sizes = sorted([len(s) for _, s in recs], reverse=True)
    rows.append((sid, str(len(recs)), str(sum(sizes)), str(sizes[0])))
    print(f'  {sid}: {len(recs)} scaffolds, {sum(sizes):,} bp total')

with open(MTDIR / 'mt_sample_manifest.tsv', 'w') as f:
    for row in rows:
        f.write('\t'.join(row) + '\n')

print(f'\n{len(samples)} sample files written.')
