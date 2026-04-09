import pandas as pd
import os
import glob

OUTDIR = '/media/data/resources/marpolbase_expression_data/repro_coexp'
TMP    = '/media/data/resources/marpolbase_expression_data/repro_coexp/tmp'

# Load node community assignments
nodes = pd.read_csv(f'{OUTDIR}/node_community_assignments.tsv', sep='\t', index_col='gene_id')

# Parse tmp file — tab-separated, columns by position
# Col 1: Category, Col 2: OG, Col 8 (0-indexed): Mp orthologs
tmp = pd.read_csv(TMP, sep='\t', header=0, dtype=str)

# Column names from the file
print("Columns:", list(tmp.columns[:12]))

# Remove old OG assignment files
old = glob.glob(f'{OUTDIR}/OG*_community_assignment.tsv')
for f in old:
    os.remove(f)
print(f"Removed {len(old)} old OG assignment files")

og_col  = tmp.columns[1]   # 'OG'
mp_col  = tmp.columns[7]   # 'Mp orthologs'

written = 0
skipped_no_mp = 0
skipped_not_in_network = 0

for _, row in tmp.iterrows():
    og  = str(row[og_col]).strip()
    mps = str(row[mp_col]).strip()

    if not og or og in ('nan', '—', ''):
        continue
    if mps in ('nan', '—', ''):
        skipped_no_mp += 1
        print(f"  {og}: no Mp orthologs — skipping")
        continue

    # Parse gene IDs, strip isoform suffix (.1, .2 etc.)
    gene_ids = []
    for g in mps.split(','):
        g = g.strip()
        if g and g != '—':
            gene_ids.append(g.split('.')[0])

    # Look up in network
    found = [g for g in gene_ids if g in nodes.index]
    missing = [g for g in gene_ids if g not in nodes.index]

    if missing:
        print(f"  {og}: not in network — {missing}")

    if not found:
        skipped_not_in_network += 1
        print(f"  {og}: no Mp genes in network — skipping file")
        continue

    subset = nodes.loc[found].reset_index()
    outpath = f'{OUTDIR}/{og}_community_assignment.tsv'
    subset.to_csv(outpath, sep='\t', index=False)
    print(f"  {og}: wrote {len(found)} gene(s) → {os.path.basename(outpath)}")
    written += 1

print(f"\nDone. Files written: {written}  |  No Mp genes: {skipped_no_mp}  |  None in network: {skipped_not_in_network}")
