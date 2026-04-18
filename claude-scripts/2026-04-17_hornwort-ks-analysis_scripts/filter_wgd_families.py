#!/usr/bin/env python3
"""
Filter wgd dmd family TSVs for Ks analysis.

Two modes:
  gametologs    - Filter paranome TSV to keep only families with ≥1 U and ≥1 V
                  gene; strips all non-UV genes so ksd only computes U-V pairs.
  interspecific - Filter MRBH TSV to remove families where the focal species
                  gene is on a sex chromosome.
"""

import sys
import argparse


def filter_gametologs(tsv_file, u_prefix, v_prefix, output_file):
    kept = 0
    total = 0
    with open(tsv_file) as f, open(output_file, 'w') as out:
        out.write(f.readline())  # header
        for line in f:
            total += 1
            parts = line.strip().split('\t')
            if len(parts) < 2:
                continue
            gf_id = parts[0]
            genes = [g.strip() for g in parts[1].split(',')]
            u_genes = [g for g in genes if g.startswith(u_prefix)]
            v_genes = [g for g in genes if g.startswith(v_prefix)]
            if u_genes and v_genes:
                kept += 1
                out.write(f"{gf_id}\t{', '.join(u_genes + v_genes)}\n")
    print(f"Kept {kept}/{total} families with both U and V genes", file=sys.stderr)


def filter_interspecific(tsv_file, focal, sex_prefixes, output_file):
    kept = 0
    total = 0
    with open(tsv_file) as f, open(output_file, 'w') as out:
        header_line = f.readline()
        out.write(header_line)
        header = header_line.strip().split('\t')
        # Column index: header row has empty first field, then species names
        try:
            header_idx = next(i for i, h in enumerate(header) if focal in h)
        except StopIteration:
            sys.exit(f"ERROR: focal species '{focal}' not found in header: {header}")
        # Data rows have GF ID in column 0; species columns are offset by 1
        col_idx = header_idx + 1
        for line in f:
            total += 1
            parts = line.strip().split('\t')
            if len(parts) <= col_idx:
                continue
            focal_gene = parts[col_idx].strip()
            if any(focal_gene.startswith(p) for p in sex_prefixes):
                continue
            kept += 1
            out.write(line)
    excluded = total - kept
    print(f"Kept {kept}/{total} families (excluded {excluded} with {focal} sex-linked gene)",
          file=sys.stderr)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = parser.add_subparsers(dest='mode', required=True)

    p1 = sub.add_parser('gametologs', help='Filter paranome for U/V gametologs')
    p1.add_argument('tsv', help='Paranome TSV from wgd dmd')
    p1.add_argument('u_prefix', help='U scaffold gene prefix e.g. Papro.U')
    p1.add_argument('v_prefix', help='V scaffold gene prefix e.g. Papro.V')
    p1.add_argument('output', help='Output filtered TSV')

    p2 = sub.add_parser('interspecific', help='Filter MRBH TSV to autosomal genes only')
    p2.add_argument('tsv', help='Interspecific TSV from wgd dmd')
    p2.add_argument('focal', help='Focal species name (substring matching column header)')
    p2.add_argument('sex_prefixes', nargs='+',
                    help='Gene prefixes to exclude e.g. Ledus.U Ledus.V')
    p2.add_argument('-o', '--output', required=True, help='Output filtered TSV')

    args = parser.parse_args()

    if args.mode == 'gametologs':
        filter_gametologs(args.tsv, args.u_prefix, args.v_prefix, args.output)
    else:
        filter_interspecific(args.tsv, args.focal, args.sex_prefixes, args.output)
