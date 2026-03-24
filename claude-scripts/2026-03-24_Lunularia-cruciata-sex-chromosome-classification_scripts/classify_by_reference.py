#!/usr/bin/env python3
"""
classify_by_reference.py

Classify contigs from four Lunularia cruciata assemblies by mapping to the
chromosome-scale reference (GCA_948567375.1_cmLunCruc20.1_genome.fa), with
repeat annotation cross-referencing to distinguish genuine alignments from
repeat-driven ones.

Each contig is assigned to the reference chromosome category with the most
aligned query bases (at mapq >= --min-mapq). For each alignment, the overlap
of the *target* (reference) interval with EDTA repeat annotations is computed.
Per-contig, we report the fraction of total aligned bases that land on
repeat-annotated reference positions.

Reference categories (from CLAUDE.md):
  autosome      OX419755.1–OX419762.1  (chr1–8)
  U_chromosome  OX419763.1             (chr9, 1.76 Mb)
                CAOYZS010000005.1      (chr9_unloc1, 1.25 Mb)
                CAOYZS010000006.1      (chr9_unloc2, 1.15 Mb)
                CAOYZS010000008.1      (chr9_unloc4, 0.14 Mb)
  organellar    OX419764.1 (mitochondrion), OX419765.1 (plastid)
  unplaced      all other CAOYZS010000*.1 scaffolds
  no_hit        no alignment at mapq >= --min-mapq

Note: CAOYZS010000007.1 (chr9_unloc3) has autosome-like coverage; treated
as 'unplaced' (see CLAUDE.md).

Repeat overlap columns (based on target/reference coordinates):
  aln_bases        total merged query bases aligned (mapq >= --min-mapq)
  aln_frac         aln_bases / contig_length
  repeat_aln_bases aligned bases whose reference position overlaps a repeat
  repeat_aln_frac  repeat_aln_bases / aln_bases  (0 if no alignment)
  repeat_class     non_repeat (<20%), mixed (20-80%), repeat_driven (>80%)

PAF files expected in --dir:
  {asm_file}_mapped_to_GCA_948567375.1_cmLunCruc20.1_genome.fa.paf

Assemblies:
  F1  Lunularia_cruciata_female.faa      Linde et al. 2021 female (fragmented)
  F2  ERR10480607.bp.p_ctg.fasta        HiFi, cmLunCruc20 (female)
  M1  Lunularia_cruciata_male.faa       Linde et al. 2021 male (fragmented)
  M2  ERR10480608.bp.p_ctg.fasta        HiFi, cmLunCruc9 (male)

Outputs (in --out):
  {ASM}_ref_classification.tsv   per-contig table
  ref_classification_summary.tsv  count and Mb per category per assembly

Usage:
  python3 classify_by_reference.py [options]

Options:
  --dir DIR        Directory with assemblies and PAF files [.]
  --repeat-gff FILE  EDTA repeat annotation GFF (reference coordinates)
  --out DIR        Output directory [./ref_classification]
  --min-mapq INT   Minimum mapping quality [1]
  --min-len INT    Minimum contig length to classify [1000]
"""

import os
import sys
import argparse
import bisect
from collections import defaultdict

# ── Reference sequence categories ─────────────────────────────────────────

AUTOSOMES = {
    'OX419755.1': 'chr1', 'OX419756.1': 'chr2', 'OX419757.1': 'chr3',
    'OX419758.1': 'chr4', 'OX419759.1': 'chr5', 'OX419760.1': 'chr6',
    'OX419761.1': 'chr7', 'OX419762.1': 'chr8',
}
U_SEQS = {
    'OX419763.1':        'chr9',
    'CAOYZS010000005.1': 'chr9_unloc1',
    'CAOYZS010000006.1': 'chr9_unloc2',
    'CAOYZS010000008.1': 'chr9_unloc4',
}
ORGANELLAR = {
    'OX419764.1': 'mitochondrion',
    'OX419765.1': 'plastid',
}

def ref_category(seq_id):
    if seq_id in AUTOSOMES:  return 'autosome'
    if seq_id in U_SEQS:     return 'U_chromosome'
    if seq_id in ORGANELLAR: return 'organellar'
    return 'unplaced'


# ── Assembly definitions ───────────────────────────────────────────────────

ASSEMBLIES = {
    'F1': 'Lunularia_cruciata_female.faa',
    'F2': 'ERR10480607.bp.p_ctg.fasta',
    'M1': 'Lunularia_cruciata_male.faa',
    'M2': 'ERR10480608.bp.p_ctg.fasta',
}
REF_SUFFIX = '_mapped_to_GCA_948567375.1_cmLunCruc20.1_genome.fa.paf'


# ── Interval utilities ─────────────────────────────────────────────────────

def merge_intervals(intervals):
    """Merge sorted list of (start, end) intervals, return merged list."""
    if not intervals:
        return []
    merged = [list(intervals[0])]
    for s, e in intervals[1:]:
        if s <= merged[-1][1]:
            merged[-1][1] = max(merged[-1][1], e)
        else:
            merged.append([s, e])
    return merged


def covered_length(intervals):
    """Total bases covered by a list of (start, end) intervals (unsorted ok)."""
    if not intervals:
        return 0
    ivs = sorted(intervals)
    total = 0
    cs, ce = ivs[0]
    for s, e in ivs[1:]:
        if s <= ce:
            ce = max(ce, e)
        else:
            total += ce - cs
            cs, ce = s, e
    total += ce - cs
    return total


def overlap_with_sorted_merged(query_s, query_e, sorted_merged):
    """
    Compute bases of [query_s, query_e) that overlap any interval in
    sorted_merged (a list of [start, end] already sorted and merged).
    Uses binary search for efficiency.
    """
    if not sorted_merged:
        return 0
    # Find first interval whose end > query_s
    starts = [iv[0] for iv in sorted_merged]
    idx = bisect.bisect_right(starts, query_s) - 1
    idx = max(idx, 0)
    total = 0
    for i in range(idx, len(sorted_merged)):
        rs, re = sorted_merged[i]
        if rs >= query_e:
            break
        ov = min(query_e, re) - max(query_s, rs)
        if ov > 0:
            total += ov
    return total


# ── GFF repeat loader ──────────────────────────────────────────────────────

def load_repeat_intervals(gff_path):
    """
    Parse EDTA GFF3 and return {chrom: sorted_merged_intervals}.
    GFF is 1-based closed; convert to 0-based half-open.
    """
    print(f"Loading repeat annotations from {gff_path}...")
    raw = defaultdict(list)
    n = 0
    with open(gff_path) as fh:
        for line in fh:
            if line.startswith('#'):
                continue
            f = line.rstrip().split('\t')
            if len(f) < 5:
                continue
            chrom = f[0]
            start = int(f[3]) - 1   # GFF 1-based → 0-based
            end   = int(f[4])        # GFF closed → half-open
            raw[chrom].append((start, end))
            n += 1
    print(f"  {n:,} repeat features on {len(raw)} sequences")
    # Sort and merge per chromosome
    merged = {}
    for chrom, ivs in raw.items():
        merged[chrom] = merge_intervals(sorted(ivs))
    return merged


# ── FASTA parsing ──────────────────────────────────────────────────────────

def parse_fasta_lengths(path):
    lengths = {}
    name, length = None, 0
    with open(path) as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith('>'):
                if name is not None:
                    lengths[name] = length
                name = line[1:].split()[0]
                length = 0
            else:
                length += len(line)
    if name is not None:
        lengths[name] = length
    return lengths


# ── PAF parsing ────────────────────────────────────────────────────────────

def parse_ref_paf(paf_path, min_mapq, repeat_ivs):
    """
    For each query contig, compute:
      - per-category merged query coverage (for classification)
      - total aligned query bases
      - repeat-overlapping aligned query bases (based on target coords)

    Returns {q_name: dict} with keys:
      best_cat, best_seq, best_seq_cov, cat_cov,
      aln_bases, repeat_aln_bases
    """
    if not os.path.exists(paf_path):
        print(f"  WARNING: {paf_path} not found")
        return {}

    # {q_name: {t_name: [(q_start, q_end)]}}  for coverage
    q_intervals  = defaultdict(lambda: defaultdict(list))
    # {q_name: [(q_start, q_end)]}  all alignments (for total aln_bases)
    all_q_ivs    = defaultdict(list)
    # {q_name: repeat overlap bases}
    rep_overlap  = defaultdict(int)

    with open(paf_path) as fh:
        for line in fh:
            f = line.rstrip().split('\t')
            if len(f) < 12:
                continue
            if int(f[11]) < min_mapq:
                continue
            q_name  = f[0]
            q_start = int(f[2])
            q_end   = int(f[3])
            t_name  = f[5]
            t_start = int(f[7])
            t_end   = int(f[8])

            q_intervals[q_name][t_name].append((q_start, q_end))
            all_q_ivs[q_name].append((q_start, q_end))

            # Repeat overlap on target side, scaled to query length of aln
            aln_len = q_end - q_start
            if aln_len > 0 and t_name in repeat_ivs:
                t_rep = overlap_with_sorted_merged(t_start, t_end, repeat_ivs[t_name])
                # Scale target repeat overlap to query alignment length
                t_aln_len = t_end - t_start
                if t_aln_len > 0:
                    rep_overlap[q_name] += int(aln_len * t_rep / t_aln_len)

    result = {}
    for q_name, t_ivs in q_intervals.items():
        cat_cov = defaultdict(int)
        seq_cov = {}
        for t_name, ivs in t_ivs.items():
            cov = covered_length(ivs)
            seq_cov[t_name] = cov
            cat_cov[ref_category(t_name)] += cov

        best_cat = max(cat_cov, key=lambda c: cat_cov[c])
        best_seq = max(seq_cov, key=lambda s: seq_cov[s])

        aln_bases = covered_length(all_q_ivs[q_name])

        result[q_name] = {
            'best_cat':       best_cat,
            'best_seq':       best_seq,
            'best_seq_cov':   seq_cov[best_seq],
            'cat_cov':        dict(cat_cov),
            'aln_bases':      aln_bases,
            'repeat_aln_bases': rep_overlap.get(q_name, 0),
        }
    return result


def repeat_class(rep_frac):
    if rep_frac >= 0.8:  return 'repeat_driven'
    if rep_frac >= 0.2:  return 'mixed'
    return 'non_repeat'


# ── Main ───────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--dir',        default='.',   help='Working directory [.]')
    parser.add_argument('--repeat-gff', default=None,  help='EDTA repeat GFF file')
    parser.add_argument('--out',        default=None)
    parser.add_argument('--min-mapq',   type=int, default=1)
    parser.add_argument('--min-len',    type=int, default=1000)
    args = parser.parse_args()

    d = args.dir
    out_dir = args.out or os.path.join(d, 'ref_classification')
    os.makedirs(out_dir, exist_ok=True)

    print(f"Parameters: min_mapq={args.min_mapq}, min_len={args.min_len}")

    # Load repeat annotations
    repeat_ivs = {}
    if args.repeat_gff:
        repeat_ivs = load_repeat_intervals(args.repeat_gff)
    else:
        print("No --repeat-gff provided; repeat overlap columns will be 0.")

    CATS = ['autosome', 'U_chromosome', 'organellar', 'unplaced', 'no_hit', 'short_contig']
    summary = {}

    for asm_id, asm_file in ASSEMBLIES.items():
        print(f"\n── {asm_id} ({asm_file}) ──")

        lengths = parse_fasta_lengths(os.path.join(d, asm_file))
        print(f"  {len(lengths):,} contigs, {sum(lengths.values())/1e6:.1f} Mb")

        paf_path = os.path.join(d, asm_file + REF_SUFFIX)
        ref_hits = parse_ref_paf(paf_path, args.min_mapq, repeat_ivs)
        print(f"  {len(ref_hits):,} contigs with ref alignment (mapq≥{args.min_mapq})")

        summary[asm_id] = defaultdict(lambda: [0, 0])
        out_tsv = os.path.join(out_dir, f"{asm_id}_ref_classification.tsv")

        with open(out_tsv, 'w') as out:
            out.write('\t'.join([
                'contig', 'length', 'ref_category', 'best_ref_seq',
                'aln_bases', 'aln_frac',
                'repeat_aln_bases', 'repeat_aln_frac', 'repeat_class',
                'autosome_cov', 'U_cov', 'organellar_cov', 'unplaced_cov',
            ]) + '\n')

            for contig, length in sorted(lengths.items(), key=lambda x: -x[1]):
                if length < args.min_len:
                    cat = 'short_contig'
                    row_extra = ['NA', 'NA', '0', '0', '0.0000', 'NA', '0', '0', '0', '0']
                elif contig in ref_hits:
                    h = ref_hits[contig]
                    cat       = h['best_cat']
                    aln       = h['aln_bases']
                    rep       = min(h['repeat_aln_bases'], aln)  # cap at aln_bases
                    rep_frac  = rep / aln if aln > 0 else 0.0
                    row_extra = [
                        str(h['best_seq']),
                        str(aln), f"{aln/length:.4f}",
                        str(rep), f"{rep_frac:.4f}", repeat_class(rep_frac),
                        str(h['cat_cov'].get('autosome',    0)),
                        str(h['cat_cov'].get('U_chromosome', 0)),
                        str(h['cat_cov'].get('organellar',  0)),
                        str(h['cat_cov'].get('unplaced',    0)),
                    ]
                else:
                    cat = 'no_hit'
                    row_extra = ['NA', '0', '0.0000', '0', '0.0000', 'NA', '0', '0', '0', '0']

                summary[asm_id][cat][0] += 1
                summary[asm_id][cat][1] += length
                out.write('\t'.join([contig, str(length), cat] + row_extra) + '\n')

        print(f"  Written: {out_tsv}")
        for cat in CATS:
            n, bp = summary[asm_id][cat]
            if n > 0:
                print(f"    {cat:<15} {n:>7,} contigs  {bp/1e6:>8.2f} Mb")

    # ── Per-category repeat breakdown for key categories ──────────────────
    print("\n── Repeat classification breakdown (U_chromosome & no_hit contigs) ──")
    for asm_id in ('F1', 'F2', 'M1', 'M2'):
        tsv = os.path.join(out_dir, f"{asm_id}_ref_classification.tsv")
        counts = defaultdict(lambda: defaultdict(lambda: [0, 0]))
        with open(tsv) as fh:
            next(fh)
            for line in fh:
                f = line.rstrip().split('\t')
                cat, rclass, length = f[2], f[8], int(f[1])
                if cat in ('U_chromosome', 'no_hit') and rclass != 'NA':
                    counts[cat][rclass][0] += 1
                    counts[cat][rclass][1] += length
        for cat in ('U_chromosome', 'no_hit'):
            total_n = sum(v[0] for v in counts[cat].values())
            if total_n == 0:
                continue
            print(f"  {asm_id} {cat}:")
            for rc in ('non_repeat', 'mixed', 'repeat_driven'):
                n, bp = counts[cat][rc]
                print(f"    {rc:<15} {n:>5,} contigs  {bp/1e6:>6.2f} Mb")

    # ── Summary TSV ────────────────────────────────────────────────────────
    summary_tsv = os.path.join(out_dir, 'ref_classification_summary.tsv')
    with open(summary_tsv, 'w') as out:
        header = ['category']
        for a in ('F1', 'F2', 'M1', 'M2'):
            header += [f'n_{a}', f'Mb_{a}']
        out.write('\t'.join(header) + '\n')
        for cat in CATS:
            row = [cat]
            for a in ('F1', 'F2', 'M1', 'M2'):
                n, bp = summary[a][cat]
                row += [str(n), f"{bp/1e6:.3f}"]
            out.write('\t'.join(row) + '\n')

    print(f"\nSummary: {summary_tsv}")


if __name__ == '__main__':
    main()
