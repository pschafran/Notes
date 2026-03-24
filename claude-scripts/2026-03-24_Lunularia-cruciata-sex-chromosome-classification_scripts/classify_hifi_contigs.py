#!/usr/bin/env python3
"""
classify_hifi_contigs.py

Classify contigs from the two Lunularia cruciata HiFi assemblies based on
reciprocal pairwise alignment:

  F2  ERR10480607.bp.p_ctg.fasta  (hifiasm primary, cmLunCruc20, female)
  M2  ERR10480608.bp.p_ctg.fasta  (hifiasm primary, cmLunCruc9,  male)

PAF files (pre-computed with minimap2 -x asm10 --secondary=no):
  ERR10480607.bp.p_ctg.fasta_mapped_to_ERR10480608.bp.p_ctg.fasta.paf  (F2→M2)
  ERR10480608.bp.p_ctg.fasta_mapped_to_ERR10480607.bp.p_ctg.fasta.paf  (M2→F2)
  ERR10480607.bp.p_ctg.fasta_mapped_to_ERR10480607.bp.p_ctg.fasta.paf  (F2 self)
  ERR10480608.bp.p_ctg.fasta_mapped_to_ERR10480608.bp.p_ctg.fasta.paf  (M2 self)

Classification per contig (based on query coverage fraction in the other assembly):
  autosome          ≥ --min-cov covered in the other assembly (mapq ≥ --min-mapq)
  sex_specific      < --low-cov covered in the other assembly at any mapq
  partial           between low-cov and min-cov; ambiguous (may be sex-linked
                    or assembly-fragmented relative to other assembly)
  repeat_ambiguous  only mapq=0 alignments in the other assembly (covered at
                    any mapq ≥ 0 but 0 at mapq ≥ --min-mapq); likely repetitive
  short_contig      below --min-len

Repeat flag (from self-alignment PAF):
  Contig is flagged if it has a non-trivial self-alignment (q_name ≠ t_name, or
  non-overlapping region of same contig) of ≥ --min-aln-len bp at any mapq.

Outputs (to --out dir):
  F2_hifi_classification.tsv   per-contig table for F2
  M2_hifi_classification.tsv   per-contig table for M2
  hifi_classification_summary.tsv  summary counts and Mb per category

Usage:
  python3 classify_hifi_contigs.py [options]

Options:
  --dir DIR          Directory containing assemblies and PAF files [.]
  --out DIR          Output directory [./hifi_classification]
  --min-len INT      Minimum contig length [1000]
  --min-mapq INT     Min mapq to call a contig present in the other assembly [1]
  --min-cov FLOAT    Coverage fraction threshold for 'autosome' [0.5]
  --low-cov FLOAT    Coverage fraction below which contig is 'sex_specific' [0.1]
  --min-aln-len INT  Min alignment length for repeat flagging [500]
"""

import os
import sys
import argparse
from collections import defaultdict


F2_FILE = 'ERR10480607.bp.p_ctg.fasta'
M2_FILE = 'ERR10480608.bp.p_ctg.fasta'


# ── Interval merge ─────────────────────────────────────────────────────────

def merge_covered(intervals):
    """Return total covered bases from a list of (start, end) intervals."""
    if not intervals:
        return 0
    intervals = sorted(intervals)
    total = 0
    cur_s, cur_e = intervals[0]
    for s, e in intervals[1:]:
        if s <= cur_e:
            cur_e = max(cur_e, e)
        else:
            total += cur_e - cur_s
            cur_s, cur_e = s, e
    total += cur_e - cur_s
    return total


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

def parse_cross_paf(paf_path, min_mapq):
    """
    Parse a cross-assembly PAF (query→target).
    Returns two dicts:
      hiq_intervals[qname] = list of (q_start, q_end) for mapq >= min_mapq
      any_intervals[qname] = list of (q_start, q_end) for mapq >= 0
    """
    hiq = defaultdict(list)
    any_ = defaultdict(list)
    if not os.path.exists(paf_path):
        return hiq, any_
    with open(paf_path) as fh:
        for line in fh:
            f = line.rstrip().split('\t')
            if len(f) < 12:
                continue
            mapq    = int(f[11])
            q_name  = f[0]
            q_start = int(f[2])
            q_end   = int(f[3])
            any_[q_name].append((q_start, q_end))
            if mapq >= min_mapq:
                hiq[q_name].append((q_start, q_end))
    return hiq, any_


def parse_self_paf(paf_path, min_aln_len=500):
    """
    Return set of contig names that have non-trivial self-alignments
    (inter-contig or non-overlapping intra-contig hits, length >= min_aln_len).
    Accepts all mapq values (self-alignments are predominantly mapq=0).
    """
    flagged = set()
    if not os.path.exists(paf_path):
        return flagged
    with open(paf_path) as fh:
        for line in fh:
            f = line.rstrip().split('\t')
            if len(f) < 12:
                continue
            q_name  = f[0]
            q_start = int(f[2])
            q_end   = int(f[3])
            t_name  = f[5]
            t_start = int(f[7])
            t_end   = int(f[8])
            if (q_end - q_start) < min_aln_len:
                continue
            if q_name != t_name:
                flagged.add(q_name)
            else:
                overlap = min(q_end, t_end) - max(q_start, t_start)
                if overlap <= 0:
                    flagged.add(q_name)
    return flagged


# ── Classification ─────────────────────────────────────────────────────────

def classify_contig(length, hiq_cov, any_cov, min_len, min_cov, low_cov):
    if length < min_len:
        return 'short_contig'
    hiq_frac = hiq_cov / length
    any_frac = any_cov / length
    if hiq_frac >= min_cov:
        return 'autosome'
    if any_frac >= min_cov:
        return 'repeat_ambiguous'   # covered at mapq=0 but not mapq>=min_mapq
    if any_frac < low_cov:
        return 'sex_specific'
    return 'partial'


# ── Main ───────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--dir',         default='.',  help='PAF/FASTA directory [.]')
    parser.add_argument('--out',         default=None, help='Output directory [./hifi_classification]')
    parser.add_argument('--min-len',     type=int,   default=1000, help='Min contig length [1000]')
    parser.add_argument('--min-mapq',    type=int,   default=1,    help='Min mapq for coverage [1]')
    parser.add_argument('--min-cov',     type=float, default=0.5,  help='Coverage threshold for autosome [0.5]')
    parser.add_argument('--low-cov',     type=float, default=0.1,  help='Max coverage for sex_specific [0.1]')
    parser.add_argument('--min-aln-len', type=int,   default=500,  help='Min aln length for repeat flag [500]')
    args = parser.parse_args()

    d = args.dir
    out_dir = args.out or os.path.join(d, 'hifi_classification')
    os.makedirs(out_dir, exist_ok=True)

    print(f"Parameters: min_len={args.min_len}, min_mapq={args.min_mapq}, "
          f"min_cov={args.min_cov}, low_cov={args.low_cov}, min_aln_len={args.min_aln_len}")

    # ── Lengths ────────────────────────────────────────────────────────────
    print("\nParsing contig lengths...")
    f2_len = parse_fasta_lengths(os.path.join(d, F2_FILE))
    m2_len = parse_fasta_lengths(os.path.join(d, M2_FILE))
    print(f"  F2: {len(f2_len):,} contigs, {sum(f2_len.values())/1e6:.1f} Mb")
    print(f"  M2: {len(m2_len):,} contigs, {sum(m2_len.values())/1e6:.1f} Mb")

    # ── Self-alignment repeat flags ────────────────────────────────────────
    print("\nParsing self-alignments for repeat flags...")
    f2_rep = parse_self_paf(os.path.join(d, f"{F2_FILE}_mapped_to_{F2_FILE}.paf"), args.min_aln_len)
    m2_rep = parse_self_paf(os.path.join(d, f"{M2_FILE}_mapped_to_{M2_FILE}.paf"), args.min_aln_len)
    print(f"  F2: {len(f2_rep):,} repeat-flagged contigs")
    print(f"  M2: {len(m2_rep):,} repeat-flagged contigs")

    # ── Cross-assembly coverage ────────────────────────────────────────────
    print("\nParsing cross-assembly alignments...")
    # F2 contigs queried against M2
    f2_hiq, f2_any = parse_cross_paf(
        os.path.join(d, f"{F2_FILE}_mapped_to_{M2_FILE}.paf"), args.min_mapq)
    # M2 contigs queried against F2
    m2_hiq, m2_any = parse_cross_paf(
        os.path.join(d, f"{M2_FILE}_mapped_to_{F2_FILE}.paf"), args.min_mapq)

    # ── Classify and write ─────────────────────────────────────────────────
    CLASSES = ['autosome', 'sex_specific', 'partial', 'repeat_ambiguous', 'short_contig']
    summary = {}  # {asm: {cls: [n, bp]}}

    for asm_id, lengths, hiq_cov, any_cov, rep_flags, sex_label in [
        ('F2', f2_len, f2_hiq, f2_any, f2_rep, 'U_candidate'),
        ('M2', m2_len, m2_hiq, m2_any, m2_rep, 'V_candidate'),
    ]:
        summary[asm_id] = defaultdict(lambda: [0, 0])
        out_tsv = os.path.join(out_dir, f"{asm_id}_hifi_classification.tsv")

        with open(out_tsv, 'w') as out:
            out.write('\t'.join([
                'contig', 'length', 'classification', 'repeat_flagged',
                'hiq_cov_bases', 'hiq_cov_frac',
                'any_cov_bases', 'any_cov_frac',
            ]) + '\n')

            for contig, length in sorted(lengths.items(), key=lambda x: -x[1]):
                hiq_bases = merge_covered(hiq_cov.get(contig, []))
                any_bases = merge_covered(any_cov.get(contig, []))
                cls = classify_contig(length, hiq_bases, any_bases,
                                      args.min_len, args.min_cov, args.low_cov)
                is_rep = 'yes' if contig in rep_flags else 'no'
                summary[asm_id][cls][0] += 1
                summary[asm_id][cls][1] += length
                out.write('\t'.join([
                    contig, str(length), cls, is_rep,
                    str(hiq_bases), f"{hiq_bases/length:.4f}",
                    str(any_bases), f"{any_bases/length:.4f}",
                ]) + '\n')

        print(f"\n  {asm_id} written to {out_tsv}")
        for cls in CLASSES:
            n, bp = summary[asm_id][cls]
            if n > 0:
                print(f"    {cls:<20} {n:>6,} contigs  {bp/1e6:>8.2f} Mb")

    # ── Summary TSV ────────────────────────────────────────────────────────
    summary_tsv = os.path.join(out_dir, 'hifi_classification_summary.tsv')
    with open(summary_tsv, 'w') as out:
        out.write('classification\tn_F2\tMb_F2\tn_M2\tMb_M2\n')
        for cls in CLASSES:
            n_f2, bp_f2 = summary['F2'][cls]
            n_m2, bp_m2 = summary['M2'][cls]
            out.write(f"{cls}\t{n_f2}\t{bp_f2/1e6:.3f}\t{n_m2}\t{bp_m2/1e6:.3f}\n")

    print(f"\nSummary written to {summary_tsv}")

    # ── Sex-specific FASTA lists ───────────────────────────────────────────
    # Write contig name lists for sex_specific contigs (for downstream use)
    for asm_id, lengths, hiq_cov, any_cov, sex_label in [
        ('F2', f2_len, f2_hiq, f2_any, 'U_candidate'),
        ('M2', m2_len, m2_hiq, m2_any, 'V_candidate'),
    ]:
        list_path = os.path.join(out_dir, f"{asm_id}_{sex_label}_contigs.txt")
        with open(list_path, 'w') as out:
            for contig, length in sorted(lengths.items(), key=lambda x: -x[1]):
                hiq_bases = merge_covered(hiq_cov.get(contig, []))
                any_bases = merge_covered(any_cov.get(contig, []))
                cls = classify_contig(length, hiq_bases, any_bases,
                                      args.min_len, args.min_cov, args.low_cov)
                if cls == 'sex_specific':
                    out.write(contig + '\n')
        n = sum(1 for _ in open(list_path))
        print(f"  {asm_id} {sex_label} list: {n} contigs → {list_path}")


if __name__ == '__main__':
    main()
