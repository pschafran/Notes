#!/usr/bin/env python3
"""
classify_contigs.py

Classify contigs from four Lunularia cruciata assemblies based on pairwise
whole-genome alignment patterns.

Assemblies:
  F1  Lunularia_cruciata_female.faa       Linde et al. 2021 female (fragmented)
  F2  ERR10480607.bp.p_ctg.fasta          HiFi assembly, cmLunCruc20 (female)
  M1  Lunularia_cruciata_male.faa         Linde et al. 2021 male (fragmented)
  M2  ERR10480608.bp.p_ctg.fasta          HiFi assembly, cmLunCruc9 (male)

PAF files were generated with:
  minimap2 -t 24 -x asm10 $target $query > {query}_mapped_to_{target}.paf

Classifications:
  autosome              present in all 4 assemblies
  U_sex_chromosome      present in both females (F1+F2), absent from both males
  V_sex_chromosome      present in both males (M1+M2), absent from both females
  present_3_missing_F1  present in F2+M1+M2 only
  present_3_missing_F2  present in F1+M1+M2 only
  present_3_missing_M1  present in F1+F2+M2 only
  present_3_missing_M2  present in F1+F2+M1 only
  cross_sex_linde_pair  present in F1+M1 only (same study, cross-sex)
  cross_sex_hifi_pair   present in F2+M2 only (same technology, cross-sex)
  cross_sex_mixed       present in F1+M2 or F2+M1 only
  assembly_specific     present in exactly one assembly
  short_contig          below --min-len threshold (not classified)

Repeat flag (applied on top of classification):
  A contig is repeat-flagged if, in the self-alignment PAF, it has alignments
  to other contigs in the same assembly (inter-contig repeats), or to
  non-overlapping regions of itself (internal tandem/inverted repeats).
  Only alignments >= 500 bp are considered (noise filter).
  minimap2 self-alignment hits are predominantly mapq=0, so --repeat-mapq
  defaults to 0 (accept all). Classifications of repeat-flagged contigs
  are unreliable.

Usage:
  python3 classify_contigs.py [options]

Options:
  --dir DIR          Directory containing assemblies and PAF files [.]
  --out DIR          Output directory [./contig_classification]
  --min-len INT      Minimum contig length to classify [1000]
  --min-mapq INT     Minimum mapping quality for cross-assembly hits [1]
  --min-cov FLOAT    Minimum query coverage fraction to call contig present [0.5]
  --repeat-mapq INT  Minimum mapq to call a self-alignment repeat hit [0]
"""

import os
import sys
import argparse
from collections import defaultdict


# ── Assembly definitions ───────────────────────────────────────────────────

ASSEMBLIES = {
    'F1': 'Lunularia_cruciata_female.faa',
    'F2': 'ERR10480607.bp.p_ctg.fasta',
    'M1': 'Lunularia_cruciata_male.faa',
    'M2': 'ERR10480608.bp.p_ctg.fasta',
}
ASM_ORDER = ['F1', 'F2', 'M1', 'M2']
FEMALE = {'F1', 'F2'}
MALE   = {'M1', 'M2'}


# ── Classification logic ───────────────────────────────────────────────────

def classify_pattern(pattern):
    """
    Classify a contig based on its presence/absence pattern across assemblies.
    pattern: dict {asm_id: bool} e.g. {'F1': True, 'F2': False, ...}
    Returns classification string.
    """
    f1, f2, m1, m2 = pattern['F1'], pattern['F2'], pattern['M1'], pattern['M2']
    n_present = sum([f1, f2, m1, m2])

    if n_present == 4:
        return 'autosome'

    if n_present == 2:
        if f1 and f2:   return 'U_sex_chromosome'
        if m1 and m2:   return 'V_sex_chromosome'
        if f1 and m1:   return 'cross_sex_linde_pair'
        if f2 and m2:   return 'cross_sex_hifi_pair'
        return 'cross_sex_mixed'            # F1+M2 or F2+M1

    if n_present == 3:
        if not f1:  return 'present_3_missing_F1'
        if not f2:  return 'present_3_missing_F2'
        if not m1:  return 'present_3_missing_M1'
        if not m2:  return 'present_3_missing_M2'

    if n_present == 1:
        return 'assembly_specific'

    # n_present == 0 should never happen; covers edge cases
    return 'unclassified'


# ── Interval utilities ─────────────────────────────────────────────────────

def merge_covered_length(intervals):
    """Merge overlapping [start, end) intervals and return total covered bases."""
    if not intervals:
        return 0
    intervals = sorted(intervals)
    merged_end = intervals[0][1]
    total = 0
    for start, end in intervals:
        if start >= merged_end:
            total += merged_end - intervals[0][0] if total == 0 else 0
            # restart
            intervals[0] = (start, end)
            merged_end = end
        else:
            merged_end = max(merged_end, end)
    # simpler approach:
    return _merge(intervals)

def _merge(intervals):
    if not intervals:
        return 0
    intervals = sorted(intervals)
    total = 0
    cur_start, cur_end = intervals[0]
    for start, end in intervals[1:]:
        if start <= cur_end:
            cur_end = max(cur_end, end)
        else:
            total += cur_end - cur_start
            cur_start, cur_end = start, end
    total += cur_end - cur_start
    return total


# ── FASTA parsing ──────────────────────────────────────────────────────────

def parse_fasta_lengths(fasta_path):
    """Return {contig_name: length} from a FASTA file."""
    lengths = {}
    name, length = None, 0
    with open(fasta_path) as fh:
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

def parse_paf_query_coverage(paf_path, min_mapq):
    """
    Parse a PAF file and return per-query-contig merged query coverage.
    Returns {query_name: covered_bases (int)}
    """
    intervals = defaultdict(list)    # query_name -> [(q_start, q_end)]
    if not os.path.exists(paf_path):
        return {}
    with open(paf_path) as fh:
        for line in fh:
            f = line.rstrip().split('\t')
            if len(f) < 12:
                continue
            mapq = int(f[11])
            if mapq < min_mapq:
                continue
            q_name  = f[0]
            q_start = int(f[2])
            q_end   = int(f[3])
            intervals[q_name].append((q_start, q_end))
    return {name: _merge(ivs) for name, ivs in intervals.items()}


def parse_self_paf_repeat_flags(paf_path, repeat_mapq, min_aln_len=500):
    """
    Parse the self-alignment PAF to identify repeat-problematic contigs.
    A contig is flagged if it has any non-trivial alignment (mapq >= repeat_mapq,
    alignment length >= min_aln_len), defined as:
      (a) query_name != target_name  (aligns to a different contig)
      (b) query_name == target_name but query and target intervals do not overlap
          (internal repeat at non-overlapping position)

    Note: minimap2 self-alignments produce predominantly mapq=0 hits, so the
    default repeat_mapq=0 accepts all alignments; min_aln_len filters noise.
    The trivial self-hit (q_name == t_name, overlapping coords) is excluded.
    Returns set of flagged contig names.
    """
    flagged = set()
    if not os.path.exists(paf_path):
        return flagged
    with open(paf_path) as fh:
        for line in fh:
            f = line.rstrip().split('\t')
            if len(f) < 12:
                continue
            mapq = int(f[11])
            if mapq < repeat_mapq:
                continue
            q_name  = f[0]
            q_len   = int(f[1])
            q_start = int(f[2])
            q_end   = int(f[3])
            t_name  = f[5]
            t_start = int(f[7])
            t_end   = int(f[8])

            # Skip alignments shorter than min_aln_len (noise / trivial overlaps)
            aln_len = q_end - q_start
            if aln_len < min_aln_len:
                continue

            if q_name != t_name:
                # (a) aligns to a different contig in the same assembly
                flagged.add(q_name)
            else:
                # (b) same contig: check for non-overlapping (internal repeat)
                overlap = min(q_end, t_end) - max(q_start, t_start)
                if overlap <= 0:
                    flagged.add(q_name)
    return flagged


# ── Main ───────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--dir',          default='.',   help='Directory with assemblies and PAFs')
    parser.add_argument('--out',          default=None,  help='Output directory [--dir/contig_classification]')
    parser.add_argument('--min-len',      type=int,   default=1000,  help='Min contig length [1000]')
    parser.add_argument('--min-mapq',     type=int,   default=1,     help='Min mapq for cross-assembly hits [1]')
    parser.add_argument('--min-cov',      type=float, default=0.5,   help='Min query coverage fraction to call present [0.5]')
    parser.add_argument('--repeat-mapq',  type=int,   default=0,     help='Min mapq for repeat flagging [0]')
    args = parser.parse_args()

    work_dir = args.dir
    out_dir  = args.out or os.path.join(work_dir, 'contig_classification')
    os.makedirs(out_dir, exist_ok=True)

    print(f"Working directory: {work_dir}")
    print(f"Output directory:  {out_dir}")
    print(f"Parameters: min_len={args.min_len}, min_mapq={args.min_mapq}, "
          f"min_cov={args.min_cov}, repeat_mapq={args.repeat_mapq}\n")

    # ── Step 1: Parse contig lengths ───────────────────────────────────────
    print("Parsing contig lengths...")
    contig_lengths = {}    # {asm_id: {contig: length}}
    for asm_id, asm_file in ASSEMBLIES.items():
        path = os.path.join(work_dir, asm_file)
        contig_lengths[asm_id] = parse_fasta_lengths(path)
        n = len(contig_lengths[asm_id])
        total_bp = sum(contig_lengths[asm_id].values())
        print(f"  {asm_id} ({asm_file}): {n:,} contigs, {total_bp/1e6:.1f} Mb")

    # ── Step 2: Parse self-alignments for repeat flags ─────────────────────
    print("\nParsing self-alignments for repeat flags...")
    repeat_flags = {}    # {asm_id: set of flagged contig names}
    for asm_id, asm_file in ASSEMBLIES.items():
        paf = os.path.join(work_dir, f"{asm_file}_mapped_to_{asm_file}.paf")
        repeat_flags[asm_id] = parse_self_paf_repeat_flags(paf, args.repeat_mapq)
        print(f"  {asm_id}: {len(repeat_flags[asm_id]):,} repeat-flagged contigs")

    # ── Step 3: Parse cross-assembly coverage ─────────────────────────────
    # query_cov[query_asm][target_asm][contig] = covered_bases
    print("\nParsing cross-assembly alignments...")
    query_cov = {q: {} for q in ASM_ORDER}
    for q_id, q_file in ASSEMBLIES.items():
        for t_id, t_file in ASSEMBLIES.items():
            paf = os.path.join(work_dir, f"{q_file}_mapped_to_{t_file}.paf")
            cov = parse_paf_query_coverage(paf, args.min_mapq)
            query_cov[q_id][t_id] = cov
            n_hits = len(cov)
            print(f"  {q_id} -> {t_id}: {n_hits:,} contigs with alignments")

    # ── Step 4: Classify each contig ──────────────────────────────────────
    print("\nClassifying contigs...")

    # Summary counters: {asm_id: {classification: [count, total_bp]}}
    summary = {a: defaultdict(lambda: [0, 0]) for a in ASM_ORDER}

    for asm_id, asm_file in ASSEMBLIES.items():
        out_tsv = os.path.join(out_dir, f"{asm_id}_contig_classification.tsv")
        lengths  = contig_lengths[asm_id]
        flagged  = repeat_flags[asm_id]

        with open(out_tsv, 'w') as out:
            header = ['contig', 'length',
                      'classification', 'repeat_flagged', 'pattern',
                      'cov_F1', 'cov_F2', 'cov_M1', 'cov_M2']
            out.write('\t'.join(header) + '\n')

            for contig, length in sorted(lengths.items(), key=lambda x: -x[1]):

                # Short contig
                if length < args.min_len:
                    classification = 'short_contig'
                    summary[asm_id][classification][0] += 1
                    summary[asm_id][classification][1] += length
                    out.write(f"{contig}\t{length}\t{classification}\t"
                              f"{'yes' if contig in flagged else 'no'}\tNA"
                              f"\tNA\tNA\tNA\tNA\n")
                    continue

                # Compute query coverage fraction in each assembly
                cov_frac = {}
                for t_id in ASM_ORDER:
                    covered = query_cov[asm_id][t_id].get(contig, 0)
                    cov_frac[t_id] = covered / length

                # Presence/absence pattern
                pattern = {t: cov_frac[t] >= args.min_cov for t in ASM_ORDER}
                # Always present in self
                pattern[asm_id] = True

                classification = classify_pattern(pattern)
                is_repeat = contig in flagged
                pattern_str = ''.join('1' if pattern[a] else '0' for a in ASM_ORDER)

                summary[asm_id][classification][0] += 1
                summary[asm_id][classification][1] += length

                out.write(
                    f"{contig}\t{length}\t{classification}\t"
                    f"{'yes' if is_repeat else 'no'}\t{pattern_str}\t"
                    f"{cov_frac['F1']:.3f}\t{cov_frac['F2']:.3f}\t"
                    f"{cov_frac['M1']:.3f}\t{cov_frac['M2']:.3f}\n"
                )

        n_classified = sum(v[0] for k, v in summary[asm_id].items() if k != 'short_contig')
        print(f"  {asm_id}: written to {out_tsv} ({n_classified:,} contigs classified)")

    # ── Step 5: Write summary ──────────────────────────────────────────────
    all_classes = [
        'autosome',
        'U_sex_chromosome', 'V_sex_chromosome',
        'present_3_missing_F1', 'present_3_missing_F2',
        'present_3_missing_M1', 'present_3_missing_M2',
        'cross_sex_linde_pair', 'cross_sex_hifi_pair', 'cross_sex_mixed',
        'assembly_specific',
        'short_contig', 'unclassified',
    ]

    summary_tsv = os.path.join(out_dir, 'classification_summary.tsv')
    with open(summary_tsv, 'w') as out:
        header = ['classification']
        for a in ASM_ORDER:
            header += [f'n_{a}', f'Mb_{a}']
        out.write('\t'.join(header) + '\n')

        for cls in all_classes:
            row = [cls]
            for a in ASM_ORDER:
                n, bp = summary[a][cls]
                row += [str(n), f"{bp/1e6:.2f}"]
            out.write('\t'.join(row) + '\n')

    print(f"\nSummary written to {summary_tsv}")

    # Print summary to stdout
    print(f"\n{'Classification':<30} {'F1':>8} {'F2':>8} {'M1':>8} {'M2':>8}  (contig counts)")
    print('-' * 70)
    for cls in all_classes:
        counts = [summary[a][cls][0] for a in ASM_ORDER]
        if any(c > 0 for c in counts):
            print(f"{cls:<30} {counts[0]:>8,} {counts[1]:>8,} {counts[2]:>8,} {counts[3]:>8,}")

    print(f"\n{'Classification':<30} {'F1':>8} {'F2':>8} {'M1':>8} {'M2':>8}  (Mb)")
    print('-' * 70)
    for cls in all_classes:
        mbs = [summary[a][cls][1]/1e6 for a in ASM_ORDER]
        if any(m > 0 for m in mbs):
            print(f"{cls:<30} {mbs[0]:>8.1f} {mbs[1]:>8.1f} {mbs[2]:>8.1f} {mbs[3]:>8.1f}")


if __name__ == '__main__':
    main()
