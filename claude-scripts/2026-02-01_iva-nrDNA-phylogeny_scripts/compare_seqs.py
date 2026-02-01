#!/usr/bin/env python3
"""
compare_seqs.py - Compare FASTA sequences for differences

Usage:
    python compare_seqs.py file1.fasta file2.fasta          # Compare two sequences
    python compare_seqs.py file1.fasta file2.fasta ...      # Compare multiple sequences

Created during Iva nrDNA phylogeny analysis session, 2026-02-01
"""
import sys
from itertools import combinations

def read_fasta(fname):
    with open(fname) as f:
        lines = f.readlines()
        return ''.join(l.strip() for l in lines if not l.startswith('>'))

def compare_two(f1, f2):
    s1 = read_fasta(f1)
    s2 = read_fasta(f2)
    if s1 == s2:
        print("Sequences are IDENTICAL")
        return

    diffs = [(i, s1[i], s2[i]) for i in range(min(len(s1), len(s2))) if s1[i] != s2[i]]
    print(f"Number of SNP differences: {len(diffs)}")
    for pos, b1, b2 in diffs[:20]:
        print(f"  Position {pos+1}: {b1} -> {b2}")
    if len(diffs) > 20:
        print(f"  ... and {len(diffs)-20} more")

def compare_multiple(files):
    seqs = {}
    for f in files:
        name = f.split('/')[-1].replace('.path_sequence.fasta', '').split('.')[-2] + '.' + f.split('/')[-1].replace('.path_sequence.fasta', '').split('.')[-1]
        seqs[name] = read_fasta(f)

    print("Pairwise comparisons:")
    for (n1, s1), (n2, s2) in combinations(seqs.items(), 2):
        if s1 == s2:
            print(f"  {n1} vs {n2}: IDENTICAL")
        else:
            diffs = [(i, s1[i], s2[i]) for i in range(len(s1)) if s1[i] != s2[i]]
            print(f"  {n1} vs {n2}: {len(diffs)} differences")
            for pos, b1, b2 in diffs[:3]:
                print(f"    Position {pos+1}: {b1} -> {b2}")

    unique = {}
    for name, seq in seqs.items():
        if seq not in unique.values():
            unique[name] = seq
    print(f"\nNumber of unique sequences: {len(unique)}")
    print(f"Unique sequence representatives: {list(unique.keys())}")

if __name__ == "__main__":
    if len(sys.argv) == 3:
        compare_two(sys.argv[1], sys.argv[2])
    else:
        compare_multiple(sys.argv[1:])
