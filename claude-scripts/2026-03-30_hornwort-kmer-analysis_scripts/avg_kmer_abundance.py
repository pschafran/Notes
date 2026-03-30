#!/usr/bin/env python3
"""
Calculate average k-mer abundance, GC content, and N (hard-masked repeat) content per contig.
Usage: python avg_kmer_abundance.py -f genome.fasta -k kmer_counts.tsv -s 21 -o output.tsv
"""
import argparse
import sys

def reverse_complement(seq):
    comp = str.maketrans('ACGTacgt', 'TGCAtgca')
    return seq.translate(comp)[::-1]

def canonical(kmer):
    rc = reverse_complement(kmer)
    return kmer if kmer <= rc else rc

def parse_fasta(fasta_path):
    contigs = {}
    name = None
    seq_parts = []
    with open(fasta_path) as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if name:
                    contigs[name] = ''.join(seq_parts).upper()
                name = line[1:].split()[0]
                seq_parts = []
            else:
                seq_parts.append(line)
        if name:
            contigs[name] = ''.join(seq_parts).upper()
    return contigs

def load_kmer_counts(tsv_path):
    print("Loading k-mer counts...", file=sys.stderr)
    counts = {}
    with open(tsv_path) as f:
        for line in f:
            kmer, count = line.strip().split('\t')
            counts[kmer] = int(count)
    print(f"  Loaded {len(counts):,} k-mers", file=sys.stderr)
    return counts

def avg_kmer_abundance(seq, k, kmer_counts):
    total = 0
    valid_kmers = 0
    n_kmers = 0
    for i in range(len(seq) - k + 1):
        kmer = seq[i:i+k]
        if 'N' in kmer:
            n_kmers += 1
            continue
        ck = canonical(kmer)
        total += kmer_counts.get(ck, 0)
        valid_kmers += 1
    return total / valid_kmers if valid_kmers > 0 else 0, valid_kmers, n_kmers

def gc_content(seq):
    """Calculate GC content, excluding N bases."""
    non_n = sum(1 for b in seq if b != 'N')
    if non_n == 0:
        return 0.0
    gc = sum(1 for b in seq if b in ('G', 'C'))
    return gc / non_n

def n_stats(seq):
    """Count N bases and their proportion (hard-masked repeats)."""
    n_bases = seq.count('N')
    return n_bases, n_bases / len(seq) if len(seq) > 0 else 0.0

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('-f', '--fasta',     required=True, help='Genome FASTA')
    ap.add_argument('-k', '--kmers',     required=True, help='Jellyfish dump TSV (kmer<TAB>count)')
    ap.add_argument('-s', '--kmer_size', type=int, default=21, help='K-mer size (default: 21)')
    ap.add_argument('-o', '--output',    required=True, help='Output TSV')
    args = ap.parse_args()

    kmer_counts = load_kmer_counts(args.kmers)
    contigs = parse_fasta(args.fasta)
    print(f"  Loaded {len(contigs):,} contigs", file=sys.stderr)

    with open(args.output, 'w') as out:
        out.write("contig\tseq_len\tnum_kmers\tn_kmers\tavg_kmer_abundance\tgc_content\tn_bases\tn_fraction\n")
        for name, seq in contigs.items():
            avg, valid_kmers, masked_kmers = avg_kmer_abundance(seq, args.kmer_size, kmer_counts)
            gc = gc_content(seq)
            n_bases, n_frac = n_stats(seq)
            out.write(f"{name}\t{len(seq)}\t{valid_kmers}\t{masked_kmers}\t{avg:.4f}\t{gc:.4f}\t{n_bases}\t{n_frac:.4f}\n")
            print(f"  {name}: avg_kmer={avg:.2f}, gc={gc:.4f}, n_frac={n_frac:.4f}", file=sys.stderr)

    print(f"Done. Results written to {args.output}", file=sys.stderr)

if __name__ == '__main__':
    main()
