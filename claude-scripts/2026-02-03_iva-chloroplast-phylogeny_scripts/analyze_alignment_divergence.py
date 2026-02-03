#!/usr/bin/env python3
"""
Analyze alignment for potential misorientation by calculating sliding window divergence.
High divergence regions may indicate inversions or assembly errors.
"""
import sys
from collections import defaultdict

def read_fasta(filepath):
    """Read aligned fasta file into dictionary."""
    sequences = {}
    current_name = None
    current_seq = []

    with open(filepath, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if current_name:
                    sequences[current_name] = ''.join(current_seq)
                current_name = line[1:].split()[0]
                current_seq = []
            else:
                current_seq.append(line.upper())
        if current_name:
            sequences[current_name] = ''.join(current_seq)

    return sequences

def calculate_divergence(seq1, seq2, start, end):
    """Calculate percent divergence between two aligned sequences in a region."""
    mismatches = 0
    valid_positions = 0

    for i in range(start, min(end, len(seq1), len(seq2))):
        c1, c2 = seq1[i], seq2[i]
        # Skip gaps and ambiguous bases
        if c1 in 'ACGT' and c2 in 'ACGT':
            valid_positions += 1
            if c1 != c2:
                mismatches += 1

    if valid_positions == 0:
        return None
    return (mismatches / valid_positions) * 100

def sliding_window_analysis(sequences, window_size=10000, step=5000, reference=None):
    """Perform sliding window divergence analysis."""
    if reference is None:
        # Use first sequence as reference
        reference = list(sequences.keys())[0]

    ref_seq = sequences[reference]
    alignment_len = len(ref_seq)

    results = {}

    for name, seq in sequences.items():
        if name == reference:
            continue

        windows = []
        for start in range(0, alignment_len - window_size, step):
            end = start + window_size
            div = calculate_divergence(ref_seq, seq, start, end)
            if div is not None:
                windows.append((start, end, div))

        results[name] = windows

    return results, reference

def find_high_divergence_regions(results, threshold=5.0):
    """Find sequences and regions with abnormally high divergence."""
    problematic = {}

    for name, windows in results.items():
        high_div_regions = [(s, e, d) for s, e, d in windows if d > threshold]
        if high_div_regions:
            # Calculate overall stats
            all_divs = [d for _, _, d in windows]
            avg_div = sum(all_divs) / len(all_divs) if all_divs else 0
            max_div = max(all_divs) if all_divs else 0

            # Check if there's a large contiguous region of high divergence
            if len(high_div_regions) >= 3:  # At least 3 consecutive windows
                problematic[name] = {
                    'avg_divergence': avg_div,
                    'max_divergence': max_div,
                    'high_div_regions': high_div_regions,
                    'num_high_windows': len(high_div_regions)
                }

    return problematic

def main():
    alignment_file = "/media/data/projects/iva_phylogeny/analysis/chloroplast_forward_aligned.fasta"

    print("Reading alignment...")
    sequences = read_fasta(alignment_file)
    print(f"Loaded {len(sequences)} sequences, alignment length: {len(list(sequences.values())[0])} bp")

    # Use Iva_frutescens_36 as reference (known good from previous session)
    reference = "Iva_frutescens_36"
    if reference not in sequences:
        reference = list(sequences.keys())[0]

    print(f"\nUsing {reference} as reference")
    print("\nPerforming sliding window analysis (10kb windows, 5kb step)...")

    results, ref = sliding_window_analysis(sequences, window_size=10000, step=5000, reference=reference)

    # Calculate summary stats for each sample
    print("\n" + "="*80)
    print("DIVERGENCE SUMMARY (vs reference)")
    print("="*80)
    print(f"{'Sample':<30} {'Avg Div %':>10} {'Max Div %':>10} {'High Div Windows':>18}")
    print("-"*80)

    all_stats = []
    for name, windows in sorted(results.items()):
        divs = [d for _, _, d in windows]
        avg_div = sum(divs) / len(divs) if divs else 0
        max_div = max(divs) if divs else 0
        high_windows = sum(1 for d in divs if d > 5.0)
        all_stats.append((name, avg_div, max_div, high_windows))
        print(f"{name:<30} {avg_div:>10.3f} {max_div:>10.3f} {high_windows:>18}")

    # Find problematic samples
    print("\n" + "="*80)
    print("POTENTIAL MISORIENTATION ISSUES (>5% divergence regions)")
    print("="*80)

    problematic = find_high_divergence_regions(results, threshold=5.0)

    if not problematic:
        print("No significant misorientation issues detected!")
    else:
        for name, info in sorted(problematic.items(), key=lambda x: -x[1]['max_divergence']):
            print(f"\n{name}:")
            print(f"  Average divergence: {info['avg_divergence']:.2f}%")
            print(f"  Maximum divergence: {info['max_divergence']:.2f}%")
            print(f"  High divergence windows: {info['num_high_windows']}")

            # Show the region spans
            regions = info['high_div_regions']
            if regions:
                min_start = min(r[0] for r in regions)
                max_end = max(r[1] for r in regions)
                print(f"  Affected region: {min_start:,} - {max_end:,} bp")

if __name__ == "__main__":
    main()
