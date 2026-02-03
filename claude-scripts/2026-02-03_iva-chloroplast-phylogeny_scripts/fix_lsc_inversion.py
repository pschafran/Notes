#!/usr/bin/env python3
"""
Fix LSC inversions by reverse complementing the LSC region.
Uses BLASTN to find the exact inversion boundary for each sample.
"""
import subprocess
import os

def reverse_complement(seq):
    """Return reverse complement of a DNA sequence."""
    complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G',
                  'a': 't', 't': 'a', 'g': 'c', 'c': 'g',
                  'N': 'N', 'n': 'n'}
    return ''.join(complement.get(base, base) for base in reversed(seq))

def read_fasta_single(filepath):
    """Read a single-sequence fasta file."""
    with open(filepath, 'r') as f:
        header = f.readline().strip()
        sequence = ''.join(line.strip() for line in f if not line.startswith('>'))
    return header, sequence

def write_fasta(filepath, header, sequence, line_width=80):
    """Write a fasta file with wrapped sequence."""
    with open(filepath, 'w') as f:
        f.write(header + '\n')
        for i in range(0, len(sequence), line_width):
            f.write(sequence[i:i+line_width] + '\n')

def find_lsc_boundary(sample_fasta, ref_fasta):
    """Use BLASTN to find where the LSC inversion ends (start of IR)."""
    result = subprocess.run(
        ['blastn', '-query', sample_fasta, '-subject', ref_fasta,
         '-outfmt', '6 qstart qend sstart send pident length'],
        capture_output=True, text=True
    )

    # Parse BLAST results to find the transition point
    # Look for where alignment switches from reverse to forward
    hits = []
    for line in result.stdout.strip().split('\n'):
        if not line:
            continue
        parts = line.split('\t')
        qstart, qend = int(parts[0]), int(parts[1])
        sstart, send = int(parts[2]), int(parts[3])
        length = int(parts[5])
        is_reverse = sstart > send
        hits.append({
            'qstart': qstart, 'qend': qend,
            'sstart': sstart, 'send': send,
            'length': length, 'is_reverse': is_reverse
        })

    # Sort by query position
    hits.sort(key=lambda x: x['qstart'])

    # Find the boundary where reverse hits end and forward hits begin
    # This should be around 83-84kb (end of LSC)
    lsc_end = None
    for i, hit in enumerate(hits):
        if hit['length'] > 10000:  # Only consider major hits
            if not hit['is_reverse'] and hit['qstart'] > 80000:
                # This is probably the start of IR (forward alignment)
                lsc_end = hit['qstart'] - 1
                break

    if lsc_end is None:
        # Fall back to looking at the pattern
        reverse_hits = [h for h in hits if h['is_reverse'] and h['length'] > 5000]
        if reverse_hits:
            lsc_end = max(h['qend'] for h in reverse_hits)

    return lsc_end

def main():
    base_dir = "/media/data/projects/iva_phylogeny/analysis"
    ref_fasta = "/tmp/claude/ref_36.fasta"

    # Samples with LSC inversions and their source files
    inverted_samples = {
        "01": f"{base_dir}/01_EREMID_240811/embplant_pt/embplant_pt.K115.complete.graph1.1.path_sequence.fasta",
        "48": f"{base_dir}/48_EREMID_240811/embplant_pt/embplant_pt.K115.complete.graph1.1.path_sequence.fasta",
        "frutescens_58": f"{base_dir}/frutescens/58_STURM_251201adq30ft/embplant_pt/embplant_pt.K115.complete.graph1.2.path_sequence.fasta",
        "hayesiana_92": f"{base_dir}/92_STURM_251201adq30ft/hayesiana_CA1_92_pt/embplant_pt.K115.complete.graph1.2.path_sequence.forward_isoform.fasta",
    }

    # Create reference file
    subprocess.run(['seqkit', 'grep', '-p', 'Iva_frutescens_36',
                   f'{base_dir}/chloroplast_forward_combined.fasta',
                   '-o', ref_fasta], check=True)

    fixed_files = {}

    for sample_name, source_file in inverted_samples.items():
        print(f"\n{'='*60}")
        print(f"Processing {sample_name}")
        print(f"Source: {source_file}")

        # Read the original sequence
        header, sequence = read_fasta_single(source_file)
        print(f"Original length: {len(sequence)} bp")

        # Write temp file for BLASTN
        temp_fasta = f"/tmp/claude/temp_{sample_name}.fasta"
        write_fasta(temp_fasta, header, sequence)

        # Find LSC boundary
        lsc_end = find_lsc_boundary(temp_fasta, ref_fasta)
        print(f"Detected LSC end position: {lsc_end}")

        if lsc_end is None or lsc_end < 70000 or lsc_end > 90000:
            print(f"WARNING: LSC boundary detection may be incorrect, using default 83700")
            lsc_end = 83700

        # Extract and reverse complement LSC
        lsc_region = sequence[:lsc_end]
        rest_of_genome = sequence[lsc_end:]

        print(f"LSC region: 1-{lsc_end} ({len(lsc_region)} bp)")
        print(f"Rest of genome: {lsc_end+1}-{len(sequence)} ({len(rest_of_genome)} bp)")

        # Reverse complement LSC
        lsc_revcomp = reverse_complement(lsc_region)

        # Reassemble
        fixed_sequence = lsc_revcomp + rest_of_genome
        print(f"Fixed sequence length: {len(fixed_sequence)} bp")

        # Write fixed file
        output_file = f"/tmp/claude/fixed_{sample_name}.fasta"
        write_fasta(output_file, f">Iva_{sample_name}", fixed_sequence)
        fixed_files[sample_name] = output_file

        print(f"Wrote: {output_file}")

    print(f"\n{'='*60}")
    print("Fixed files created:")
    for name, path in fixed_files.items():
        print(f"  {name}: {path}")

    return fixed_files

if __name__ == "__main__":
    main()
