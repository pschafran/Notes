#!/usr/bin/env python3
"""
Fix circular rotation of chloroplast genome to match reference start position.
Uses rbcL position as anchor point.
"""

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

def rotate_sequence(sequence, new_start):
    """Rotate a circular sequence to start at a new position."""
    return sequence[new_start:] + sequence[:new_start]

def main():
    # Reference rbcL start position (from tblastn): 28967
    # imbricata_69 rbcL start position: 108532

    # Read imbricata_69
    header, sequence = read_fasta_single("/tmp/claude/imbricata_69.fasta")
    print(f"Original sequence length: {len(sequence)}")

    # Calculate rotation
    # We want rbcL to end up at approximately position 28967
    # Currently rbcL is at 108532
    # Rotation amount = 108532 - 28967 = 79565

    ref_rbcl_pos = 28967
    current_rbcl_pos = 108532
    rotation = current_rbcl_pos - ref_rbcl_pos

    print(f"Reference rbcL position: {ref_rbcl_pos}")
    print(f"Current rbcL position: {current_rbcl_pos}")
    print(f"Rotation amount: {rotation}")

    # Rotate the sequence
    rotated = rotate_sequence(sequence, rotation)
    print(f"Rotated sequence length: {len(rotated)}")

    # Write output
    output_file = "/tmp/claude/fixed_imbricata_69.fasta"
    write_fasta(output_file, ">Iva_imbricata_69", rotated)
    print(f"Wrote: {output_file}")

    # Verify by checking rbcL position in rotated sequence
    # (Would need to run tblastn to confirm)

if __name__ == "__main__":
    main()
