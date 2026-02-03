#!/usr/bin/env python3
"""
Analyze chloroplast genome orientations using ndhD and rbcL gene orientation.
"""
import subprocess
import os
import sys
from pathlib import Path
import re

BASE_DIR = Path("/media/data/projects/iva_phylogeny/analysis")
NDHD_FAA = BASE_DIR / "ndhD.faa"
RBCL_FAA = BASE_DIR / "rbcL.faa"

def get_orientation(query_faa, subject_fasta):
    """Use tblastn to determine gene orientation. Returns 'forward' if sstart < send, 'reverse' otherwise."""
    try:
        result = subprocess.run(
            ["tblastn", "-query", str(query_faa), "-subject", str(subject_fasta),
             "-outfmt", "6 sstart send evalue bitscore", "-max_target_seqs", "1"],
            capture_output=True, text=True, timeout=30
        )
        if result.returncode != 0 or not result.stdout.strip():
            return None
        # Take the first (best) hit
        line = result.stdout.strip().split('\n')[0]
        parts = line.split('\t')
        sstart, send = int(parts[0]), int(parts[1])
        return "forward" if sstart < send else "reverse"
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        return None

def extract_sample_name(pt_dir):
    """Extract a meaningful sample name from the _pt directory path."""
    pt_dir = Path(pt_dir)
    dir_name = pt_dir.name  # e.g., embplant_pt or hayesiana_CA2_93_pt
    parent_name = pt_dir.parent.name  # e.g., 93_STURM_251201adq30ft
    grandparent = pt_dir.parent.parent.name  # e.g., analysis, frutescens, imbricata

    # Extract sample number from parent directory
    sample_match = re.match(r'^(\d+[-\d]*)', parent_name)
    sample_num = sample_match.group(1) if sample_match else parent_name.split('_')[0]

    # Check if the _pt directory has species info
    if dir_name != "embplant_pt":
        # Extract species from dir name like "hayesiana_CA2_93_pt"
        species_match = re.match(r'([a-z]+)_', dir_name)
        species = species_match.group(1) if species_match else ""
    else:
        species = ""

    # Add subdirectory context if applicable
    if grandparent in ["frutescens", "imbricata"]:
        species = grandparent if not species else species

    if species:
        return f"{species}_{sample_num}"
    return sample_num

def find_fasta_files(pt_dir):
    """Find fasta files in a _pt directory, returning pairs of potential isoforms."""
    pt_path = Path(pt_dir)
    fastas = list(pt_path.glob("*.fasta"))

    # Filter to get main isoform files (exclude revcomp files)
    main_fastas = [f for f in fastas if "revcomp" not in f.name]

    # Look for labeled forward/reverse files first
    forward_files = [f for f in main_fastas if "forward_isoform" in f.name]
    reverse_files = [f for f in main_fastas if "reverse_isoform" in f.name]

    if forward_files and reverse_files:
        return {"forward": forward_files[0], "reverse": reverse_files[0], "labeled": True}

    # Otherwise, look for graph1.1 and graph1.2 files
    graph1_files = sorted([f for f in main_fastas if "graph1.1" in f.name or "graph1.2" in f.name],
                          key=lambda x: x.name)

    # Handle cases with multiple isoforms (like sample 67)
    if len(graph1_files) >= 2:
        # Take only the first two (1.1 and 1.2)
        graph1_files = [f for f in graph1_files if ".graph1.1." in f.name or ".graph1.2." in f.name]
        if len(graph1_files) >= 2:
            return {"1": graph1_files[0], "2": graph1_files[1], "labeled": False}

    return None

def main():
    # Find all _pt directories
    pt_dirs = []
    for root, dirs, files in os.walk(BASE_DIR):
        for d in dirs:
            if d.endswith("_pt"):
                pt_dirs.append(os.path.join(root, d))

    pt_dirs.sort()

    print("Sample\tDirectory\tndhD_graph1.1\trbcL_graph1.1\tndhD_graph1.2\trbcL_graph1.2\tForward_isoform\tFasta_file")

    for pt_dir in pt_dirs:
        sample_name = extract_sample_name(pt_dir)
        fastas = find_fasta_files(pt_dir)

        if not fastas:
            print(f"{sample_name}\t{pt_dir}\tNO_FASTA\t-\t-\t-\t-\t-", file=sys.stderr)
            continue

        if fastas.get("labeled"):
            # Already labeled
            forward_file = fastas["forward"]
            ndhd_fwd = get_orientation(NDHD_FAA, forward_file)
            rbcl_fwd = get_orientation(RBCL_FAA, forward_file)
            print(f"{sample_name}\t{pt_dir}\t{ndhd_fwd or 'NA'}\t{rbcl_fwd or 'NA'}\t-\t-\tpre-labeled\t{forward_file}")
        else:
            # Need to determine orientation
            fasta1 = fastas.get("1")
            fasta2 = fastas.get("2")

            ndhd_1 = get_orientation(NDHD_FAA, fasta1) if fasta1 else None
            rbcl_1 = get_orientation(RBCL_FAA, fasta1) if fasta1 else None
            ndhd_2 = get_orientation(NDHD_FAA, fasta2) if fasta2 else None
            rbcl_2 = get_orientation(RBCL_FAA, fasta2) if fasta2 else None

            # Determine which is forward based on ndhD (primary) and rbcL (secondary confirmation)
            forward_isoform = None
            forward_file = None

            if ndhd_1 == "forward":
                forward_isoform = "1"
                forward_file = fasta1
            elif ndhd_2 == "forward":
                forward_isoform = "2"
                forward_file = fasta2

            print(f"{sample_name}\t{pt_dir}\t{ndhd_1 or 'NA'}\t{rbcl_1 or 'NA'}\t{ndhd_2 or 'NA'}\t{rbcl_2 or 'NA'}\t{forward_isoform or 'UNKNOWN'}\t{forward_file or 'NA'}")

if __name__ == "__main__":
    main()
