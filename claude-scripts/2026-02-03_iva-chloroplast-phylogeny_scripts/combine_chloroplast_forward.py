#!/usr/bin/env python3
"""
Combine forward-oriented chloroplast genomes into a single file.
Select unique samples, preferring species-labeled directories over generic embplant_pt.
"""
import os
import re
from pathlib import Path

# Forward isoform files from analysis
forward_isoforms = {
    "01": "/media/data/projects/iva_phylogeny/analysis/01_EREMID_240811/embplant_pt/embplant_pt.K115.complete.graph1.1.path_sequence.fasta",
    "13-3-33": "/media/data/projects/iva_phylogeny/analysis/13-3-33_EREMID_240811/embplant_pt/embplant_pt.K115.complete.graph1.2.path_sequence.fasta",
    "25": "/media/data/projects/iva_phylogeny/analysis/25_EREMID_240811/embplant_pt/embplant_pt.K115.complete.graph1.2.path_sequence.fasta",
    "46": "/media/data/projects/iva_phylogeny/analysis/46_EREMID_240811/embplant_pt/embplant_pt.K115.complete.graph1.1.path_sequence.fasta",
    "48": "/media/data/projects/iva_phylogeny/analysis/48_EREMID_240811/embplant_pt/embplant_pt.K115.complete.graph1.1.path_sequence.fasta",
    "50": "/media/data/projects/iva_phylogeny/analysis/50_EREMID_240811/embplant_pt/embplant_pt.K115.complete.graph1.2.path_sequence.fasta",
    "annua_51": "/media/data/projects/iva_phylogeny/analysis/51_STURM_251201adq30ft/annua_TXS_51_pt/embplant_pt.K115.complete.graph1.2.path_sequence.forward_isoform.fasta",
    "54": "/media/data/projects/iva_phylogeny/analysis/54_EREMID_240811/embplant_pt/embplant_pt.K115.complete.graph1.1.path_sequence.fasta",
    "55": "/media/data/projects/iva_phylogeny/analysis/55_EREMID_240811/embplant_pt/embplant_pt.K115.complete.graph1.2.path_sequence.fasta",
    "57": "/media/data/projects/iva_phylogeny/analysis/57_EREMID_240811/embplant_pt/embplant_pt.K115.complete.graph1.2.path_sequence.fasta",
    "asperifolia_60": "/media/data/projects/iva_phylogeny/analysis/60_STURM_251201adq30ft/asperifolia_FL_60_pt/embplant_pt.K115.complete.graph1.1.path_sequence.forward_isoform.fasta",
    "xanthiifolia_64": "/media/data/projects/iva_phylogeny/analysis/64_STURM_251201adq30ft/xanthiifolia_CO1_64_pt/embplant_pt.K115.complete.graph1.2.path_sequence.forward_isoform.fasta",
    "frutescens_66": "/media/data/projects/iva_phylogeny/analysis/66_STURM_251201adq30ft/frutescens_FLE_66_pt/embplant_pt.K115.complete.graph1.2.path_sequence.forward_isoform.fasta",
    "angustifolia_77": "/media/data/projects/iva_phylogeny/analysis/77_STURM_251201adq30ft/angustifolia_LA_77_pt/embplant_pt.K115.complete.graph1.1.path_sequence.forward_isoform.fasta",
    "frutescens_83": "/media/data/projects/iva_phylogeny/analysis/83_STURM_251201adq30ft/frutescens_MA_83_pt/embplant_pt.K115.complete.graph1.1.path_sequence.forward_isoform.fasta",
    "xanthifolia_88": "/media/data/projects/iva_phylogeny/analysis/88_STURM_251201adq30ft/xanthifolia_CO2_88_pt/embplant_pt.K115.complete.graph1.2.path_sequence.forward_isoform.fasta",
    "axillaris_89": "/media/data/projects/iva_phylogeny/analysis/89_STURM_251201adq30ft/axillaris_co_89_pt/embplant_pt.K115.complete.graph1.2.path_sequence.forward_isoform.fasta",
    "axillaris_91": "/media/data/projects/iva_phylogeny/analysis/91_STURM_251201adq30ft/axillaris_CA_91_pt/embplant_pt.K115.complete.graph1.2.path_sequence.forward_isoform.fasta",
    "hayesiana_92": "/media/data/projects/iva_phylogeny/analysis/92_STURM_251201adq30ft/hayesiana_CA1_92_pt/embplant_pt.K115.complete.graph1.2.path_sequence.forward_isoform.fasta",
    "hayesiana_93": "/media/data/projects/iva_phylogeny/analysis/93_STURM_251201adq30ft/hayesiana_CA2_93_pt/embplant_pt.K115.complete.graph1.1.path_sequence.forward_isoform.fasta",
    "frutescens_36": "/media/data/projects/iva_phylogeny/analysis/frutescens/36_STURM_251201adq30ft/embplant_pt/embplant_pt.K115.complete.graph1.1.path_sequence.fasta",
    "frutescens_37": "/media/data/projects/iva_phylogeny/analysis/frutescens/37_STURM_251201adq30ft/embplant_pt/embplant_pt.K115.complete.graph1.1.path_sequence.fasta",
    "frutescens_49": "/media/data/projects/iva_phylogeny/analysis/frutescens/49_STURM_251201adq30ft/embplant_pt/embplant_pt.K115.complete.graph1.2.path_sequence.fasta",
    "frutescens_58": "/media/data/projects/iva_phylogeny/analysis/frutescens/58_STURM_251201adq30ft/embplant_pt/embplant_pt.K115.complete.graph1.2.path_sequence.fasta",
    "frutescens_63": "/media/data/projects/iva_phylogeny/analysis/frutescens/63_STURM_251201adq30ft/embplant_pt/embplant_pt.K115.complete.graph1.1.path_sequence.fasta",
    "frutescens_68": "/media/data/projects/iva_phylogeny/analysis/frutescens/68_STURM_251201adq30ft/embplant_pt/embplant_pt.K115.complete.graph1.2.path_sequence.fasta",
    "frutescens_71": "/media/data/projects/iva_phylogeny/analysis/frutescens/71_STURM_251201adq30ft/embplant_pt/embplant_pt.K115.complete.graph1.1.path_sequence.fasta",
    "frutescens_86": "/media/data/projects/iva_phylogeny/analysis/frutescens/86_STURM_251201adq30ft/embplant_pt/embplant_pt.K115.complete.graph1.2.path_sequence.fasta",
    "frutescens_94": "/media/data/projects/iva_phylogeny/analysis/frutescens/94_STURM_251201adq30ft/embplant_pt/embplant_pt.K115.complete.graph1.1.path_sequence.fasta",
    "imbricata_61": "/media/data/projects/iva_phylogeny/analysis/imbricata/61_STURM_251201adq30ft/embplant_pt/embplant_pt.K115.complete.graph1.2.path_sequence.fasta",
    "imbricata_67": "/media/data/projects/iva_phylogeny/analysis/imbricata/67_STURM_251201adq30ft/embplant_pt/embplant_pt.K55.complete.graph1.1.repeat_pattern1.path_sequence.fasta",
    "imbricata_69": "/media/data/projects/iva_phylogeny/analysis/imbricata/69_STURM_251201adq30ft/embplant_pt/embplant_pt.K115.scaffolds.graph1.1.path_sequence.fasta",
    "imbricata_72": "/media/data/projects/iva_phylogeny/analysis/imbricata/72_STURM_251201adq30ft/embplant_pt/embplant_pt.K115.complete.graph1.2.path_sequence.fasta",
    "imbricata_85": "/media/data/projects/iva_phylogeny/analysis/imbricata/85_STURM_251201adq30ft/embplant_pt/embplant_pt.K115.complete.graph1.2.path_sequence.fasta",
}

def read_fasta_sequence(filepath):
    """Read a single sequence from a fasta file."""
    with open(filepath, 'r') as f:
        header = f.readline().strip()
        sequence = ""
        for line in f:
            if line.startswith('>'):
                break
            sequence += line.strip()
    return header, sequence

def main():
    output_file = "/media/data/projects/iva_phylogeny/analysis/chloroplast_forward_combined.fasta"

    with open(output_file, 'w') as out:
        for sample_name, fasta_path in sorted(forward_isoforms.items()):
            if not os.path.exists(fasta_path):
                print(f"WARNING: File not found: {fasta_path}")
                continue

            header, sequence = read_fasta_sequence(fasta_path)
            # Write with new sample name
            out.write(f">Iva_{sample_name}\n")
            # Write sequence in 80-character lines
            for i in range(0, len(sequence), 80):
                out.write(sequence[i:i+80] + "\n")

            print(f"Added: Iva_{sample_name} ({len(sequence)} bp)")

    print(f"\nOutput written to: {output_file}")
    print(f"Total samples: {len(forward_isoforms)}")

if __name__ == "__main__":
    main()
