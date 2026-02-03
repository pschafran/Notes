#!/usr/bin/env python3
"""
Regenerate combined chloroplast file with fixed sequences.
"""
import os

BASE_DIR = "/media/data/projects/iva_phylogeny/analysis"

# Fixed sequences (from /tmp/claude/)
fixed_sequences = {
    "Iva_01": "/tmp/claude/fixed_01.fasta",
    "Iva_48": "/tmp/claude/fixed_48.fasta",
    "Iva_frutescens_58": "/tmp/claude/fixed_frutescens_58.fasta",
    "Iva_hayesiana_92": "/tmp/claude/fixed_hayesiana_92.fasta",
    "Iva_imbricata_69": "/tmp/claude/fixed_imbricata_69.fasta",
}

# Original forward isoform files for unchanged samples
original_files = {
    "Iva_13-3-33": f"{BASE_DIR}/13-3-33_EREMID_240811/embplant_pt/embplant_pt.K115.complete.graph1.2.path_sequence.fasta",
    "Iva_25": f"{BASE_DIR}/25_EREMID_240811/embplant_pt/embplant_pt.K115.complete.graph1.2.path_sequence.fasta",
    "Iva_46": f"{BASE_DIR}/46_EREMID_240811/embplant_pt/embplant_pt.K115.complete.graph1.1.path_sequence.fasta",
    "Iva_50": f"{BASE_DIR}/50_EREMID_240811/embplant_pt/embplant_pt.K115.complete.graph1.2.path_sequence.fasta",
    "Iva_annua_51": f"{BASE_DIR}/51_STURM_251201adq30ft/annua_TXS_51_pt/embplant_pt.K115.complete.graph1.2.path_sequence.forward_isoform.fasta",
    "Iva_54": f"{BASE_DIR}/54_EREMID_240811/embplant_pt/embplant_pt.K115.complete.graph1.1.path_sequence.fasta",
    "Iva_55": f"{BASE_DIR}/55_EREMID_240811/embplant_pt/embplant_pt.K115.complete.graph1.2.path_sequence.fasta",
    "Iva_57": f"{BASE_DIR}/57_EREMID_240811/embplant_pt/embplant_pt.K115.complete.graph1.2.path_sequence.fasta",
    "Iva_asperifolia_60": f"{BASE_DIR}/60_STURM_251201adq30ft/asperifolia_FL_60_pt/embplant_pt.K115.complete.graph1.1.path_sequence.forward_isoform.fasta",
    "Iva_xanthiifolia_64": f"{BASE_DIR}/64_STURM_251201adq30ft/xanthiifolia_CO1_64_pt/embplant_pt.K115.complete.graph1.2.path_sequence.forward_isoform.fasta",
    "Iva_frutescens_66": f"{BASE_DIR}/66_STURM_251201adq30ft/frutescens_FLE_66_pt/embplant_pt.K115.complete.graph1.2.path_sequence.forward_isoform.fasta",
    "Iva_angustifolia_77": f"{BASE_DIR}/77_STURM_251201adq30ft/angustifolia_LA_77_pt/embplant_pt.K115.complete.graph1.1.path_sequence.forward_isoform.fasta",
    "Iva_frutescens_83": f"{BASE_DIR}/83_STURM_251201adq30ft/frutescens_MA_83_pt/embplant_pt.K115.complete.graph1.1.path_sequence.forward_isoform.fasta",
    "Iva_xanthifolia_88": f"{BASE_DIR}/88_STURM_251201adq30ft/xanthifolia_CO2_88_pt/embplant_pt.K115.complete.graph1.2.path_sequence.forward_isoform.fasta",
    "Iva_axillaris_89": f"{BASE_DIR}/89_STURM_251201adq30ft/axillaris_co_89_pt/embplant_pt.K115.complete.graph1.2.path_sequence.forward_isoform.fasta",
    "Iva_axillaris_91": f"{BASE_DIR}/91_STURM_251201adq30ft/axillaris_CA_91_pt/embplant_pt.K115.complete.graph1.2.path_sequence.forward_isoform.fasta",
    "Iva_hayesiana_93": f"{BASE_DIR}/93_STURM_251201adq30ft/hayesiana_CA2_93_pt/embplant_pt.K115.complete.graph1.1.path_sequence.forward_isoform.fasta",
    "Iva_frutescens_36": f"{BASE_DIR}/frutescens/36_STURM_251201adq30ft/embplant_pt/embplant_pt.K115.complete.graph1.1.path_sequence.fasta",
    "Iva_frutescens_37": f"{BASE_DIR}/frutescens/37_STURM_251201adq30ft/embplant_pt/embplant_pt.K115.complete.graph1.1.path_sequence.fasta",
    "Iva_frutescens_49": f"{BASE_DIR}/frutescens/49_STURM_251201adq30ft/embplant_pt/embplant_pt.K115.complete.graph1.2.path_sequence.fasta",
    "Iva_frutescens_63": f"{BASE_DIR}/frutescens/63_STURM_251201adq30ft/embplant_pt/embplant_pt.K115.complete.graph1.1.path_sequence.fasta",
    "Iva_frutescens_68": f"{BASE_DIR}/frutescens/68_STURM_251201adq30ft/embplant_pt/embplant_pt.K115.complete.graph1.2.path_sequence.fasta",
    "Iva_frutescens_71": f"{BASE_DIR}/frutescens/71_STURM_251201adq30ft/embplant_pt/embplant_pt.K115.complete.graph1.1.path_sequence.fasta",
    "Iva_frutescens_86": f"{BASE_DIR}/frutescens/86_STURM_251201adq30ft/embplant_pt/embplant_pt.K115.complete.graph1.2.path_sequence.fasta",
    "Iva_frutescens_94": f"{BASE_DIR}/frutescens/94_STURM_251201adq30ft/embplant_pt/embplant_pt.K115.complete.graph1.1.path_sequence.fasta",
    "Iva_imbricata_61": f"{BASE_DIR}/imbricata/61_STURM_251201adq30ft/embplant_pt/embplant_pt.K115.complete.graph1.2.path_sequence.fasta",
    "Iva_imbricata_67": f"{BASE_DIR}/imbricata/67_STURM_251201adq30ft/embplant_pt/embplant_pt.K55.complete.graph1.1.repeat_pattern1.path_sequence.fasta",
    "Iva_imbricata_72": f"{BASE_DIR}/imbricata/72_STURM_251201adq30ft/embplant_pt/embplant_pt.K115.complete.graph1.2.path_sequence.fasta",
    "Iva_imbricata_85": f"{BASE_DIR}/imbricata/85_STURM_251201adq30ft/embplant_pt/embplant_pt.K115.complete.graph1.2.path_sequence.fasta",
}

def read_fasta_sequence(filepath):
    """Read sequence from fasta file (handles both fixed and original format)."""
    with open(filepath, 'r') as f:
        header = f.readline().strip()
        sequence = ''.join(line.strip() for line in f if not line.startswith('>'))
    return sequence

def main():
    output_file = f"{BASE_DIR}/chloroplast_forward_combined_fixed.fasta"

    # Combine all sequences
    all_samples = {**fixed_sequences, **original_files}

    with open(output_file, 'w') as out:
        for sample_name in sorted(all_samples.keys()):
            filepath = all_samples[sample_name]

            if not os.path.exists(filepath):
                print(f"WARNING: Missing file for {sample_name}: {filepath}")
                continue

            sequence = read_fasta_sequence(filepath)

            # Write with sample name header
            out.write(f">{sample_name}\n")
            for i in range(0, len(sequence), 80):
                out.write(sequence[i:i+80] + "\n")

            source = "FIXED" if sample_name in fixed_sequences else "original"
            print(f"Added: {sample_name} ({len(sequence)} bp) [{source}]")

    print(f"\nOutput: {output_file}")
    print(f"Total samples: {len(all_samples)}")

if __name__ == "__main__":
    main()
