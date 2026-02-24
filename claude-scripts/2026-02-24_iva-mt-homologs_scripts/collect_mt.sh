#!/bin/bash
# collect_mt.sh - Collect and label all mt scaffold sequences for phylogeny paper
#
# Handles 23 usable samples:
#   - 20 samples with single graph1.1.path_sequence.fasta
#   - 3 samples with only repeat_pattern files (axillaris_89, frutescens_36, frutescens_49)
#     → use repeat_pattern1.1 as representative
#
# Output:
#   all_mt_scaffolds.fasta   — combined labelled sequences (>SAMPLEID|original_header)
#   mt_manifest.tsv          — sample, file, n_seqs, total_bp, notes

set -euo pipefail

OUTDIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BASE="/media/data/projects/iva_phylogeny/analysis"

COMBINED="$OUTDIR/all_mt_scaffolds.fasta"
MANIFEST="$OUTDIR/mt_manifest.tsv"

> "$COMBINED"
echo -e "sample_id\tfile\tn_seqs\ttotal_bp\tnotes" > "$MANIFEST"

add_sample() {
    local sample_id="$1"
    local fasta_file="$2"
    local notes="${3:-}"

    if [ ! -f "$fasta_file" ]; then
        echo "ERROR: file not found: $fasta_file" >&2
        echo -e "${sample_id}\t${fasta_file}\tNA\tNA\tFILE_NOT_FOUND" >> "$MANIFEST"
        return 1
    fi

    local n_seqs
    n_seqs=$(grep -c "^>" "$fasta_file")
    local total_bp
    total_bp=$(grep -v "^>" "$fasta_file" | tr -d '\n\r' | wc -c)

    # Relabel headers: >SAMPLEID|original_header
    awk -v sid="$sample_id" '
        /^>/ { print ">" sid "|" substr($0, 2) }
        !/^>/ { print }
    ' "$fasta_file" >> "$COMBINED"

    echo -e "${sample_id}\t${fasta_file}\t${n_seqs}\t${total_bp}\t${notes}" >> "$MANIFEST"
    printf "  %-20s  %3d seqs  %9d bp  %s\n" "$sample_id" "$n_seqs" "$total_bp" "$notes"
}

echo "=== Collecting mt scaffold sequences ==="
echo ""

echo "--- EREMID batch (single graph1.1 file) ---"
add_sample "angustifolia_46" "$BASE/46_EREMID_240811/embplant_mt/embplant_mt.K115.scaffolds.graph1.1.path_sequence.fasta"
add_sample "annua_01"        "$BASE/annua/01_EREMID_240811/embplant_mt/embplant_mt.K115.scaffolds.graph1.1.path_sequence.fasta"
add_sample "cultigen_133"    "$BASE/cultigen/13-3-33_EREMID_240811/embplant_mt/embplant_mt.K115.scaffolds.graph1.1.path_sequence.fasta"
add_sample "frutescens_55"   "$BASE/frutescens/55_EREMID_240811/embplant_mt/embplant_mt.K115.scaffolds.graph1.1.path_sequence.fasta"
add_sample "microcephala_57" "$BASE/microcephala/57_EREMID_240811/embplant_mt/embplant_mt.K115.scaffolds.graph1.1.path_sequence.fasta"
add_sample "texensis_48"     "$BASE/texensis/48_EREMID_240811/embplant_mt/embplant_mt.K115.scaffolds.graph1.1.path_sequence.fasta"

echo ""
echo "--- STURM batch (single graph1.1 file) ---"
# Species-prefixed mt directories
add_sample "angustifolia_77" "$BASE/angustifolia/77_STURM_251201adq30ft/angustifolia_LA_77_mt/embplant_mt.K115.scaffolds.graph1.1.path_sequence.fasta"
add_sample "asperifolia_60"  "$BASE/asperifolia/60_STURM_251201adq30ft/asperifolia_FL_60_mt/embplant_mt.K115.scaffolds.graph1.1.path_sequence.fasta"
add_sample "hayesiana_92"    "$BASE/hayesiana/92_STURM_251201adq30ft/hayesiana_CA1_92_mt/embplant_mt.K115.scaffolds.graph1.1.path_sequence.fasta"
add_sample "hayesiana_93"    "$BASE/hayesiana/93_STURM_251201adq30ft/hayesiana_CA2_93_mt/embplant_mt.K115.scaffolds.graph1.1.path_sequence.fasta"
add_sample "xanthifolia_64"  "$BASE/xanthifolia/64_STURM_251201adq30ft/xanthiifolia_CO1_64_mt/embplant_mt.K115.scaffolds.graph1.1.path_sequence.fasta"
# embplant_mt directories
add_sample "frutescens_58"   "$BASE/frutescens/58_STURM_251201adq30ft/embplant_mt/embplant_mt.K85.scaffolds.graph1.1.path_sequence.fasta"
add_sample "frutescens_63"   "$BASE/frutescens/63_STURM_251201adq30ft/embplant_mt/embplant_mt.K115.scaffolds.graph1.1.path_sequence.fasta"
add_sample "frutescens_68"   "$BASE/frutescens/68_STURM_251201adq30ft/embplant_mt/embplant_mt.K115.scaffolds.graph1.1.path_sequence.fasta"
add_sample "frutescens_71"   "$BASE/frutescens/71_STURM_251201adq30ft/embplant_mt/embplant_mt.K115.scaffolds.graph1.1.path_sequence.fasta"
add_sample "frutescens_86"   "$BASE/frutescens/86_STURM_251201adq30ft/embplant_mt/embplant_mt.K115.scaffolds.graph1.1.path_sequence.fasta"
add_sample "frutescens_94"   "$BASE/frutescens/94_STURM_251201adq30ft/embplant_mt/embplant_mt.K115.scaffolds.graph1.1.path_sequence.fasta"
add_sample "imbricata_67"    "$BASE/imbricata/67_STURM_251201adq30ft/embplant_mt/embplant_mt.K115.scaffolds.graph1.1.path_sequence.fasta"
add_sample "imbricata_72"    "$BASE/imbricata/72_STURM_251201adq30ft/embplant_mt/embplant_mt.K115.scaffolds.graph1.1.path_sequence.fasta"
add_sample "imbricata_85"    "$BASE/imbricata/85_STURM_251201adq30ft/embplant_mt/embplant_mt.K115.scaffolds.graph1.1.path_sequence.fasta"

echo ""
echo "--- Repeat-pattern samples (using repeat_pattern1.1 as representative) ---"
add_sample "axillaris_89"    "$BASE/axillaris/89_STURM_251201adq30ft/embplant_mt/embplant_mt.K115.scaffolds.graph1.repeat_pattern1.1.path_sequence.fasta" "repeat_pattern_representative"
add_sample "frutescens_36"   "$BASE/frutescens/36_STURM_251201adq30ft/embplant_mt/embplant_mt.K115.scaffolds.graph1.repeat_pattern1.1.path_sequence.fasta" "repeat_pattern_representative"
add_sample "frutescens_49"   "$BASE/frutescens/49_STURM_251201adq30ft/embplant_mt/embplant_mt.K115.scaffolds.graph1.repeat_pattern1.1.path_sequence.fasta" "repeat_pattern_representative"

echo ""
echo "=== Done ==="
echo "Combined FASTA : $COMBINED"
echo "Manifest       : $MANIFEST"
TOTAL_SEQS=$(grep -c "^>" "$COMBINED")
TOTAL_SAMPLES=$(grep -v "^sample_id" "$MANIFEST" | grep -v "FILE_NOT_FOUND" | wc -l)
echo "Samples        : $TOTAL_SAMPLES"
echo "Total sequences: $TOTAL_SEQS"
