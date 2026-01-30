#!/bin/bash
# Concatenate chloroplast genome isoforms by ndhD orientation
# Created: 2026-01-30
# Project: Iva frutescens chloroplast phylogeny

# Define forward and reverse isoform assignments based on ndhD orientation
# Forward in isoform 1: 36, 37, 63, 71, 94
# Forward in isoform 2: 49, 58, 68, 86

declare -A forward_isoform
forward_isoform[36]=1
forward_isoform[37]=1
forward_isoform[49]=2
forward_isoform[58]=2
forward_isoform[63]=1
forward_isoform[68]=2
forward_isoform[71]=1
forward_isoform[86]=2
forward_isoform[94]=1

WORKDIR="/media/data/projects/iva_phylogeny/analysis/frutescens"
FORWARD_OUT="$WORKDIR/chloroplast_forward.fasta"
REVERSE_OUT="$WORKDIR/chloroplast_reverse.fasta"

# Clear output files
> "$FORWARD_OUT"
> "$REVERSE_OUT"

for sample_dir in "$WORKDIR"/*/embplant_pt; do
    sample=$(dirname "$sample_dir" | xargs basename | cut -d'_' -f1)

    fwd_iso=${forward_isoform[$sample]}

    for fasta in "$sample_dir"/*.fasta; do
        isoform=$(basename "$fasta" | grep -oP 'graph1\.\K[12]')

        # Get sequence (skip header, join lines)
        seq=$(grep -v "^>" "$fasta" | tr -d '\n')

        if [ "$isoform" -eq "$fwd_iso" ]; then
            echo ">frutescens_${sample}_forward" >> "$FORWARD_OUT"
            echo "$seq" >> "$FORWARD_OUT"
        else
            echo ">frutescens_${sample}_reverse" >> "$REVERSE_OUT"
            echo "$seq" >> "$REVERSE_OUT"
        fi
    done
done

echo "Created: $FORWARD_OUT"
echo "Created: $REVERSE_OUT"
