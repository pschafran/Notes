#!/bin/bash
# run_edta_liverworts.sh
# Run EDTA on 5 new liverwort genomes sequentially.
# Detects and auto-fixes known EDTA bugs; retries on failure.
# Usage: bash run_edta_liverworts.sh [&]

BASE="/media/data/projects/hornwort_sex_chromosomes/analysis/liverwort_genomes"
MASTER_LOG="$BASE/edta_run_liverworts.log"
THREADS=24
MAX_RETRIES=3

# Activate conda EDTA environment
source ~/miniconda3/etc/profile.d/conda.sh
conda activate EDTA

EDTA_BIN="$(which EDTA.pl)"
if [ -z "$EDTA_BIN" ]; then
    echo "ERROR: EDTA.pl not found after conda activate EDTA. Exiting." | tee -a "$MASTER_LOG"
    exit 1
fi

# Species array: "subdir prefix"
SPECIES=(
    "Blasia_pusilla scaffold_890-Blasia_pusilla"
    "Plagiochasma_appendiculatum scaffold_952-Plagiochasma_appendiculatum"
    "Jungermannia_erectum scaffold_955-Jungermannia_erectum"
    "Conocephalum_conicum scaffold_664-Conocephalum_conicum"
    "Acrolejeunea_sandvicensis scaffold_888-Acrolejeunea_sandvicensis"
)

log() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $1" | tee -a "$MASTER_LOG"
}

# Known bug fix: LTR_retriever writes *.mod.mod.* but EDTA_raw.pl expects *.mod.*
# Fix: copy *.mod.mod.* -> *.mod.* in LTR/ and EDTA.raw/ dirs
fix_ltr_retriever_naming() {
    local dir="$1"
    local prefix="$2"
    local ltr_dir="${dir}/${prefix}.fa.mod.EDTA.raw/LTR"
    local raw_dir="${dir}/${prefix}.fa.mod.EDTA.raw"
    local fixed=0

    for src in $(find "$ltr_dir" "$raw_dir" -maxdepth 1 -name "*.mod.mod.*" 2>/dev/null); do
        dst=$(echo "$src" | sed 's/\.mod\.mod\./\.mod\./')
        if [ "$src" != "$dst" ] && [ ! -f "$dst" ]; then
            cp "$src" "$dst"
            log "  LTR fix: $src -> $dst"
            fixed=1
        fi
    done

    return $fixed
}

# Known bug fix: TIR-Learner TypeError from absl-py incompatibility
# Fix: downgrade absl-py to 0.9.0
fix_absl_py() {
    log "  absl-py fix: downgrading to 0.9.0"
    pip install absl-py==0.9.0 -q
}

run_edta() {
    local subdir="$1"
    local prefix="$2"
    local dir="$BASE/$subdir"
    local genome="${dir}/${prefix}.fa"
    local cds="${dir}/${prefix}.cds.fa"
    local final_output="${dir}/${prefix}.fa.mod.EDTA.TEanno.gff3"
    local species_log="${dir}/edta.log"

    # Skip if already complete
    if [ -f "$final_output" ]; then
        log "[$subdir] Already complete (found $final_output). Skipping."
        return 0
    fi

    log "[$subdir] === Starting EDTA ==="
    log "[$subdir]   genome: $genome"
    log "[$subdir]   cds:    $cds"

    cd "$dir" || { log "[$subdir] ERROR: cannot cd to $dir"; return 1; }

    absl_fixed=0

    for attempt in $(seq 1 $MAX_RETRIES); do
        log "[$subdir] Attempt $attempt/$MAX_RETRIES"

        EDTA.pl \
            --genome "$genome" \
            --cds "$cds" \
            --anno 1 \
            --sensitive 1 \
            --threads $THREADS \
            >> "$species_log" 2>&1
        edta_exit=$?

        # Success check
        if [ -f "$final_output" ]; then
            log "[$subdir] EDTA completed successfully on attempt $attempt."
            return 0
        fi

        log "[$subdir] EDTA exited with code $edta_exit; checking for known errors..."

        # --- Bug 1: absl-py / TIR-Learner TypeError ---
        if grep -qE "TypeError|absl" "$species_log" 2>/dev/null && [ $absl_fixed -eq 0 ]; then
            log "[$subdir] Detected TIR-Learner/absl-py error. Applying fix..."
            fix_absl_py
            absl_fixed=1
            log "[$subdir] absl-py fixed. Retrying EDTA (will resume from checkpoint)..."
            continue
        fi

        # --- Bug 2: LTR_retriever *.mod.mod.* naming ---
        if find "${dir}/${prefix}.fa.mod.EDTA.raw" -name "*.mod.mod.*" 2>/dev/null | grep -q .; then
            log "[$subdir] Detected LTR_retriever naming bug. Applying fix..."
            fix_ltr_retriever_naming "$dir" "$prefix"
            log "[$subdir] LTR naming fixed. Retrying EDTA (will resume from checkpoint)..."
            continue
        fi

        # --- Unknown error ---
        if [ $attempt -lt $MAX_RETRIES ]; then
            log "[$subdir] Unknown error. Last 20 lines of log:"
            tail -20 "$species_log" | while read line; do log "[$subdir]   $line"; done
            log "[$subdir] Retrying..."
        else
            log "[$subdir] FAILED after $MAX_RETRIES attempts. Check $species_log"
            return 1
        fi
    done

    log "[$subdir] FAILED: max retries exhausted."
    return 1
}

# -------------------------------------------------------
log "======================================================="
log "EDTA liverwort batch run starting"
log "EDTA binary: $EDTA_BIN"
log "Threads per run: $THREADS"
log "Species: ${#SPECIES[@]}"
log "======================================================="

overall_exit=0
for entry in "${SPECIES[@]}"; do
    subdir=$(echo "$entry" | awk '{print $1}')
    prefix=$(echo "$entry" | awk '{print $2}')
    run_edta "$subdir" "$prefix" || overall_exit=1
    log ""  # blank line between species
done

log "======================================================="
log "Batch run complete. Exit status: $overall_exit"
log "======================================================="
exit $overall_exit
