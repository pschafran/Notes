#!/bin/bash
# check_odp_files.sh
# Validates protein/chrom/genome file consistency for ODP bryophyte species.
# For each species: fixes stop codons, checks protein<->chrom name matching,
# duplicates, and chrom chromosome names against genome FASTA sequence names.
#
# Stop codon handling:
#   trailing '.' (legitimate stop codon) -> truncated, sequence kept
#   internal '.' (in-frame stop codon)   -> sequence removed from protein + chrom
# Original files are backed up to .bak before any modification.
#
# Usage: bash check_odp_files.sh [BASE_DIR]
# Default BASE_DIR: /media/data/projects/hornwort_sex_chromosomes/analysis/synteny/odp_bryophytes_20260318

BASE=${1:-/media/data/projects/hornwort_sex_chromosomes/analysis/synteny/odp_bryophytes_20260318}

SPECIES="Cepur Encon Ensed Hycur Phafr Phpat Ptkno Spfal Sycan Syrur Mapol Maqua Ricav Riflu Rinat"

for sp in $SPECIES; do
  prot=$BASE/proteins/$sp.fa
  chrom=$BASE/chrom/$sp.chrom
  genome=$BASE/genomes/$sp.fa

  # Check files exist
  for f in "$prot" "$chrom" "$genome"; do
    if [ ! -f "$f" ]; then
      echo "$sp: MISSING FILE: $f"
      continue 2
    fi
  done

  # Fix stop codons in protein sequences.
  # Python prints summary to stderr; removed sequence names to stdout.
  removed_names=$(python3 - "$prot" 2>/tmp/${sp}_stop_report.txt <<'PYEOF'
import sys, shutil

prot_file = sys.argv[1]

# Parse FASTA, preserving original line structure
records = []
with open(prot_file) as f:
    header = None
    seq_lines = []
    for line in f:
        line = line.rstrip('\r\n')
        if line.startswith('>'):
            if header is not None:
                records.append((header, seq_lines))
            header = line
            seq_lines = []
        else:
            seq_lines.append(line)
    if header is not None:
        records.append((header, seq_lines))

trailing_seqs, internal_seqs = [], []
trailing_total, internal_total = 0, 0
keep_records = []
remove_names = []

for header, seq_lines in records:
    name = header[1:].split()[0]
    new_lines = list(seq_lines)

    # Find last non-empty sequence line
    last_idx = next((i for i in range(len(new_lines)-1, -1, -1) if new_lines[i]), None)

    # Strip trailing stop first
    had_trailing = False
    if last_idx is not None and new_lines[last_idx].endswith('.'):
        new_lines[last_idx] = new_lines[last_idx][:-1]
        had_trailing = True

    # Count internal stops in remaining sequence
    n_internal = sum(line.count('.') for line in new_lines)

    if n_internal > 0:
        # Internal stops -> remove sequence entirely
        internal_total += n_internal
        internal_seqs.append(name)
        remove_names.append(name)
    else:
        # No internal stops -> keep (possibly with trailing stop truncated)
        if had_trailing:
            trailing_total += 1
            trailing_seqs.append(name)
        keep_records.append((header, new_lines))

# Report to stderr
if trailing_seqs:
    print(f"  TRAILING_STOP: {len(trailing_seqs)} seqs, {trailing_total} stops truncated", file=sys.stderr)
    for s in trailing_seqs[:5]: print(f"    {s}", file=sys.stderr)
if internal_seqs:
    print(f"  INTERNAL_STOP: {len(internal_seqs)} seqs, {internal_total} IFS — sequences removed", file=sys.stderr)
    for s in internal_seqs[:5]: print(f"    {s}", file=sys.stderr)

# Backup and rewrite protein file if any changes
if trailing_seqs or internal_seqs:
    shutil.copy2(prot_file, prot_file + '.bak')
    with open(prot_file, 'w') as f:
        for header, seq_lines in keep_records:
            f.write(header + '\n')
            for line in seq_lines:
                f.write(line + '\n')

# Print removed names to stdout for bash to filter chrom
for name in remove_names:
    print(name)
PYEOF
)

  stop_report=$(cat /tmp/${sp}_stop_report.txt)
  rm -f /tmp/${sp}_stop_report.txt

  if [ -n "$stop_report" ]; then
    echo "$sp: STOP_CODONS fixed:" >&2
    echo "$stop_report" >&2
  fi

  # If sequences were removed, backup and filter the chrom file
  if [ -n "$removed_names" ]; then
    cp "$chrom" "$chrom.bak"
    grep -vFwf <(echo "$removed_names") "$chrom.bak" > "$chrom"
    echo "WARNING: $sp — protein and chrom files updated; originals saved to .bak" >&2
  elif [ -n "$stop_report" ]; then
    echo "WARNING: $sp — protein file updated (trailing stops truncated); original saved to .bak" >&2
  fi

  # Consistency checks (run on updated files)
  prot_names=$(grep "^>" $prot | awk '{print $1}' | sed 's/>//' | sort)
  chrom_genes=$(awk '{print $1}' $chrom | sort)
  chrom_chroms=$(awk '{print $2}' $chrom | sort -u)
  genome_seqs=$(grep "^>" $genome | awk '{print $1}' | sed 's/>//' | sort -u)

  np=$(echo "$prot_names" | wc -l)
  nc=$(echo "$chrom_genes" | wc -l)

  dup_prot=$(echo "$prot_names" | uniq -d)
  dup_chrom=$(echo "$chrom_genes" | uniq -d)
  only_prot=$(comm -23 <(echo "$prot_names") <(echo "$chrom_genes"))
  only_chrom=$(comm -13 <(echo "$prot_names") <(echo "$chrom_genes"))
  missing_chr=$(comm -23 <(echo "$chrom_chroms") <(echo "$genome_seqs"))

  errors=""
  [ -n "$dup_prot" ]    && errors+="  DUP_PROT: $(echo "$dup_prot" | wc -l) duplicates\n"
  [ -n "$dup_chrom" ]   && errors+="  DUP_CHROM: $(echo "$dup_chrom" | wc -l) duplicates\n"
  [ -n "$only_prot" ]   && errors+="  PROT_NOT_IN_CHROM: $(echo "$only_prot" | wc -l) entries\n$(echo "$only_prot" | head -5 | sed 's/^/    /')\n"
  [ -n "$only_chrom" ]  && errors+="  CHROM_NOT_IN_PROT: $(echo "$only_chrom" | wc -l) entries\n$(echo "$only_chrom" | head -5 | sed 's/^/    /')\n"
  [ -n "$missing_chr" ] && errors+="  CHR_NOT_IN_GENOME: $(echo "$missing_chr" | wc -l) chromosomes\n$(echo "$missing_chr" | head -5 | sed 's/^/    /')\n"

  if [ -z "$errors" ]; then
    echo "$sp: OK ($np proteins, $nc chrom entries)"
  else
    echo "$sp: ERRORS ($np proteins, $nc chrom entries)"
    printf "$errors"
  fi
done
