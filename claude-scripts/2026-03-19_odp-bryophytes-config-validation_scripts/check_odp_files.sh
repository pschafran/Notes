#!/bin/bash
# check_odp_files.sh
# Validates protein/chrom/genome file consistency for ODP bryophyte species.
# For each species: checks protein<->chrom name matching, duplicates, and
# chrom chromosome names against genome FASTA sequence names.
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
