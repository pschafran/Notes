#! /bin/bash
#
# Peter W. Schafran ps997@cornell.edu
#
# Usage: parse_eggnog.sh QUERY_SEQ_NAME[.txt] EGGNOG_COLUMN1 EGGNOG_COLUMN2...EGGNOG_COLUMNn
#
# Accepted EGGNOG_COLUMNS:
# SEED_ORTHOLOG, SEED_EVALUE, SEED_SCORE, BEST_TAX_LEVEL, PREFERRED_NAME, GO,
# EC, KEGG_KO, KEGG_PATHWAY, KEGG_MODULE, KEGG_REACTION, KEGG_RCLASS, BRITE,
# KEGG_TC, CAZY, BIGG, TAXONOMY, EGGNOG_OG, BEST_EGGNOG_OG, COG, FREE_TEXT
# or use ALL for whole line

QUERY=$1
PARAMS=()

# Determine input is file or string
if [ -f $1 ]
  then GREP_TYPE="-f"
else
  GREP_TYPE="-w"
fi

if [ $2 = "ALL" ]
  then grep -F $GREP_TYPE "$QUERY" /home/ps997/HornwortBase_20210503/Hornworts.emapper.annotations
else
  # Convert text command line args to corresponding column number in eggnog file
  for i in "${@:2}"
    do if [ $i = "SEED_ORTHOLOG" ]
      then PARAMS+=(2)
    elif [ $i = "SEED_EVALUE" ]
      then PARAMS+=(3)
    elif [ $i = "SEED_SCORE" ]
      then PARAMS+=(4)
    elif [ $i = "BEST_TAX_LEVEL" ]
      then PARAMS+=(5)
    elif [ $i = "PREFERRED_NAME" ]
      then PARAMS+=(6)
    elif [ $i = "GO" ]
      then PARAMS+=(7)
    elif [ $i = "EC" ]
      then PARAMS+=(8)
    elif [ $i = "KEGG_KO" ]
      then PARAMS+=(9)
    elif [ $i = "KEGG_PATHWAY" ]
      then PARAMS+=(10)
    elif [ $i = "KEGG_MODULE" ]
      then PARAMS+=(11)
    elif [ $i = "KEGG_REACTION" ]
      then PARAMS+=(12)
    elif [ $i = "KEGG_RCLASS" ]
      then PARAMS+=(13)
    elif [ $i = "BRITE" ]
      then PARAMS+=(14)
    elif [ $i = "KEGG_TC" ]
      then PARAMS+=(15)
    elif [ $i = "CAZY" ]
      then PARAMS+=(16)
    elif [ $i = "BIGG" ]
      then PARAMS+=(17)
    elif [ $i = "TAXONOMY" ]
      then PARAMS+=(18)
    elif [ $i = "EGGNOG_OG" ]
      then PARAMS+=(19)
    elif [ $i = "BEST_EGGNOG_OG" ]
      then PARAMS+=(20)
    elif [ $i = "COG" ]
      then PARAMS+=(21)
    elif [ $i = "FREE_TEXT" ]
      then PARAMS+=(22)
    fi
  done

  # Remove any double entries in input
  UNIQ_PARAMS=($(printf "%s\n" "${PARAMS[@]}" | sort | uniq))

  # Convert array to comma separated list to pass to cut
  DELIM=""
  FIELDS=""
  for i in "${UNIQ_PARAMS[@]}"
    do FIELDS="$FIELDS$DELIM$i"
    DELIM=","
  done

  grep -F $GREP_TYPE "$QUERY" /home/ps997/HornwortBase_20210503/Hornworts.emapper.annotations | cut -f 1,"$FIELDS"
fi
