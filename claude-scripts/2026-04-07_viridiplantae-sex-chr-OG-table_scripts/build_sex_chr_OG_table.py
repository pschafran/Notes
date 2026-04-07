#! /usr/bin/env python3

# Builds sex_chr_near_universal_OGs.tsv — annotated table of near-universal
# hornwort sex chromosome orthogroups with model organism orthologs.
#
# Input:  OrthoFinder/Results_Mar30_1/Orthogroups/Orthogroups.tsv
# Output: OrthoFinder/Results_Mar30_1/Orthogroups/sex_chr_near_universal_OGs.tsv
#
# Orthogroups covered (17 total):
#   All U+V:        OG0000002, OG0000031
#   All U, no V:    OG0002583, OG0008579
#   All but 1 sample: OG0000179 (PhphyV), OG0000824 (LedusV)
#   All but Anang:  OG0000460, OG0000985, OG0001152, OG0003154, OG0003206
#   All but Ledus:  OG0000020, OG0000064, OG0000617, OG0000660, OG0002735
#   All but Phphy:  OG0001301

import csv

OGS_TSV = "OrthoFinder/Results_Mar30_1/Orthogroups/Orthogroups.tsv"
OUTPATH = "OrthoFinder/Results_Mar30_1/Orthogroups/sex_chr_near_universal_OGs.tsv"

# --- Load orthogroups ---
ogs = {}
with open(OGS_TSV) as f:
    reader = csv.reader(f, delimiter='\t')
    header = next(reader)
    at_col = next(i for i, h in enumerate(header) if 'Arabidopsis_thaliana' in h)
    mp_col = next(i for i, h in enumerate(header) if 'Marchantia_polymorpha' in h)
    pp_col = next(i for i, h in enumerate(header) if 'Physcomitrium_patens' in h)
    print(f"Column indices — At: {at_col}, Mp: {mp_col}, Pp: {pp_col}")
    for row in reader:
        og = row[0]
        def parse(col):
            if col >= len(row) or not row[col].strip():
                return []
            return [g.strip().replace('A_thaliana_', '') for g in row[col].split(',') if g.strip()]
        ogs[og] = {'at': parse(at_col), 'mp': parse(mp_col), 'pp': parse(pp_col)}

def fmt(genes, max_n=10, priority_prefixes=None):
    """Format gene list: list all if <=max_n, else first max_n with MpUg/MpVg prioritised."""
    if not genes:
        return "—"
    if priority_prefixes and len(genes) > max_n:
        priority = [g for g in genes if any(g.startswith(p) for p in priority_prefixes)]
        others   = [g for g in genes if not any(g.startswith(p) for p in priority_prefixes)]
        selected = priority + others
    else:
        selected = genes
    total = len(genes)
    shown = selected[:max_n]
    result = ", ".join(shown)
    if total > max_n:
        result += f" (+{total - max_n} more)"
    return result

# --- Table rows: (Category, OG, Missing, Preferred name, eggNOG description, Web-searched function) ---
data = [
    ("All U+V",         "OG0000002", "—",           "—",       "transposition, RNA-mediated",
     "Retroelement-derived (ATMG mitochondrial ORFs); likely transposon-related"),
    ("All U+V",         "OG0000031", "—",           "—",       "—",
     "Unknown"),
    ("All U, no V",     "OG0002583", "All V absent","FGMYB",   "transcription factor",
     "R2R3-MYB TF (FGMYB); promotes female/U gametophyte identity in Marchantia"),
    ("All U, no V",     "OG0008579", "All V absent","—",       "—",
     "Unknown; MpUg00040 is U-chromosome-linked in Marchantia"),
    ("All but Phphy V", "OG0000179", "PhphyV",      "HK3",     "His Kinase A (phospho-acceptor) domain",
     "HK3 — cytokinin receptor histidine kinase"),
    ("All but Ledus V", "OG0000824", "LedusV",      "TPL",     "Topless-related protein",
     "TPL/TPR — Topless transcriptional co-repressor"),
    ("All but Anang",   "OG0000460", "Anang U+V",   "HDG2",    "homeobox-leucine zipper protein",
     "HDG2 — Homeodomain GLABROUS 2, HD-ZIP IV transcription factor"),
    ("All but Anang",   "OG0000985", "Anang U+V",   "TCP23",   "Transcription factor",
     "TCP23 — class II TCP transcription factor"),
    ("All but Anang",   "OG0001152", "Anang U+V",   "LF4",     "Protein tyrosine kinase",
     "MAK/MOK-type RCK kinase (Chlamydomonas LF4 ortholog); implicated in cilia/flagella length control and male germ cell function"),
    ("All but Anang",   "OG0003154", "Anang U+V",   "EDR1",    "serine/threonine-protein kinase",
     "EDR1 — Enhanced Disease Resistance 1, MAPKKK"),
    ("All but Anang",   "OG0003206", "Anang U+V",   "DSK2a/b", "Ubiquitin family",
     "DSK2a/b — UBL-UBA ubiquitin shuttle receptor"),
    ("All but Ledus",   "OG0000020", "Ledus U+V",   "xbaIM",   "N(4)/N(6)-methyltransferase family",
     "Adenine N6-methyltransferase (restriction-modification system origin)"),
    ("All but Ledus",   "OG0000064", "Ledus U+V",   "—",       "—",
     "Unknown"),
    ("All but Ledus",   "OG0000617", "Ledus U+V",   "UBC31",   "ubiquitin-conjugating enzyme family",
     "UBC31 — E2 ubiquitin-conjugating enzyme"),
    ("All but Ledus",   "OG0000660", "Ledus U+V",   "USP24",   "peptidase C19 family",
     "USP24 — ubiquitin-specific protease (deubiquitinase)"),
    ("All but Ledus",   "OG0002735", "Ledus U+V",   "—",       "DEAD box helicase family",
     "DEAD-box RNA helicase"),
    ("All but Phphy",   "OG0001301", "Phphy U+V",   "CDF6",    "Cyclic dof factor",
     "CDF6 — Cycling DOF Factor 6, circadian clock regulation of flowering time"),
]

# --- Write output ---
header_row = ["Category", "OG", "Missing", "Preferred name", "eggNOG description",
              "Web-searched function", "At orthologs", "Mp orthologs", "Pp orthologs"]

with open(OUTPATH, "w", newline="") as f:
    writer = csv.writer(f, delimiter="\t", quoting=csv.QUOTE_NONNUMERIC)
    writer.writerow(header_row)
    for row in data:
        cat, og, missing, pref, eggnogdesc, webfunc = row
        d = ogs[og]
        at_str = fmt(d['at'], max_n=10)
        mp_str = fmt(d['mp'], max_n=10, priority_prefixes=['MpUg', 'MpVg'])
        pp_str = fmt(d['pp'], max_n=10)
        writer.writerow([cat, og, missing, pref, eggnogdesc, webfunc, at_str, mp_str, pp_str])

print(f"Saved {OUTPATH} ({len(data)} rows)")
