#!/usr/bin/env python3
"""
Annotate new Isoetes LFY sequences for GenBank submission.
Uses AY541781.1 (I. engelmannii) as reference for exon/intron boundary.

Reference structure (AY541781.1, 1041 bp):
  exon 2:   <1..61   (partial at 5' end; encodes RFLEEVQHICRERGEKCPTK with codon_start=2)
  intron 2:  62..>1041 (GT..AG; partial at 3' end)
"""

import re
import sys
from Bio import SeqIO, pairwise2
from Bio.Seq import Seq

# ── Constants ────────────────────────────────────────────────────────────────
FASTA_FILE = "Kizzort_Branch_paper.fasta"
GB_FILE    = "I_engelmannii_AY541781.1.gb"
OUT_ANNO   = "new_sequences_annotation.tsv"
OUT_TBL    = "new_sequences.tbl"
OUT_FSA    = "new_sequences.fsa"
OUT_ISSUES = "annotation_issues.tsv"

REF_ID              = "Isoetes_engelmannii_AY541781.1"
REF_INTRON2_START   = 62          # 1-based in AY541781
REF_EXON2_END       = 61          # 1-based in AY541781
REF_CODON_START     = 2
REF_FIRST_CODON_POS = 2           # 1-based position in ref where first complete codon starts
REF_PROTEIN         = "RFLEEVQHICRERGEKCPTK"
# Expected C-terminal tail used for translation validation (last 5 aa before intron)
EXPECTED_CTAIL      = "CPTK"      # last 4 aa of ref protein (avoid stop codon edge effects)

ACCESSION_PAT = re.compile(r'[A-Z]{2}\d{6}\.\d')

# ── Load sequences ────────────────────────────────────────────────────────────
all_seqs = {}
with open(FASTA_FILE) as fh:
    for rec in SeqIO.parse(fh, "fasta"):
        all_seqs[rec.id] = rec

ref_seq = str(all_seqs[REF_ID].seq).upper()

new_seqs = {sid: rec for sid, rec in all_seqs.items()
            if not ACCESSION_PAT.search(sid)}

print(f"Reference: {REF_ID} ({len(ref_seq)} bp)")
print(f"New sequences to annotate: {len(new_seqs)}")

# ── Helper: pairwise alignment → map ref position to query ───────────────────
def map_ref_pos(ref, query, ref_pos, flank=30):
    """
    Global pairwise align ref[:ref_pos+flank] vs query[:ref_pos+flank+20].
    Returns the 1-based position in query that corresponds to ref_pos (1-based).
    Returns None if that reference position falls in a gap.
    """
    r = ref[:ref_pos + flank].replace('N', 'A')
    q = query[:ref_pos + flank + 20].replace('N', 'A')
    alns = pairwise2.align.globalms(r, q, 2, -1, -8, -0.5, one_alignment_only=True)
    if not alns:
        return None
    ref_aln, qry_aln = alns[0].seqA, alns[0].seqB
    rc = qc = 0
    for rb, qb in zip(ref_aln, qry_aln):
        if rb != '-':
            rc += 1
        if qb != '-':
            qc += 1
        if rc == ref_pos:
            return None if qb == '-' else qc
    return None

# ── Helper: find GT splice donor near an approximate position ─────────────────
def find_gt(seq, approx, window=8):
    """
    Search for 'GT' within ±window of approx (1-based).
    Returns 1-based position of the G in 'GT', or None.
    Prefers position closest to approx.
    """
    seq = seq.upper()
    candidates = []
    for offset in range(-window, window + 1):
        pos0 = approx - 1 + offset   # 0-based
        if 0 <= pos0 < len(seq) - 1 and seq[pos0:pos0+2] == 'GT':
            candidates.append((abs(offset), pos0 + 1))  # (distance, 1-based)
    if not candidates:
        return None
    candidates.sort()
    return candidates[0][1]

# ── Helper: verify CDS translation ───────────────────────────────────────────
def translate_cds(seq, exon_end, first_codon_start):
    """
    Translate from first_codon_start (1-based) through exon_end (1-based).
    Returns translated protein string.
    """
    coding = seq[first_codon_start - 1 : exon_end].upper()
    # Trim to multiple of 3
    coding = coding[:len(coding) - len(coding) % 3]
    return str(Seq(coding).translate())

def check_translation(seq, exon_end, codon_start):
    """
    Try all three codon_start values (1,2,3); return the one whose C-terminal
    tail best matches EXPECTED_CTAIL.  Returns (best_codon_start, protein, ok).
    """
    best = None
    for cs in [codon_start, 1, 2, 3]:  # prefer alignment-derived cs first
        first = cs  # first complete codon starts at position cs (1-based)
        protein = translate_cds(seq, exon_end, first)
        if protein.endswith(EXPECTED_CTAIL):
            return cs, protein, True
    # None matched: return alignment-derived codon_start with its protein
    protein = translate_cds(seq, exon_end, codon_start)
    return codon_start, protein, False

# ── Process each new sequence ─────────────────────────────────────────────────
results  = {}
issues   = {}

for sid in sorted(new_seqs):
    seq = str(new_seqs[sid].seq).upper()
    seq_len = len(seq)

    # 1. Map reference intron start (pos 62) → approximate position in this seq
    approx = map_ref_pos(ref_seq, seq, REF_INTRON2_START)
    if approx is None:
        issues[sid] = "Could not map ref intron-start position into sequence"
        continue

    # 2. Find GT splice donor
    intron_start = find_gt(seq, approx, window=8)
    if intron_start is None:
        issues[sid] = (f"No GT splice donor found within ±8 bp of mapped pos {approx}; "
                       f"context: ...{seq[max(0,approx-6):approx+6]}...")
        continue

    exon_end = intron_start - 1

    # 3. Derive codon_start from alignment:
    #    find where ref position 2 (first codon of ref CDS) falls in this seq
    first_codon_in_query = map_ref_pos(ref_seq, seq, REF_FIRST_CODON_POS)
    if first_codon_in_query is None or first_codon_in_query < 1:
        # Fallback: infer from intron position (assume same reading frame as ref)
        first_codon_in_query = ((intron_start - REF_INTRON2_START) % 3) + REF_CODON_START
    # Clamp to 1..3 (handles rare cases where alignment puts it at 4+)
    codon_start = ((first_codon_in_query - 1) % 3) + 1

    # 4. Verify / adjust translation; try other frames if C-tail doesn't match
    codon_start, protein, prot_ok = check_translation(seq, exon_end, codon_start)

    # 5. Confirm GT at splice site
    gt_bases = seq[intron_start - 1: intron_start + 1]

    results[sid] = {
        'seq_len':              seq_len,
        'approx_intron':        approx,
        'intron_start':         intron_start,
        'exon_end':             exon_end,
        'first_codon_in_query': first_codon_in_query,
        'codon_start':          codon_start,
        'gt_bases':             gt_bases,
        'protein':              protein,
        'prot_ok':              prot_ok,
    }

# ── Report results ────────────────────────────────────────────────────────────
ok_count   = sum(1 for r in results.values() if r['prot_ok'])
bad_prot   = {s: r for s, r in results.items() if not r['prot_ok']}
non_gt     = {s: r for s, r in results.items() if r['gt_bases'] != 'GT'}

print(f"\n{'='*60}")
print(f"Annotated:            {len(results)}")
print(f"  Correct translation: {ok_count}")
print(f"  Wrong translation:   {len(bad_prot)}")
print(f"  Non-GT splice site:  {len(non_gt)}")
print(f"Could not annotate:   {len(issues)}")
print(f"{'='*60}\n")

if bad_prot:
    print("PROTEIN MISMATCHES:")
    for sid, r in sorted(bad_prot.items()):
        print(f"  {sid}")
        print(f"    got:      {r['protein']}")
        print(f"    expected: {REF_PROTEIN}")
    print()

if non_gt:
    print("NON-GT SPLICE SITES:")
    for sid, r in sorted(non_gt.items()):
        print(f"  {sid}: splice site = '{r['gt_bases']}' at pos {r['intron_start']}")
    print()

if issues:
    print("COULD NOT ANNOTATE:")
    for sid, msg in sorted(issues.items()):
        print(f"  {sid}: {msg}")
    print()

# ── Write annotation TSV ──────────────────────────────────────────────────────
with open(OUT_ANNO, 'w') as fh:
    fh.write("seq_id\tseq_len\texon2_end\tintron2_start\tfirst_codon_in_query\t"
             "codon_start\tgt_bases\tprotein\tprot_ok\n")
    for sid, r in sorted(results.items()):
        fh.write(f"{sid}\t{r['seq_len']}\t{r['exon_end']}\t{r['intron_start']}\t"
                 f"{r['first_codon_in_query']}\t{r['codon_start']}\t"
                 f"{r['gt_bases']}\t{r['protein']}\t{r['prot_ok']}\n")
    for sid, msg in sorted(issues.items()):
        fh.write(f"{sid}\tERROR\t\t\t\t\t\t\t{msg}\n")

print(f"Annotation table → {OUT_ANNO}")

# ── Write feature table (.tbl) ────────────────────────────────────────────────
# tbl2asn feature table format:
#   >Feature SeqID
#   <start  end  feature_key
#           qualifier  value
# Partial: <N = extends before N; >N = extends beyond N

def tbl_start(pos, partial):
    return f"<{pos}" if partial else str(pos)

def tbl_end(pos, partial):
    return f">{pos}" if partial else str(pos)

INTRON_ONLY_SEQS = set()  # populated below

with open(OUT_TBL, 'w') as fh:
    # ── Sequences with exon + intron ──────────────────────────────────────────
    for sid, r in sorted(results.items()):
        ee  = r['exon_end']
        ist = r['intron_start']
        sl  = r['seq_len']
        cs  = r['codon_start']

        fh.write(f">Feature {sid}\n")

        # gene (partial both ends: exon + intron both truncated)
        fh.write(f"<1\t>{sl}\tgene\n")
        fh.write(f"\t\t\tgene\tLFY\n")

        # mRNA: only the exon 2 portion is present in this sequence
        fh.write(f"<1\t>{ee}\tmRNA\n")
        fh.write(f"\t\t\tproduct\tleafy\n")

        # CDS (partial at both ends)
        fh.write(f"<1\t>{ee}\tCDS\n")
        fh.write(f"\t\t\tcodon_start\t{cs}\n")
        fh.write(f"\t\t\tnote\tLFY\n")
        fh.write(f"\t\t\tproduct\tleafy\n")

        # exon 2 (partial at 5' end)
        fh.write(f"<1\t{ee}\texon\n")
        fh.write(f"\t\t\tnumber\t2\n")

        # intron 2 (partial at 3' end)
        fh.write(f"{ist}\t>{sl}\tintron\n")
        fh.write(f"\t\t\tnumber\t2\n")

    # ── Sequences that are entirely intronic (no exon 2 present) ─────────────
    # MS-08 and NC-05 start after the exon; annotate as partial gene + partial intron only
    for sid in sorted(issues):
        if sid not in new_seqs:
            continue
        sl = len(str(new_seqs[sid].seq))
        INTRON_ONLY_SEQS.add(sid)
        fh.write(f">Feature {sid}\n")
        fh.write(f"<1\t>{sl}\tgene\n")
        fh.write(f"\t\t\tgene\tLFY\n")
        fh.write(f"<1\t>{sl}\tintron\n")
        fh.write(f"\t\t\tnumber\t2\n")

print(f"Feature table       → {OUT_TBL}")

# ── Write new sequences FASTA (.fsa) ─────────────────────────────────────────
# Deflines: [organism=...] [mol_type=genomic DNA] — other source qualifiers
# use placeholders; user fills from metadata spreadsheet.
def parse_organism(sid):
    parts = sid.replace(';','').split('_')
    if parts[0] == 'Isoetes':
        ep = parts[1]
        if ep in ('sp', 'boomiiORmicrovela', 'Leary'):
            return "Isoetes sp."
        return f"Isoetes {ep}"
    elif parts[0] == 'Brunton':
        return "Isoetes sp."   # species unknown; user must fill
    elif parts[0] in ('MS', 'NC'):
        # e.g. MS-08_I_bigdawg  or NC-05_Isoetes_silvaticaE
        # Try to find "I_" or "Isoetes_" pattern
        for i, p in enumerate(parts):
            if p in ('I', 'Isoetes') and i + 1 < len(parts):
                ep = parts[i+1]
                return f"Isoetes {ep}" if p == 'Isoetes' else f"Isoetes sp."
        return "Isoetes sp."
    return "Isoetes sp."

all_tbl_seqs = set(results.keys()) | INTRON_ONLY_SEQS
with open(OUT_FSA, 'w') as fh:
    for sid in sorted(all_tbl_seqs):
        rec = new_seqs[sid]
        organism = parse_organism(sid)
        defline = (f">{sid} "
                   f"[organism={organism}] "
                   f"[mol_type=genomic DNA] "
                   f"[specimen-voucher=PLACEHOLDER] "
                   f"[country=PLACEHOLDER] "
                   f"[collection-date=PLACEHOLDER]")
        fh.write(defline + "\n")
        seq = str(rec.seq)
        for i in range(0, len(seq), 60):
            fh.write(seq[i:i+60] + "\n")

print(f"FASTA for submission → {OUT_FSA}")

# ── Issues file ───────────────────────────────────────────────────────────────
with open(OUT_ISSUES, 'w') as fh:
    fh.write("seq_id\tissue\n")
    for sid, msg in sorted(issues.items()):
        fh.write(f"{sid}\t{msg}\n")
    for sid, r in sorted(bad_prot.items()):
        fh.write(f"{sid}\tPROT_MISMATCH: got {r['protein']} expected {REF_PROTEIN}\n")
    for sid, r in sorted(non_gt.items()):
        if sid not in bad_prot:
            fh.write(f"{sid}\tNON_GT_SPLICE: {r['gt_bases']} at pos {r['intron_start']}\n")

print(f"Issues log          → {OUT_ISSUES}")
