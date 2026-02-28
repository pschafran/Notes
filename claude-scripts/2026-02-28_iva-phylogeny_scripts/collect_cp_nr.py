#!/usr/bin/env python3
"""
collect_cp_nr.py — Collect and orient cp/nr sequences for phylogeny paper samples

CP orientation convention:
  1. Select isoform where ndhD (SSC marker) is forward (sstart < send via tblastn)
  2. If rbcL (LSC marker) is also forward in that isoform → LSC is inverted → revcomp whole seq

NR orientation convention:
  18S rRNA on + strand (detected by barrnap --kingdom euk)

Output:
  phylogeny_paper/cp/{sample_id}_cp.fasta   — one seq per file, header >{sample_id}
  phylogeny_paper/nr/{sample_id}_nr.fasta   — one seq per file, header >{sample_id}
  phylogeny_paper/cp/cp_manifest.tsv
  phylogeny_paper/nr/nr_manifest.tsv
  phylogeny_paper/cp/all_cp.fasta           — combined
  phylogeny_paper/nr/all_nr.fasta           — combined
"""

import os
import sys
import glob
import subprocess
import tempfile
from pathlib import Path

BASE    = Path("/media/data/projects/iva_phylogeny/analysis")
OUTCP   = BASE / "phylogeny_paper" / "cp"
OUTNR   = BASE / "phylogeny_paper" / "nr"
NDHDREF = str(BASE / "ndhD.faa")
RBCLREF = str(BASE / "rbcL.faa")

# ─── Sample definitions ────────────────────────────────────────────────────
# (sample_id, sample_dir, cp_subdir, nr_subdir)
SAMPLES = [
    # EREMID batch
    ("angustifolia_46", BASE / "46_EREMID_240811",                    "embplant_pt", "embplant_nr"),
    ("annua_01",        BASE / "annua/01_EREMID_240811",              "embplant_pt", "embplant_nr"),
    ("cultigen_133",    BASE / "cultigen/13-3-33_EREMID_240811",      "embplant_pt", "embplant_nr"),
    ("frutescens_55",   BASE / "frutescens/55_EREMID_240811",         "embplant_pt", "embplant_nr"),
    ("microcephala_57", BASE / "microcephala/57_EREMID_240811",       "embplant_pt", "embplant_nr"),
    ("texensis_48",     BASE / "texensis/48_EREMID_240811",           "embplant_pt", "embplant_nr"),
    # STURM batch – species-prefixed dirs
    ("angustifolia_77", BASE / "angustifolia/77_STURM_251201adq30ft", "angustifolia_LA_77_pt",  "angustifolia_LA_77_nr"),
    ("asperifolia_60",  BASE / "asperifolia/60_STURM_251201adq30ft",  "asperifolia_FL_60_pt",   "asperifolia_FL_60_nr"),
    ("hayesiana_92",    BASE / "hayesiana/92_STURM_251201adq30ft",    "hayesiana_CA1_92_pt",    "hayesiana_CA1_92_nr"),
    ("hayesiana_93",    BASE / "hayesiana/93_STURM_251201adq30ft",    "hayesiana_CA2_93_pt",    "hayesiana_CA2_93_nr"),
    ("xanthifolia_64",  BASE / "xanthifolia/64_STURM_251201adq30ft",  "xanthiifolia_CO1_64_pt", "xanthiifolia_CO1_64_nr"),
    # STURM batch – embplant dirs
    ("axillaris_89",    BASE / "axillaris/89_STURM_251201adq30ft",    "embplant_pt", "embplant_nr"),
    ("frutescens_36",   BASE / "frutescens/36_STURM_251201adq30ft",   "embplant_pt", "embplant_nr"),
    ("frutescens_49",   BASE / "frutescens/49_STURM_251201adq30ft",   "embplant_pt", "embplant_nr"),
    ("frutescens_58",   BASE / "frutescens/58_STURM_251201adq30ft",   "embplant_pt", "embplant_nr"),
    ("frutescens_63",   BASE / "frutescens/63_STURM_251201adq30ft",   "embplant_pt", "embplant_nr"),
    ("frutescens_68",   BASE / "frutescens/68_STURM_251201adq30ft",   "embplant_pt", "embplant_nr"),
    ("frutescens_71",   BASE / "frutescens/71_STURM_251201adq30ft",   "embplant_pt", "embplant_nr"),
    ("frutescens_86",   BASE / "frutescens/86_STURM_251201adq30ft",   "embplant_pt", "embplant_nr"),
    ("frutescens_94",   BASE / "frutescens/94_STURM_251201adq30ft",   "embplant_pt", "embplant_nr"),
    ("imbricata_67",    BASE / "imbricata/67_STURM_251201adq30ft",    "embplant_pt", "embplant_nr"),
    ("imbricata_72",    BASE / "imbricata/72_STURM_251201adq30ft",    "embplant_pt", "embplant_nr"),
    ("imbricata_85",    BASE / "imbricata/85_STURM_251201adq30ft",    "embplant_pt", "embplant_nr"),
]

# ─── Sequence helpers ──────────────────────────────────────────────────────
COMP = str.maketrans("ACGTacgtNnRrYyKkMmSsWwBbDdHhVv",
                     "TGCAtgcaNnYyRrMmKkSsWwVvHhDdBb")

def revcomp(seq):
    return seq.translate(COMP)[::-1]

def read_fasta(path):
    seqs = []
    with open(path) as f:
        header, chunks = None, []
        for line in f:
            line = line.rstrip()
            if line.startswith(">"):
                if header is not None:
                    seqs.append((header, "".join(chunks)))
                header, chunks = line[1:], []
            else:
                chunks.append(line)
        if header is not None:
            seqs.append((header, "".join(chunks)))
    return seqs

def write_fasta(path, records, wrap=80):
    with open(path, "w") as f:
        for header, seq in records:
            f.write(f">{header}\n")
            for i in range(0, len(seq), wrap):
                f.write(seq[i:i+wrap] + "\n")

# ─── BLAST helper ──────────────────────────────────────────────────────────
def tblastn_strand(query_faa, seq):
    """
    tblastn query_faa against seq.
    Returns True  (forward, sstart < send),
            False (reverse, sstart > send),
            None  (no hit at e<=1e-10).
    """
    with tempfile.NamedTemporaryFile(mode="w", suffix=".fasta", delete=False) as tmp:
        tmp.write(f">subject\n{seq}\n")
        tmp_path = tmp.name
    try:
        r = subprocess.run(
            ["tblastn", "-query", query_faa, "-subject", tmp_path,
             "-outfmt", "6 sstart send bitscore",
             "-max_hsps", "1", "-evalue", "1e-10"],
            capture_output=True, text=True
        )
        best_score, best_fwd = 0.0, None
        for line in r.stdout.splitlines():
            parts = line.split("\t")
            if len(parts) < 3:
                continue
            sstart, send, score = int(parts[0]), int(parts[1]), float(parts[2])
            if score > best_score:
                best_score = score
                best_fwd = (sstart < send)
        return best_fwd
    finally:
        os.unlink(tmp_path)

# ─── CP collection ─────────────────────────────────────────────────────────
def collect_cp(sample_id, pt_dir):
    """
    Returns (oriented_seq, source_basename, notes) or (None, None, error).
    """
    # Collect candidates – skip .revcomp. files created by earlier manual work
    candidates = sorted(
        f for f in glob.glob(str(pt_dir / "*.fasta"))
        if ".revcomp." not in f
    )
    if not candidates:
        return None, None, "NO_FASTA"

    # imbricata_67 has 10 paths across 3 repeat_patterns;
    # restrict to repeat_pattern1 (simplest structure)
    if sample_id == "imbricata_67":
        rp1 = [f for f in candidates if "repeat_pattern1" in f]
        if rp1:
            candidates = rp1

    # ── Step 1: find isoform where ndhD is forward ──
    chosen, seq, notes = None, None, ""
    for f in candidates:
        recs = read_fasta(f)
        if not recs:
            continue
        s = recs[0][1]
        if tblastn_strand(NDHDREF, s) is True:
            chosen, seq = f, s
            notes = "ndhD_fwd"
            break

    if chosen is None:
        # Try reverse complement of each candidate
        for f in candidates:
            recs = read_fasta(f)
            if not recs:
                continue
            s_rc = revcomp(recs[0][1])
            if tblastn_strand(NDHDREF, s_rc) is True:
                chosen, seq = f, s_rc
                notes = "ndhD_fwd_after_revcomp"
                break

    if chosen is None:
        # Give up: use first candidate as-is
        recs = read_fasta(candidates[0])
        chosen = candidates[0]
        seq = recs[0][1] if recs else ""
        notes = "ndhD_not_detected;first_candidate_as_is"

    # ── Step 2: check rbcL — if also forward, LSC is inverted → revcomp ──
    rbcL = tblastn_strand(RBCLREF, seq)
    if rbcL is True:
        seq = revcomp(seq)
        notes += ";rbcL_fwd→LSC_inverted→revcomp"
    elif rbcL is False:
        notes += ";rbcL_rev(ok)"
    else:
        notes += ";rbcL_not_detected"

    return seq, os.path.basename(chosen), notes

# ─── NR collection ─────────────────────────────────────────────────────────
def pick_nr_file(sample_id, nr_dir):
    """
    Pick a representative nr fasta file.
    Returns (filepath, selection_note) or (None, skip_reason).
    """
    if sample_id == "hayesiana_93":
        return None, "SKIPPED:1000_nr_paths_unresolvable"

    all_f = sorted(glob.glob(str(nr_dir / "*.fasta")))
    if not all_f:
        return None, "NO_FASTA"

    if len(all_f) == 1:
        return all_f[0], "single_file"

    # microcephala_57: graph1 is only 63 bp; graph2 is the real assembly
    if sample_id == "microcephala_57":
        g2 = [f for f in all_f if "graph2" in f]
        if g2:
            return g2[0], "graph2_chosen(graph1=63bp)"

    # Prefer a non-repeat-pattern file if available
    non_rp = [f for f in all_f if "repeat_pattern" not in f]
    if non_rp:
        return non_rp[0], f"first_non-rp(of {len(all_f)} files)"

    # Otherwise use repeat_pattern1
    rp1 = [f for f in all_f if "repeat_pattern1" in f]
    if rp1:
        return rp1[0], f"repeat_pattern1(of {len(all_f)} files)"

    return all_f[0], f"first_of_{len(all_f)}_files"


def barrnap_18S_strand(seq):
    """Run barrnap on seq; return '+', '-', or None."""
    with tempfile.NamedTemporaryFile(mode="w", suffix=".fasta", delete=False) as tmp:
        tmp.write(f">seq\n{seq}\n")
        tmp_path = tmp.name
    try:
        r = subprocess.run(
            ["barrnap", "--kingdom", "euk", "--quiet", tmp_path],
            capture_output=True, text=True
        )
        for line in r.stdout.splitlines():
            if line.startswith("#") or not line.strip():
                continue
            parts = line.split("\t")
            if len(parts) >= 9 and "18S" in parts[8]:
                return parts[6]   # '+' or '-'
        return None
    finally:
        os.unlink(tmp_path)


def collect_nr(sample_id, nr_dir):
    """
    Returns (oriented_seq, source_basename, notes) or (None, None, error).
    """
    chosen, pick_note = pick_nr_file(sample_id, nr_dir)
    if chosen is None:
        return None, None, pick_note

    recs = read_fasta(chosen)
    if not recs:
        return None, os.path.basename(chosen), "EMPTY_FILE"

    seq = recs[0][1]
    notes = pick_note

    strand = barrnap_18S_strand(seq)
    if strand == "+":
        notes += ";18S_fwd(ok)"
    elif strand == "-":
        seq = revcomp(seq)
        notes += ";18S_rev→revcomp"
    else:
        notes += ";18S_not_detected"

    return seq, os.path.basename(chosen), notes

# ─── Main ──────────────────────────────────────────────────────────────────
def main():
    OUTCP.mkdir(exist_ok=True)
    OUTNR.mkdir(exist_ok=True)

    cp_rows = [("sample_id", "seq_len_bp", "source_file", "notes")]
    nr_rows = [("sample_id", "seq_len_bp", "source_file", "notes")]
    cp_combined, nr_combined = [], []

    for sample_id, sample_dir, cp_subdir, nr_subdir in SAMPLES:
        pt_dir  = Path(sample_dir) / cp_subdir
        nr_dir  = Path(sample_dir) / nr_subdir

        print(f"\n{'='*62}\n  {sample_id}\n{'='*62}")

        # ── CP ──────────────────────────────────────────────────────
        if not pt_dir.exists():
            print(f"  CP  MISSING DIR: {pt_dir}")
            cp_rows.append((sample_id, "NA", "MISSING_DIR", "NO_CP_DIR"))
        else:
            seq, src, notes = collect_cp(sample_id, pt_dir)
            if seq is None:
                print(f"  CP  FAILED — {notes}")
                cp_rows.append((sample_id, "NA", str(src or ""), notes))
            else:
                out = OUTCP / f"{sample_id}_cp.fasta"
                write_fasta(out, [(sample_id, seq)])
                print(f"  CP  {len(seq):,} bp  →  {out.name}  [{notes}]")
                cp_rows.append((sample_id, str(len(seq)), src, notes))
                cp_combined.append((sample_id, seq))

        # ── NR ──────────────────────────────────────────────────────
        if not nr_dir.exists():
            print(f"  NR  MISSING DIR: {nr_dir}")
            nr_rows.append((sample_id, "NA", "MISSING_DIR", "NO_NR_DIR"))
        else:
            seq, src, notes = collect_nr(sample_id, nr_dir)
            if seq is None:
                print(f"  NR  FAILED/SKIPPED — {notes}")
                nr_rows.append((sample_id, "NA", str(src or ""), notes))
            else:
                out = OUTNR / f"{sample_id}_nr.fasta"
                write_fasta(out, [(sample_id, seq)])
                print(f"  NR  {len(seq):,} bp  →  {out.name}  [{notes}]")
                nr_rows.append((sample_id, str(len(seq)), src, notes))
                nr_combined.append((sample_id, seq))

    # Combined files
    write_fasta(OUTCP / "all_cp.fasta", cp_combined)
    write_fasta(OUTNR / "all_nr.fasta", nr_combined)

    # Manifests
    for path, rows in [(OUTCP / "cp_manifest.tsv", cp_rows),
                       (OUTNR / "nr_manifest.tsv", nr_rows)]:
        with open(path, "w") as f:
            for row in rows:
                f.write("\t".join(row) + "\n")

    print(f"\n{'='*62}")
    print(f"  Done.")
    print(f"  CP: {len(cp_combined)} sequences  →  {OUTCP}")
    print(f"  NR: {len(nr_combined)} sequences  →  {OUTNR}")
    print(f"{'='*62}")


if __name__ == "__main__":
    main()
