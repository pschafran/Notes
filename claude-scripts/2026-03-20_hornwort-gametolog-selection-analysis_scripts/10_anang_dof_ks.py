#!/usr/bin/env python3
"""
Ka/Ks computation for the Anthoceros angustus DOF zinc-finger U-V paralog pair
(AnangRef2_chrU_g0091 vs AnangRef2_chrV_g0007).

Tests all 4 transcript combinations (2 U transcripts × 2 V transcripts):
  chrU_g0091.t1, chrU_g0091.t2  x  chrV_g0007.t1, chrV_g0007.t2

For each combination:
  1. blastp to get pident
  2. mafft protein alignment
  3. CDS threading onto protein alignment
  4. yn00 (PAML) to estimate Ka, Ks, omega

Outputs results to stdout and
  output/anang_dof_ks/anang_dof_results.tsv
"""

import subprocess, tempfile, os, re, itertools
import numpy as np

WDIR     = "/media/data/projects/hornwort_sex_chromosomes/analysis/synteny/circos_hornwort_sex_acc_chr/output/anang_dof_ks"
PROT_FA  = f"{WDIR}/anang_dof_prot_all.fa"
CDS_FA   = f"{WDIR}/anang_dof_cds_all.fa"

U_TRANSCRIPTS = ["chrU_g0091.t1", "chrU_g0091.t2"]
V_TRANSCRIPTS = ["chrV_g0007.t1", "chrV_g0007.t2"]

# ── Helpers ────────────────────────────────────────────────────────────────────

def read_fasta(path):
    seqs = {}
    name = None
    with open(path) as f:
        for line in f:
            line = line.rstrip()
            if line.startswith(">"):
                name = line[1:].split()[0]
                seqs[name] = []
            elif name:
                seqs[name].append(line)
    return {k: "".join(v) for k, v in seqs.items()}


def write_fasta(path, seqs):
    with open(path, "w") as f:
        for name, seq in seqs.items():
            f.write(f">{name}\n{seq}\n")


def run_blastp(prot_dict, name_u, name_v, workd):
    """Return pident for a single U vs V protein pair."""
    q = os.path.join(workd, "query.fa")
    db = os.path.join(workd, "db.fa")
    write_fasta(q,  {name_u: prot_dict[name_u]})
    write_fasta(db, {name_v: prot_dict[name_v]})
    r = subprocess.run(
        ["blastp", "-query", q, "-subject", db,
         "-outfmt", "6 pident qlen slen length", "-max_hsps", "1"],
        capture_output=True, text=True)
    if r.stdout.strip():
        parts = r.stdout.strip().split("\t")
        return float(parts[0])
    return np.nan


def run_mafft(prot_dict, name_u, name_v, workd):
    """Align two proteins with mafft; return dict {name: aln_seq}."""
    inp = os.path.join(workd, "mafft_in.fa")
    out = os.path.join(workd, "mafft_out.fa")
    write_fasta(inp, {name_u: prot_dict[name_u], name_v: prot_dict[name_v]})
    subprocess.run(
        ["mafft", "--quiet", "--auto", inp],
        stdout=open(out, "w"), stderr=subprocess.DEVNULL, check=True)
    return read_fasta(out)


def thread_cds_onto_aln(aln_seq, cds_seq):
    """
    Given a gapped protein alignment and a CDS, produce a gapped codon alignment.
    Returns None if CDS length is inconsistent.
    """
    cds_clean = cds_seq.replace("-", "")
    # Strip trailing stop if present
    if len(cds_clean) % 3 == 0:
        last_codon = cds_clean[-3:].upper()
        if last_codon in ("TAA", "TAG", "TGA"):
            cds_clean = cds_clean[:-3]
    n_codons_expected = sum(1 for c in aln_seq if c != "-")
    n_codons_cds = len(cds_clean) // 3
    if n_codons_cds != n_codons_expected:
        return None
    result = []
    codon_i = 0
    for aa in aln_seq:
        if aa == "-":
            result.append("---")
        else:
            result.append(cds_clean[codon_i*3 : codon_i*3+3])
            codon_i += 1
    return "".join(result)


def build_codon_aln(aln, cds_dict, name_u, name_v):
    """Thread both sequences; remove gap columns; return (codon_aln_u, codon_aln_v) or None."""
    ca_u = thread_cds_onto_aln(aln[name_u], cds_dict[name_u])
    ca_v = thread_cds_onto_aln(aln[name_v], cds_dict[name_v])
    if ca_u is None or ca_v is None:
        return None
    # Remove columns where either side is a gap codon
    codons_u = [ca_u[i:i+3] for i in range(0, len(ca_u), 3)]
    codons_v = [ca_v[i:i+3] for i in range(0, len(ca_v), 3)]
    keep_u, keep_v = [], []
    for cu, cv in zip(codons_u, codons_v):
        if cu != "---" and cv != "---":
            keep_u.append(cu)
            keep_v.append(cv)
    return "".join(keep_u), "".join(keep_v)


def run_yn00(seq_u, seq_v, name_u, name_v, workd):
    """Write yn00 input/control, run yn00, parse Ks/Ka/omega."""
    aln_path = os.path.join(workd, "yn00.aln")
    ctl_path = os.path.join(workd, "yn00.ctl")
    out_path = os.path.join(workd, "yn00.out")
    n = len(seq_u) // 3
    with open(aln_path, "w") as f:
        f.write(f" 2 {len(seq_u)}\n")
        f.write(f"{name_u[:30]:<30s}  {seq_u}\n")
        f.write(f"{name_v[:30]:<30s}  {seq_v}\n")
    with open(ctl_path, "w") as f:
        f.write(f"seqfile  = {aln_path}\n")
        f.write(f"outfile  = {out_path}\n")
        f.write("verbose  = 0\n")
    subprocess.run(["yn00", ctl_path], capture_output=True, cwd=workd)
    return parse_yn00(out_path)


def parse_yn00(out_path):
    """Parse YN section from yn00 output; return (Ks, Ka, omega) or (nan, nan, nan)."""
    try:
        content = open(out_path).read()
    except FileNotFoundError:
        return np.nan, np.nan, np.nan
    in_yn = False
    for line in content.splitlines():
        if "Yang & Nielsen" in line:
            in_yn = True
            continue
        if in_yn:
            nums = re.findall(r"[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?", line)
            if len(nums) >= 9:
                try:
                    # Columns: seq1 seq2 S N t kappa omega dN±SE dS±SE
                    # indices:  0    1   2 3 4 5     6     7  8   9  10
                    omega = float(nums[6])
                    Ka    = float(nums[7])
                    Ks    = float(nums[9])
                    return Ks, Ka, omega
                except (ValueError, IndexError):
                    pass
    return np.nan, np.nan, np.nan


# ── Main ───────────────────────────────────────────────────────────────────────

prot_dict = read_fasta(PROT_FA)
cds_dict  = read_fasta(CDS_FA)

# Print sequence lengths
print("Transcript lengths (aa / nt):")
for t in U_TRANSCRIPTS + V_TRANSCRIPTS:
    aa = len(prot_dict.get(t, ""))
    nt = len(cds_dict.get(t, ""))
    print(f"  {t}: {aa} aa, {nt} nt")

print()
print(f"{'U transcript':<22} {'V transcript':<22} {'pident':>7} {'Ks':>7} {'Ka':>7} {'omega':>7} {'n_codons':>9}")
print("-" * 85)

results = []
for name_u, name_v in itertools.product(U_TRANSCRIPTS, V_TRANSCRIPTS):
    with tempfile.TemporaryDirectory() as workd:
        pident = run_blastp(prot_dict, name_u, name_v, workd)
        aln    = run_mafft(prot_dict, name_u, name_v, workd)
        ca     = build_codon_aln(aln, cds_dict, name_u, name_v)
        if ca is None:
            print(f"  {name_u:<22} {name_v:<22}  CDS/protein length mismatch — skipped")
            continue
        seq_u, seq_v = ca
        n_codons = len(seq_u) // 3
        Ks, Ka, omega = run_yn00(seq_u, seq_v, name_u, name_v, workd)
        results.append(dict(U=name_u, V=name_v, pident=pident,
                            Ks=Ks, Ka=Ka, omega=omega, n_codons=n_codons))
        print(f"  {name_u:<22} {name_v:<22} {pident:>7.1f} {Ks:>7.3f} {Ka:>7.3f} {omega:>7.3f} {n_codons:>9}")

# Write TSV
out_tsv = f"{WDIR}/anang_dof_results.tsv"
with open(out_tsv, "w") as f:
    f.write("U_transcript\tV_transcript\tpident\tKs\tKa\tomega\tn_codons\n")
    for r in results:
        f.write(f"{r['U']}\t{r['V']}\t{r['pident']:.2f}\t{r['Ks']:.4f}\t{r['Ka']:.4f}\t{r['omega']:.4f}\t{r['n_codons']}\n")
print(f"\nSaved: {out_tsv}")

# Identify best combination (lowest Ks among valid)
valid = [r for r in results if not np.isnan(r["Ks"]) and r["Ks"] > 0]
if valid:
    best = min(valid, key=lambda r: r["Ks"])
    print(f"\nBest combination: {best['U']} x {best['V']}")
    print(f"  pident={best['pident']:.1f}%  Ks={best['Ks']:.4f}  Ka={best['Ka']:.4f}  omega={best['omega']:.4f}")
    # OrthoFinder used primary transcripts: chrU_g0091.t1 and chrV_g0007.t2
    primary = next((r for r in results if r["U"] == "chrU_g0091.t1" and r["V"] == "chrV_g0007.t2"), None)
    if primary:
        print(f"\nOrthoFinder primary (t1 x t2): Ks={primary['Ks']:.4f}  Ka={primary['Ka']:.4f}  omega={primary['omega']:.4f}")
        if best["U"] != "chrU_g0091.t1" or best["V"] != "chrV_g0007.t2":
            delta = best["Ks"] - primary["Ks"]
            print(f"Improvement vs primary: ΔKs={delta:+.4f}")
