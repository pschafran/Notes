#!/usr/bin/env python3
"""
Ka/Ks computation for all Anthoceros angustus U-V gametolog pairs identified
by ODP RBH analysis (odp_hornworts_20260310, AnangF vs AnangM).

For each of the 14 RBH pairs, tests all transcript combinations (from
AnangRef2_PROT.fa), selects the best alignment (lowest Ks), and computes
Ka/Ks via mafft + yn00.

Also maps each pair to an OrthoFinder orthogroup (using the AnangRef2 primary
transcript column).

Output:
  output/anang_dof_ks/anang_gametolog_results.tsv
"""

import subprocess, tempfile, os, re, itertools, csv
import numpy as np
import pandas as pd
from collections import defaultdict

# ── Paths ──────────────────────────────────────────────────────────────────────
RBH_FILE  = ("/media/data/projects/hornwort_sex_chromosomes/analysis/synteny"
             "/odp_hornworts_20260310/odp/step2-figures/synteny_nocolor"
             "/AnangF_AnangM_xy_reciprocal_best_hits.plotted.rbh")
PROT_FA   = ("/media/data/projects/hornwort_sex_chromosomes/analysis"
             "/Anthoceros_angustus/Ref2/AnangRef2_PROT.fa")
CDS_FA    = ("/media/data/projects/hornwort_sex_chromosomes/analysis/synteny"
             "/circos_hornwort_sex_acc_chr/output/anang_dof_ks/anang_all_cds.fa")
PROT_ALL_FA = ("/media/data/projects/hornwort_sex_chromosomes/analysis/synteny"
               "/circos_hornwort_sex_acc_chr/output/anang_dof_ks/anang_all_prot.fa")
OG_TSV    = ("/media/data/projects/hornwort_sex_chromosomes/analysis/orthofinder"
             "/hornworts_20260302/OrthoFinder/Results_Mar05/Orthogroups/Orthogroups.tsv")
OG_MERGE  = ("/media/data/projects/hornwort_sex_chromosomes/analysis/sex_chr_gene_content"
             "/merged_ogs/og_merge_map.tsv")
SPLIT_OG  = ("/media/data/projects/hornwort_sex_chromosomes/analysis/sex_chr_gene_content"
             "/split_og_candidates.tsv")
EGGNOG    = ("/media/data/projects/hornwort_sex_chromosomes/analysis/functional_annotations"
             "/AnangRef2_PROT.emapper.annotations")
OUTDIR    = ("/media/data/projects/hornwort_sex_chromosomes/analysis/synteny"
             "/circos_hornwort_sex_acc_chr/output/anang_dof_ks")

KS_MAX = 3.0
KS_MIN = 0.05

# ── Load RBH pairs (chrU vs chrV only) ────────────────────────────────────────
pairs = []   # list of (gene_U, gene_V)  — base names with AnangRef2_ prefix
with open(RBH_FILE) as f:
    header = next(f).rstrip("\n").split("\t")
    f_gene_col  = header.index("AnangF_gene")
    m_gene_col  = header.index("AnangM_gene")
    f_scaf_col  = header.index("AnangF_scaf")
    m_scaf_col  = header.index("AnangM_scaf")
    for line in f:
        parts = line.rstrip("\n").split("\t")
        fg, mg = parts[f_gene_col], parts[m_gene_col]
        fs, ms = parts[f_scaf_col], parts[m_scaf_col]
        if (fs == "chrU" and ms == "chrV") or (fs == "chrV" and ms == "chrU"):
            # Normalise so U is always first
            if fs == "chrU":
                gene_u = "AnangRef2_" + fg.split("AnangRef2_")[-1] if fg.startswith("AnangRef2_") else "AnangRef2_" + fg
                gene_v = "AnangRef2_" + mg.split("AnangRef2_")[-1] if mg.startswith("AnangRef2_") else "AnangRef2_" + mg
            else:
                gene_u = "AnangRef2_" + mg.split("AnangRef2_")[-1] if mg.startswith("AnangRef2_") else "AnangRef2_" + mg
                gene_v = "AnangRef2_" + fg.split("AnangRef2_")[-1] if fg.startswith("AnangRef2_") else "AnangRef2_" + fg
            # Strip transcript suffix to get gene base
            def base(t): return re.sub(r'\.t\d+$', '', t)
            pairs.append((base(gene_u), base(gene_v)))

print(f"Loaded {len(pairs)} chrU-chrV RBH pairs")

# ── Load all sequences ─────────────────────────────────────────────────────────
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

prot_all  = read_fasta(PROT_ALL_FA)
cds_all   = read_fasta(CDS_FA)

# gffread strips AnangRef2_ prefix → keys are like "chrU_g0022.t1"
# AnangRef2_PROT.fa keys are like "AnangRef2_chrU_g0022.t1"
# Build a unified lookup that handles both
def seq_lookup(d, name):
    """Try name, then name with/without AnangRef2_ prefix."""
    if name in d:
        return d[name]
    short = name.replace("AnangRef2_", "")
    if short in d:
        return d[short]
    return None

# List transcripts for a gene base name
def get_transcripts(d, gene_base):
    """Return all transcript names matching gene_base (with or without prefix)."""
    short_base = gene_base.replace("AnangRef2_", "")
    ts = [k for k in d if k == gene_base or
          k.startswith(gene_base + ".") or
          k == short_base or
          k.startswith(short_base + ".")]
    return sorted(set(ts))

# ── OrthoFinder OG lookup ─────────────────────────────────────────────────────
print("Loading orthogroups...")
gene_to_og = {}
with open(OG_TSV) as fh:
    reader = csv.reader(fh, delimiter="\t")
    header = next(reader)
    # Find Anang primary transcripts column
    anang_col = next((i for i, c in enumerate(header)
                      if c == "AnangRef2_PROT_primary_transcripts"), None)
    if anang_col is None:
        print("  WARNING: AnangRef2 column not found in OrthoFinder TSV")
    else:
        for row in reader:
            og = row[0]
            if anang_col < len(row) and row[anang_col].strip():
                for gene in row[anang_col].split(", "):
                    gene = gene.strip()
                    if gene:
                        gene_to_og[gene] = og

og_merge = {}
with open(OG_MERGE) as fh:
    next(fh)
    for line in fh:
        parts = line.rstrip("\n").split("\t")
        og_merge[parts[0]] = parts[1]

def resolve_og(og):
    return og_merge.get(og, og)

uv_split_canonical = {}
uv_split_desc = {}
with open(SPLIT_OG) as fh:
    hdr = next(fh).rstrip("\n").split("\t")
    pt_i = hdr.index("pair_type")
    la_i = hdr.index("label_a")
    lb_i = hdr.index("label_b")
    for line in fh:
        parts = line.rstrip("\n").split("\t")
        if parts[pt_i] != "UV_paralog":
            continue
        og_a, og_b = resolve_og(parts[0]), resolve_og(parts[1])
        canon = "+".join(sorted([og_a, og_b]))
        uv_split_canonical[og_a] = canon
        uv_split_canonical[og_b] = canon
        la = parts[la_i] if la_i < len(parts) else "-"
        lb = parts[lb_i] if lb_i < len(parts) else "-"
        uv_split_desc[canon] = la if la not in ("-", "") else lb

def pair_og(gene_u, gene_v):
    """Resolve orthogroup for a U-V pair, handling UV_paralog splits."""
    # Try with .t1 suffix (primary transcript used in OrthoFinder)
    og_u_raw = gene_to_og.get(gene_u + ".t1") or gene_to_og.get(gene_u + ".t2")
    og_v_raw = gene_to_og.get(gene_v + ".t1") or gene_to_og.get(gene_v + ".t2")
    og_u = resolve_og(og_u_raw) if og_u_raw else None
    og_v = resolve_og(og_v_raw) if og_v_raw else None
    if og_u and og_v and og_u != og_v:
        cf = uv_split_canonical.get(og_u)
        cm = uv_split_canonical.get(og_v)
        if cf and cf == cm:
            return cf
        return og_u
    return og_u or og_v

# ── eggNOG annotation ─────────────────────────────────────────────────────────
gene_to_desc = {}
try:
    with open(EGGNOG) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 21:
                continue
            gene = parts[0]
            desc = parts[7].strip()
            pref = parts[8].strip()
            label = pref if pref and pref != "-" else desc
            if label and label != "-":
                gene_to_desc[gene] = label
except FileNotFoundError:
    print(f"  WARNING: eggNOG file not found: {EGGNOG}")

# ── yn00 helpers (same as script 10) ─────────────────────────────────────────
def write_fasta(path, seqs):
    with open(path, "w") as f:
        for name, seq in seqs.items():
            f.write(f">{name}\n{seq}\n")

def run_blastp(seq_u, seq_v, name_u, name_v, workd):
    q  = os.path.join(workd, "q.fa")
    db = os.path.join(workd, "db.fa")
    write_fasta(q,  {name_u: seq_u})
    write_fasta(db, {name_v: seq_v})
    r = subprocess.run(["blastp", "-query", q, "-subject", db,
                        "-outfmt", "6 pident", "-max_hsps", "1"],
                       capture_output=True, text=True)
    return float(r.stdout.strip().split()[0]) if r.stdout.strip() else np.nan

def run_mafft(seq_u, seq_v, name_u, name_v, workd):
    inp = os.path.join(workd, "in.fa")
    out = os.path.join(workd, "out.fa")
    write_fasta(inp, {name_u: seq_u, name_v: seq_v})
    subprocess.run(["mafft", "--quiet", "--auto", inp],
                   stdout=open(out, "w"), stderr=subprocess.DEVNULL, check=True)
    return read_fasta(out)

def thread_cds(aln_seq, cds_seq):
    cds_clean = cds_seq.replace("-", "")
    if len(cds_clean) % 3 == 0:
        last = cds_clean[-3:].upper()
        if last in ("TAA", "TAG", "TGA"):
            cds_clean = cds_clean[:-3]
    n_exp = sum(1 for c in aln_seq if c != "-")
    if len(cds_clean) // 3 != n_exp:
        return None
    result = []
    ci = 0
    for aa in aln_seq:
        if aa == "-":
            result.append("---")
        else:
            result.append(cds_clean[ci*3:ci*3+3])
            ci += 1
    return "".join(result)

def build_codon_aln(aln, cds_u, cds_v, name_u, name_v):
    ca_u = thread_cds(aln[name_u], cds_u)
    ca_v = thread_cds(aln[name_v], cds_v)
    if ca_u is None or ca_v is None:
        return None
    codons_u = [ca_u[i:i+3] for i in range(0, len(ca_u), 3)]
    codons_v = [ca_v[i:i+3] for i in range(0, len(ca_v), 3)]
    ku, kv = [], []
    for cu, cv in zip(codons_u, codons_v):
        if cu != "---" and cv != "---":
            ku.append(cu)
            kv.append(cv)
    return "".join(ku), "".join(kv)

def run_yn00(seq_u, seq_v, name_u, name_v, workd):
    aln = os.path.join(workd, "yn.aln")
    ctl = os.path.join(workd, "yn.ctl")
    out = os.path.join(workd, "yn.out")
    with open(aln, "w") as f:
        f.write(f" 2 {len(seq_u)}\n")
        f.write(f"{name_u[:30]:<30s}  {seq_u}\n")
        f.write(f"{name_v[:30]:<30s}  {seq_v}\n")
    with open(ctl, "w") as f:
        f.write(f"seqfile = {aln}\noutfile = {out}\nverbose = 0\n")
    subprocess.run(["yn00", ctl], capture_output=True, cwd=workd)
    try:
        content = open(out).read()
    except FileNotFoundError:
        return np.nan, np.nan, np.nan
    in_yn = False
    for line in content.splitlines():
        if "Yang & Nielsen" in line:
            in_yn = True
            continue
        if in_yn:
            nums = re.findall(r"[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?", line)
            if len(nums) >= 10:
                try:
                    return float(nums[9]), float(nums[7]), float(nums[6])  # Ks, Ka, omega
                except (ValueError, IndexError):
                    pass
    return np.nan, np.nan, np.nan

# ── Main loop ──────────────────────────────────────────────────────────────────
results = []
print(f"\n{'Pair':<45} {'RBH transcripts':<36} {'Best combo':<36} "
      f"{'pident':>7} {'Ks':>7} {'Ka':>7} {'omega':>7} {'n_cod':>6}")
print("-" * 155)

for gene_u, gene_v in pairs:
    ts_u = get_transcripts(prot_all, gene_u)
    ts_v = get_transcripts(prot_all, gene_v)

    if not ts_u or not ts_v:
        print(f"  {gene_u} x {gene_v}: sequences not found — skipped")
        continue

    best = None
    combos_tried = 0
    for tu, tv in itertools.product(ts_u, ts_v):
        seq_u_prot = seq_lookup(prot_all, tu)
        seq_v_prot = seq_lookup(prot_all, tv)
        seq_u_cds  = seq_lookup(cds_all, tu)
        seq_v_cds  = seq_lookup(cds_all, tv)
        if any(s is None for s in [seq_u_prot, seq_v_prot, seq_u_cds, seq_v_cds]):
            continue
        combos_tried += 1
        with tempfile.TemporaryDirectory() as workd:
            pident = run_blastp(seq_u_prot, seq_v_prot, tu, tv, workd)
            aln    = run_mafft(seq_u_prot, seq_v_prot, tu, tv, workd)
            ca     = build_codon_aln(aln, seq_u_cds, seq_v_cds, tu, tv)
            if ca is None:
                continue
            su, sv = ca
            Ks, Ka, omega = run_yn00(su, sv, tu, tv, workd)
            n_cod = len(su) // 3
            if not np.isnan(Ks) and Ks > 0:
                if best is None or Ks < best["Ks"]:
                    best = dict(tu=tu, tv=tv, pident=pident,
                                Ks=Ks, Ka=Ka, omega=omega, n_cod=n_cod)

    og = pair_og(gene_u, gene_v)
    desc = "-"
    for suffix in [".t1", ".t2"]:
        for prefix in ["AnangRef2_", ""]:
            key = prefix + gene_u.replace("AnangRef2_", "") + suffix
            desc = gene_to_desc.get(key, desc)
            if desc != "-":
                break
        if desc != "-":
            break
    if desc == "-" and og in uv_split_desc:
        desc = uv_split_desc[og]

    if best:
        # Note if RBH transcript differs from best
        rbh_u = [t for t in ts_u if any(t.endswith(s) for s in [".t1",".t2",".t3"])][0] if ts_u else "?"
        rbh_v = [t for t in ts_v if any(t.endswith(s) for s in [".t1",".t2",".t3"])][0] if ts_v else "?"
        print(f"  {gene_u:<22} x {gene_v:<22}  "
              f"{rbh_u.split('_')[-1]+'|'+rbh_v.split('_')[-1]:<20}  "
              f"{best['tu'].split('_')[-1]+'|'+best['tv'].split('_')[-1]:<20}  "
              f"{best['pident']:>7.1f} {best['Ks']:>7.3f} {best['Ka']:>7.3f} "
              f"{best['omega']:>7.3f} {best['n_cod']:>6}")
        results.append(dict(
            species          = "A. angustus",
            gene_U           = gene_u,
            gene_V           = gene_v,
            orthogroup       = og or "-",
            transcript_U     = best["tu"],
            transcript_V     = best["tv"],
            n_combos         = combos_tried,
            blastp_pident    = round(best["pident"], 2),
            Ks               = round(best["Ks"], 4),
            Ka               = round(best["Ka"], 4),
            omega            = round(best["omega"], 4),
            n_codons         = best["n_cod"],
            description      = desc,
        ))
    else:
        print(f"  {gene_u:<22} x {gene_v:<22}  no valid Ka/Ks result")
        results.append(dict(
            species="A. angustus", gene_U=gene_u, gene_V=gene_v,
            orthogroup=og or "-", transcript_U="", transcript_V="",
            n_combos=combos_tried, blastp_pident=np.nan,
            Ks=np.nan, Ka=np.nan, omega=np.nan, n_codons=0,
            description=desc,
        ))

# ── Save ───────────────────────────────────────────────────────────────────────
out_tsv = f"{OUTDIR}/anang_gametolog_results.tsv"
df = pd.DataFrame(results)
df.to_csv(out_tsv, sep="\t", index=False)
print(f"\nSaved: {out_tsv}")

passed = df[df["Ks"].notna() & (df["Ks"] >= KS_MIN) & (df["Ks"] <= KS_MAX)]
print(f"\nPassed Ks filter ({KS_MIN}–{KS_MAX}): {len(passed)}/{len(df)}")
for _, r in passed.iterrows():
    flag = "purifying" if r["omega"] < 0.5 else ("diversifying" if r["omega"] > 1.0 else "neutral")
    print(f"  {r['gene_U']} x {r['gene_V']}  Ks={r['Ks']:.3f}  Ka={r['Ka']:.3f}  "
          f"omega={r['omega']:.3f}  {flag}  [{r['orthogroup']}]  {r['description'][:50]}")
