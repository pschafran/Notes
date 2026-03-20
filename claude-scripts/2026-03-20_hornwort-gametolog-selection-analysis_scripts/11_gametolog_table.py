#!/usr/bin/env python3
"""
Comprehensive table of all sex-chromosome gametolog Ka/Ks estimates.

Includes all pairs from the four focal species with paired U/V genomes,
plus the Anthoceros angustus DOF zinc-finger pair computed separately.
Rows that did not pass filters are included with an exclusion reason.
Transcript substitutions made to improve Ka/Ks estimates (script 09) are noted.

Output: output/gametolog_kaks_table.tsv
"""

import csv, re
import numpy as np
import pandas as pd
from collections import defaultdict

# ── Paths ──────────────────────────────────────────────────────────────────────
OG_TSV      = ("/media/data/projects/hornwort_sex_chromosomes/analysis/orthofinder"
               "/hornworts_20260302/OrthoFinder/Results_Mar05/Orthogroups/Orthogroups.tsv")
OG_MERGE_MAP = ("/media/data/projects/hornwort_sex_chromosomes/analysis/sex_chr_gene_content"
                "/merged_ogs/og_merge_map.tsv")
SPLIT_OG_TSV = ("/media/data/projects/hornwort_sex_chromosomes/analysis/sex_chr_gene_content"
                "/split_og_candidates.tsv")
EGGNOG_DIR  = "/media/data/projects/hornwort_sex_chromosomes/analysis/functional_annotations"
TRANSCRIPT_CHECK_TSV = ("/media/data/projects/hornwort_sex_chromosomes/analysis/synteny"
                        "/circos_hornwort_sex_acc_chr/output/transcript_check_results.tsv")
OUTDIR      = ("/media/data/projects/hornwort_sex_chromosomes/analysis/synteny"
               "/circos_hornwort_sex_acc_chr/output")
ANANG_GAMETOLOG_TSV = ("/media/data/projects/hornwort_sex_chromosomes/analysis/synteny"
                       "/circos_hornwort_sex_acc_chr/output/anang_dof_ks/anang_gametolog_results.tsv")

KS_MIN  = 0.05
KS_MAX  = 3.0

# ── Species ────────────────────────────────────────────────────────────────────
SPECIES = [
    dict(label="L. dussii",
         ks_tsv=("/media/data/projects/hornwort_sex_chromosomes/analysis"
                 "/Leiosporoceros/sex_chromosome_analyses/rbh_ks/LedusF_LedusM_RBH_Ks.tsv"),
         f_sex="LedusF.S3", m_sex="LedusM.S5",
         og_col_f="LedusF_PROT_primary", og_col_m="LedusM_PROT_primary",
         egg_f=f"{EGGNOG_DIR}/LedusF_PROT.emapper.annotations",
         egg_m=f"{EGGNOG_DIR}/LedusM_PROT.emapper.annotations"),
    dict(label="P. proskaueri",
         ks_tsv=("/media/data/projects/hornwort_sex_chromosomes/analysis"
                 "/Paraphymatoceros_proskaueri/sex_chromosome_analyses/rbh_ks/PaproF_PaproM_RBH_Ks.tsv"),
         f_sex="PaproF.S5", m_sex="PaproM.S5",
         og_col_f="PaproF_PROT_primary", og_col_m="PaproM_PROT_primary",
         egg_f=f"{EGGNOG_DIR}/PaproF_PROT.emapper.annotations",
         egg_m=f"{EGGNOG_DIR}/PaproM_PROT.emapper.annotations"),
    dict(label="P. phymatodes",
         ks_tsv=("/media/data/projects/hornwort_sex_chromosomes/analysis"
                 "/Phymatoceros_phymatodes/sex_chromosome_analyses/rbh_ks/PhphyF_PhphyM_RBH_Ks.tsv"),
         f_sex="PhphyF.S5", m_sex="PhphyM.S5",
         og_col_f="PhphyF_PROT_primary", og_col_m="PhphyM_PROT_primary",
         egg_f=f"{EGGNOG_DIR}/PhphyF_PROT.emapper.annotations",
         egg_m=f"{EGGNOG_DIR}/PhphyM_PROT.emapper.annotations"),
    dict(label="P. chiloensis",
         ks_tsv=("/media/data/projects/hornwort_sex_chromosomes/analysis"
                 "/Phaeomegaceros_fimbriatus/sex_chromosome_analyses/rbh_ks/PhchiF_PhchiM_RBH_Ks.tsv"),
         f_sex="PhchiF.S4", m_sex="PhchiM.S6",
         og_col_f="PhchiF_PROT_primary", og_col_m="PhchiM_PROT_primary",
         egg_f=f"{EGGNOG_DIR}/PhchiF_PROT.emapper.annotations",
         egg_m=f"{EGGNOG_DIR}/PhchiM_PROT.emapper.annotations"),
]

# ── Helper ─────────────────────────────────────────────────────────────────────
def gene_base(transcript):
    parts = transcript.split(".")
    return ".".join(parts[:2]) if len(parts) >= 2 else transcript

# ── Load OG maps ───────────────────────────────────────────────────────────────
print("Loading orthogroups...")
gene_to_og = {}
all_og_cols = {sp[c] for sp in SPECIES for c in ("og_col_f", "og_col_m")}
with open(OG_TSV) as fh:
    reader = csv.reader(fh, delimiter="\t")
    header = next(reader)
    col_idx = {c: i for i, c in enumerate(header) if c in all_og_cols}
    for row in reader:
        og = row[0]
        for col, idx in col_idx.items():
            if idx < len(row) and row[idx].strip():
                for gene in row[idx].split(", "):
                    gene = gene.strip()
                    if gene:
                        gene_to_og[gene] = og

og_merge = {}
with open(OG_MERGE_MAP) as fh:
    next(fh)
    for line in fh:
        parts = line.rstrip("\n").split("\t")
        og_merge[parts[0]] = parts[1]

def resolve_og(og):
    return og_merge.get(og, og)

uv_split_canonical = {}
uv_split_desc = {}
with open(SPLIT_OG_TSV) as fh:
    header = next(fh).rstrip("\n").split("\t")
    pt_idx = header.index("pair_type")
    la_idx = header.index("label_a")
    lb_idx = header.index("label_b")
    for line in fh:
        parts = line.rstrip("\n").split("\t")
        if parts[pt_idx] != "UV_paralog":
            continue
        og_a, og_b = resolve_og(parts[0]), resolve_og(parts[1])
        canon = "+".join(sorted([og_a, og_b]))
        uv_split_canonical[og_a] = canon
        uv_split_canonical[og_b] = canon
        la = parts[la_idx] if la_idx < len(parts) else "-"
        lb = parts[lb_idx] if lb_idx < len(parts) else "-"
        uv_split_desc[canon] = la if la not in ("-", "") else lb

# ── Load eggNOG ────────────────────────────────────────────────────────────────
print("Loading eggNOG annotations...")
gene_to_desc = {}
for sp in SPECIES:
    for egg_path in [sp["egg_f"], sp["egg_m"]]:
        with open(egg_path) as fh:
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

# ── Load transcript override data (script 09) ──────────────────────────────────
print("Loading transcript override data...")
tc = pd.read_csv(TRANSCRIPT_CHECK_TSV, sep="\t")
# key: (species, gene_F_base, gene_M_base)
override = {}
for _, row in tc.iterrows():
    override[(row["species"], row["gene_F"], row["gene_M"])] = row

# ── Resolve OG for a gene pair ─────────────────────────────────────────────────
def resolve_pair_og(gene_f, gene_m):
    og_f_raw = gene_to_og.get(gene_f)
    og_m_raw = gene_to_og.get(gene_m)
    og_f = resolve_og(og_f_raw) if og_f_raw else None
    og_m = resolve_og(og_m_raw) if og_m_raw else None
    if og_f and og_m and og_f != og_m:
        cf = uv_split_canonical.get(og_f)
        cm = uv_split_canonical.get(og_m)
        if cf and cf == cm:
            return cf
        return og_f
    return og_f or og_m

# ── Collect all rows ───────────────────────────────────────────────────────────
print("Collecting pairs...")
rows = []

for sp in SPECIES:
    df = pd.read_csv(sp["ks_tsv"], sep="\t")
    df["Ka_Ks"] = pd.to_numeric(df["Ka_Ks"], errors="coerce")
    df["Ka"]    = pd.to_numeric(df["Ka"],    errors="coerce")
    df["Ks"]    = pd.to_numeric(df["Ks"],    errors="coerce")
    sex = df[(df["scaffold_F"] == sp["f_sex"]) & (df["scaffold_M"] == sp["m_sex"])].copy()

    for _, row in sex.iterrows():
        gf = row["gene_F"]   # transcript name
        gm = row["gene_M"]
        gf_base = gene_base(gf)
        gm_base = gene_base(gm)

        og = resolve_pair_og(gf, gm)
        desc = gene_to_desc.get(gf) or gene_to_desc.get(gm) or "-"
        if desc == "-" and og in uv_split_desc:
            desc = uv_split_desc[og]

        primary_Ks    = row["Ks"]
        primary_Ka    = row["Ka"]
        primary_omega = row["Ka_Ks"]
        if "aln_length_codons" in row.index and pd.notna(row["aln_length_codons"]):
            aln_len = int(row["aln_length_codons"])
        elif "n_codons" in row.index and pd.notna(row["n_codons"]):
            aln_len = int(row["n_codons"])
        else:
            aln_len = pd.NA
        blastp_pident = (row["blastp_pident"] if "blastp_pident" in row.index else pd.NA)

        # Check for transcript override
        ov = override.get((sp["label"], gf_base, gm_base))
        if ov is not None and ov["flag"] == "IMPROVED":
            used_F    = ov["best_F"]
            used_M    = ov["best_M"]
            used_Ks   = ov["best_Ks"]
            used_Ka   = ov["best_Ka"]
            used_omega = ov["best_omega"]
            if pd.notna(ov["best_pident"]):
                blastp_pident = ov["best_pident"]
            primary_F = ov["primary_F"]
            primary_M = ov["primary_M"]
            transcript_note = (f"Transcripts replaced: OrthoFinder used "
                               f"{primary_F}/{primary_M} (Ks={ov['primary_Ks']:.3f}); "
                               f"best alignment uses {used_F}/{used_M} "
                               f"(Ks={ov['best_Ks']:.3f})")
        else:
            used_F     = gf
            used_M     = gm
            used_Ks    = primary_Ks
            used_Ka    = primary_Ka
            used_omega = primary_omega
            transcript_note = ""

        # Filter status (based on best/used values)
        if pd.isna(used_Ks) or pd.isna(used_omega) or used_omega <= 0:
            status = "excluded: omega missing"
        elif used_Ks < KS_MIN:
            status = "excluded: Ks < 0.05 (too recent)"
        elif used_Ks > KS_MAX:
            status = f"excluded: Ks > {KS_MAX} (saturated)"
        else:
            status = "passed_Ks_filter"   # n_species check applied later

        rows.append(dict(
            species         = sp["label"],
            gene_F          = gf_base,
            gene_M          = gm_base,
            orthogroup      = og or "-",
            transcript_F    = used_F,
            transcript_M    = used_M,
            Ks              = round(used_Ks, 4) if pd.notna(used_Ks) else "",
            Ka              = round(used_Ka, 4) if pd.notna(used_Ka) else "",
            omega           = round(used_omega, 4) if pd.notna(used_omega) else "",
            aln_length_codons = aln_len,
            blastp_pident   = (round(blastp_pident, 1) if pd.notna(blastp_pident) else ""),
            filter_status   = status,
            transcript_note = transcript_note,
            description     = desc,
        ))

# ── A. angustus pairs (scripts 10 + 12) ───────────────────────────────────────
# ODP RBH pairs (script 12)
anang_df = pd.read_csv(ANANG_GAMETOLOG_TSV, sep="\t")
for _, ar in anang_df.iterrows():
    og = str(ar["orthogroup"]) if pd.notna(ar.get("orthogroup")) else "-"
    desc = str(ar["description"]) if pd.notna(ar.get("description")) else "-"
    ks_val = ar["Ks"]
    if pd.isna(ks_val) or ks_val < KS_MIN:
        status = "excluded: Ks < 0.05 (too recent)" if (not pd.isna(ks_val) and ks_val < KS_MIN) else "excluded: omega missing"
    elif ks_val > KS_MAX:
        status = f"excluded: Ks > {KS_MAX} (saturated)"
    else:
        status = "passed_Ks_filter"
    rows.append(dict(
        species           = "A. angustus",
        gene_F            = str(ar["gene_U"]),
        gene_M            = str(ar["gene_V"]),
        orthogroup        = og,
        transcript_F      = str(ar["transcript_U"]),
        transcript_M      = str(ar["transcript_V"]),
        Ks                = round(float(ks_val), 4) if pd.notna(ks_val) else "",
        Ka                = round(float(ar["Ka"]), 4) if pd.notna(ar.get("Ka")) else "",
        omega             = round(float(ar["omega"]), 4) if pd.notna(ar.get("omega")) else "",
        aln_length_codons = int(ar["n_codons"]) if pd.notna(ar.get("n_codons")) else pd.NA,
        blastp_pident     = round(float(ar["blastp_pident"]), 1) if pd.notna(ar.get("blastp_pident")) else "",
        filter_status     = status,
        transcript_note   = "Ka/Ks computed from ODP RBH pair (intra-genomic U-V paralog, no Ks TSV)",
        description       = desc,
    ))

# DOF pair from script 10 (OrthoFinder-identified, not in ODP RBH)
rows.append(dict(
    species           = "A. angustus",
    gene_F            = "AnangRef2_chrU_g0091",
    gene_M            = "AnangRef2_chrV_g0007",
    orthogroup        = "OG0011295+OG0011420",
    transcript_F      = "chrU_g0091.t1",
    transcript_M      = "chrV_g0007.t1",
    Ks                = 2.5727,
    Ka                = 0.6995,
    omega             = 0.2719,
    aln_length_codons = 415,
    blastp_pident     = 38.1,
    filter_status     = "passed_Ks_filter",
    transcript_note   = ("OrthoFinder primary was chrU_g0091.t1/chrV_g0007.t2; "
                         "best alignment uses chrU_g0091.t1/chrV_g0007.t1 (Ks=2.811→2.573); "
                         "Ka/Ks computed via script 10 (OrthoFinder UV_paralog, not in ODP RBH)"),
    description       = "Cyclic dof factor",
))

# ── Determine n_species per OG (for passed rows) and update filter status ──────
og_passed_species = defaultdict(set)
for r in rows:
    if r["filter_status"] == "passed_Ks_filter" and r["orthogroup"] != "-":
        og_passed_species[r["orthogroup"]].add(r["species"])

for r in rows:
    if r["filter_status"] == "passed_Ks_filter":
        og = r["orthogroup"]
        n_sp = len(og_passed_species.get(og, set()))
        if n_sp < 2:
            r["filter_status"] = "excluded: n_species < 2"
        else:
            r["filter_status"] = "included"

# ── Write output ───────────────────────────────────────────────────────────────
out_tsv = f"{OUTDIR}/gametolog_kaks_table.tsv"
cols = ["species", "gene_F", "gene_M", "orthogroup",
        "transcript_F", "transcript_M",
        "Ks", "Ka", "omega",
        "aln_length_codons", "blastp_pident",
        "filter_status", "transcript_note", "description"]

df_out = pd.DataFrame(rows, columns=cols)

# Sort: by species order, then by Ks
sp_order = {sp["label"]: i for i, sp in enumerate(SPECIES)}
sp_order["A. angustus"] = len(SPECIES)
df_out["_sp_ord"] = df_out["species"].map(sp_order)
df_out = df_out.sort_values(["_sp_ord", "Ks"]).drop(columns="_sp_ord")

df_out.to_csv(out_tsv, sep="\t", index=False)
print(f"\nSaved: {out_tsv}")
print(f"Total rows: {len(df_out)}")
print(f"  included:           {(df_out.filter_status=='included').sum()}")
print(f"  excluded (Ks>3.0):  {df_out.filter_status.str.startswith('excluded: Ks >').sum()}")
print(f"  excluded (Ks<0.05): {df_out.filter_status.str.startswith('excluded: Ks <').sum()}")
print(f"  excluded (n_sp<2):  {(df_out.filter_status=='excluded: n_species < 2').sum()}")
print(f"  excluded (omega):   {(df_out.filter_status=='excluded: omega missing').sum()}")
print(f"\nRows with transcript substitution notes: "
      f"{(df_out.transcript_note != '').sum()}")
