#!/usr/bin/env python3
"""
Generate MCMCtree input files (PHYLIP + tree + control) for each HOG.
Option A: focal species only (no non-focal calibration taxa).

Calibrations from Penaloza-Bojaca et al. 2025 (secondary calibrations).
Root calibration (U/V divergence) from hornwort crown age bounds.
All times in units of 100 Ma for MCMCtree.
"""

import os

BASE = "/media/data/projects/hornwort_sex_chromosomes/analysis/gene_analyses/gametolog_ancestry_tests"
FAMSA_DIR = os.path.join(BASE, "hog_alignments")
OUT_DIR = os.path.join(BASE, "mcmctree")

# --- Calibrations (100 Ma units, soft bounds B(min,max)) ---
CAL = {
    "root": "B(3.49,5.05)",      # U/V divergence: min=hornwort crown (348.56 Ma), max=bryophyte split (~505 Ma)
    "A":    "B(3.3071,3.4856)",  # Hornwort crown: Leiosporoceros vs rest (330.71-348.56 Ma)
    "B":    "B(2.9197,3.0748)",  # Anthoceros vs (Notothyladaceae+PhymatDendro) (291.97-307.48 Ma)
    "C":    "B(1.7075,1.9355)",  # Notothyladaceae+PhymatDendro crown: Papro vs Phphy clade (170.75-193.55 Ma)
    "D":    "B(1.2010,1.3949)",  # PhymatDendro crown: Phphy vs Dendrocerotaceae (120.10-139.49 Ma)
    "E":    "B(0.7598,0.8719)",  # Dendrocerotaceae crown: Phchi vs Noaen (75.98-87.19 Ma)
}

# --- Focal sequences per HOG ---
# short_name -> exact sequence name in FAMSA file
HOG_SEQS = {
    "OG0000460": {
        "LedusF": "LedusF.S3G012600.t2",
        "LedusM": "LedusM.S5G005100.t1",
        "PhphyF": "PhphyF.S5G024600.t1",
        "PhphyM": "PhphyM.S5G054900.t1",
        "NoaenF": "NoaenF.S5G001400.t2",
        "PhchiF": "PhchiF.S4G011700.t2",
        "PhchiM": "PhchiM.S6G007400.t1",
        "PaproF": "PaproF.S5G042600.t2",
        "PaproM": "PaproM.S5G006700.t3",
    },
    "OG0000985": {
        "LedusF": "LedusF.S3G014200.t1",
        "LedusM": "LedusM.S5G003400.t2",
        "PhphyF": "PhphyF.S5G055300.t2",
        "PhphyM": "PhphyM.S5G025900.t2",
        "NoaenF": "NoaenF.S5G051400.t2",
        "PhchiF": "PhchiF.S4G008200.t1",
        "PhchiM": "PhchiM.S6G034100.t2",
        "PaproF": "PaproF.S5G021100.t2",
        "PaproM": "PaproM.S5G003800.t2",
    },
    "OG0001152": {
        "LedusF": "LedusF.S3G012200.t2",
        "LedusM": "LedusM.S5G019400.t2",
        # PhphyF excluded (intruder in gene tree; groups with V sequences)
        "PhphyM": "PhphyM.S5G004600.t4",
        "NoaenF": "NoaenF.S5G017100.t1",
        "PhchiF": "PhchiF.S4G001500.t2",
        "PhchiM": "PhchiM.S6G018100.t2",
        "PaproF": "PaproF.S5G030500.t1",
        "PaproM": "PaproM.S5G005600.t1",
    },
    "OG0001301": {
        "LedusF": "LedusF.S3G000300.t1",
        "LedusM": "LedusM.S5G012000.t1",
        # No PhphyF or PhphyM in this OG
        "NoaenF": "NoaenF.S5G011000.t2",
        "PhchiF": "PhchiF.S4G001100.t1",  # tandem dup; using S4G001100 (first copy)
        "PhchiM": "PhchiM.S6G032200.t1",
        "PaproF": "PaproF.S5G037400.t1",
        "PaproM": "PaproM.S5G007000.t1",
        "AnangF": "AnangRef2_g18449.t3",  # confirmed U-linked
        "AnangM": "AnangRef2_g18616.t1",  # confirmed V-linked
    },
    "OG0002583": {  # U-only HOG
        "LedusF": "LedusF.S3G001600.t1",
        "PhphyF": "PhphyF.S5G006900.t1",
        "NoaenF": "NoaenF.S5G048500.t2",
        "PhchiF": "PhchiF.S4G002300.t1",
        "PaproF": "PaproF.S5G031700.t2",
        "AnangF": "AnangRef2_g18296.t1",  # U-only HOG; included as U
    },
    "OG0002735": {  # LedusM absent from OG
        "LedusF": "LedusF.S4G015500.t1",
        # LedusM absent from this OG
        "PhphyF": "PhphyF.S5G059500.t1",
        "PhphyM": "PhphyM.S5G005400.t1",
        "NoaenF": "NoaenF.S5G056800.t1",
        "PhchiF": "PhchiF.S4G010100.t1",
        "PhchiM": "PhchiM.S6G005900.t1",
        "PaproF": "PaproF.S5G035700.t1",
        "PaproM": "PaproM.S5G006100.t1",
        "AnangF": "AnangRef2_g18376.t1",  # confirmed U-linked
        "AnangM": "AnangRef2_g18716.t1",  # confirmed V-linked
    },
    "OG0003154": {  # PhchiM absent from OG
        "LedusF": "LedusF.S3G010100.t1",
        "LedusM": "LedusM.S5G019500.t1",
        "PhphyF": "PhphyF.S5G053000.t1",
        "PhphyM": "PhphyM.S5G028500.t1",
        "NoaenF": "NoaenF.S5G005600.t1",
        "PhchiF": "PhchiF.S4G000500.t1",
        # PhchiM absent from this OG
        "PaproF": "PaproF.S5G036300.t2",
        "PaproM": "PaproM.S5G024800.t1",
    },
    "OG0008579": {  # U-only HOG
        "LedusF": "LedusF.S3G013600.t1",
        "PhphyF": "PhphyF.S5G052100.t1",
        "NoaenF": "NoaenF.S5G077200.t2",
        "PhchiF": "PhchiF.S4G016400.t2",
        "PaproF": "PaproF.S5G031100.t1",
        "AnangF": "AnangRef2_g18310.t1",  # U-only HOG; included as U
    },
}

def read_fasta(path):
    seqs, order, cur = {}, [], None
    with open(path) as f:
        for line in f:
            line = line.rstrip()
            if line.startswith('>'):
                cur = line[1:].split()[0]
                order.append(cur)
                seqs[cur] = ''
            elif cur:
                seqs[cur] += line
    return order, seqs

def strip_gap_cols(seqs_dict):
    names = list(seqs_dict.keys())
    length = len(seqs_dict[names[0]])
    keep = [i for i in range(length) if any(seqs_dict[n][i] != '-' for n in names)]
    return {n: ''.join(seqs_dict[n][i] for i in keep) for n in names}, len(keep)

def build_clade(taxa_set, is_u, cal):
    """
    Build calibrated newick string for one clade (U or V).
    Species tree topology: (Ledus, (Anang, (Papro, (Phphy, (Noaen, Phchi)))))
    Noaen only in U-clade; calibrations applied where taxa are present.
    """
    ledus = "LedusF" if is_u else "LedusM"
    anang = "AnangF" if is_u else "AnangM"
    papro = "PaproF" if is_u else "PaproM"
    phphy = "PhphyF" if is_u else "PhphyM"
    phchi = "PhchiF" if is_u else "PhchiM"
    noaen = "NoaenF"

    has = {t: t in taxa_set for t in [ledus, anang, papro, phphy, phchi]}
    has[noaen] = (noaen in taxa_set) and is_u

    # Innermost: Dendrocerotaceae (Noaen+Phchi in U; just Phchi in V)
    if has[noaen] and has[phchi]:
        inner = f"({noaen},{phchi})'{cal['E']}'"
    elif has[phchi]:
        inner = phchi
    elif has[noaen]:
        inner = noaen
    else:
        inner = None

    # PhymatDendro crown: Phphy + Dendrocerotaceae
    if has[phphy] and inner:
        inner = f"({phphy},{inner})'{cal['D']}'"
    elif has[phphy]:
        inner = phphy

    # Notothyladaceae+PhymatDendro crown: Papro + above
    if has[papro] and inner:
        inner = f"({papro},{inner})'{cal['C']}'"
    elif has[papro]:
        inner = papro

    # Anthoceros node: Anang + above
    if has[anang] and inner:
        inner = f"({anang},{inner})'{cal['B']}'"
    elif has[anang]:
        inner = anang

    # Hornwort crown: Ledus + above
    if has[ledus] and inner:
        inner = f"({ledus},{inner})'{cal['A']}'"
    elif has[ledus]:
        inner = ledus

    return inner

def build_tree(taxa_set, cal):
    """Build the complete gametolog tree with U and V clades."""
    v_taxa = {t for t in taxa_set if t.endswith('M')}
    u_clade = build_clade(taxa_set, is_u=True, cal=cal)

    if not v_taxa:
        # U-only HOG: single clade, no UV root calibration
        return f"{u_clade};"
    else:
        v_clade = build_clade(taxa_set, is_u=False, cal=cal)
        return f"({u_clade},{v_clade})'{cal['root']}';"

# MCMCtree control file template (Step 1: approximate likelihood)
CTL_STEP1 = """\
          seqfile = {seqfile}
         treefile = {treefile}
          outfile = step1_out

          seqtype = 2      * 2 = amino acid
        cleandata = 0
          usedata = 3      * 3 = compute gradient/Hessian for approx likelihood

            model = 2      * 2 = Empirical (rate matrix from aaRatefile)
       aaRatefile = /usr/lib/paml/data/dat/lg.dat   * LG model
            alpha = 0.5
           ncatG = 4

            clock = 2      * 2 = independent log-normal rates (uncorrelated)

      BDparas = 1 1 0    * birth, death, sampling (for time prior)
  rgene_gamma = 2 4      * gamma prior on mean rate (mean=0.5 sub/site/100Ma)
 sigma2_gamma = 1 10     * gamma prior on sigma2 (mean=0.1, for ILN variance)

         burnin = 2000
       sampfreq = 5
        nsample = 20000
"""

# MCMCtree control file template (Step 2: MCMC dating)
CTL_STEP2 = """\
          seqfile = {seqfile}
         treefile = {treefile}
          outfile = step2_out

          seqtype = 2
        cleandata = 0
          usedata = 2      * 2 = use approximate likelihood (reads in.BV from step1 out.BV)

            model = 2
       aaRatefile = /usr/lib/paml/data/dat/lg.dat
            alpha = 0.5
           ncatG = 4

            clock = 2

      BDparas = 1 1 0
  rgene_gamma = 2 4
 sigma2_gamma = 1 10

         print = 1      * 1 = write MCMC samples to mcmc.txt

         burnin = 50000
       sampfreq = 5
        nsample = 100000
"""

# --- Main ---
os.makedirs(OUT_DIR, exist_ok=True)

print(f"{'HOG':<15} {'taxa':>5} {'sites':>6}  Tree")
print("-" * 120)

for hog, seq_map in HOG_SEQS.items():
    famsa_path = os.path.join(FAMSA_DIR, f"HOG_{hog}.FAMSA.fa")
    hog_dir = os.path.join(OUT_DIR, f"HOG_{hog}")
    os.makedirs(os.path.join(hog_dir, "step1"), exist_ok=True)
    os.makedirs(os.path.join(hog_dir, "step2"), exist_ok=True)

    # Read FAMSA alignment
    _, seqs = read_fasta(famsa_path)

    # Map focal sequences
    focal = {}
    for short, famsa_name in seq_map.items():
        if famsa_name in seqs:
            focal[short] = seqs[famsa_name]
        else:
            print(f"  WARNING: {famsa_name} not in {hog}")

    # Strip gap-only columns from focal subset
    stripped, aln_len = strip_gap_cols(focal)
    n_taxa = len(stripped)

    # Write PHYLIP alignment
    phy_name = f"HOG_{hog}_focal.phy"
    phy_path = os.path.join(hog_dir, phy_name)
    with open(phy_path, 'w') as f:
        f.write(f"  {n_taxa}  {aln_len}\n")
        for name, seq in stripped.items():
            f.write(f"{name:<10}  {seq}\n")

    # Build and write tree
    tree_str = build_tree(set(stripped.keys()), CAL)
    tre_name = f"HOG_{hog}.tre"
    tre_path = os.path.join(hog_dir, tre_name)
    with open(tre_path, 'w') as f:
        f.write(f"  {n_taxa}  1\n")
        f.write(f"{tree_str}\n")

    # Write control files (relative paths so mcmctree is run from step1/ or step2/)
    phy_rel = f"../HOG_{hog}_focal.phy"
    tre_rel = f"../HOG_{hog}.tre"

    with open(os.path.join(hog_dir, "step1", "mcmctree.ctl"), 'w') as f:
        f.write(CTL_STEP1.format(seqfile=phy_rel, treefile=tre_rel))
    with open(os.path.join(hog_dir, "step2", "mcmctree.ctl"), 'w') as f:
        f.write(CTL_STEP2.format(seqfile=phy_rel, treefile=tre_rel))

    print(f"{hog:<15} {n_taxa:>5} {aln_len:>6}  {tree_str}")
