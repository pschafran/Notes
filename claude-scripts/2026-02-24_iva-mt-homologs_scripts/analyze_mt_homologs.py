#!/usr/bin/env python3
"""
analyze_mt_homologs.py  — Cluster mt scaffolds into homolog groups via BLAST graph,
extract per-cluster representatives, align with MAFFT, and report summary statistics.

Usage:
    python3 analyze_mt_homologs.py

Inputs  (expected in same directory as this script):
    all_mt_scaffolds.fasta
    mt_allvsall.blastn

Outputs:
    mt_clusters.tsv          — per-scaffold cluster assignment
    mt_cluster_stats.tsv     — per-cluster statistics
    mt_summary.txt           — printed summary report
    cluster_NN/              — one subdirectory per top cluster
        cluster_NN_seqs.fasta    — one representative per sample
        cluster_NN_aligned.fasta — MAFFT alignment

Filters applied to BLAST hits before graph construction:
    - exclude same-sequence hits (qseqid == sseqid)
    - exclude same-sample hits
    - pident >= 70 (already filtered by blastn -perc_identity)
    - alignment length >= 500 bp
    - alignment length >= 50% of the shorter sequence length
"""

import os
import sys
import math
import subprocess
import collections
from pathlib import Path

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------
SCRIPT_DIR   = Path(__file__).parent
FASTA_FILE   = SCRIPT_DIR / "all_mt_scaffolds.fasta"
BLAST_FILE   = SCRIPT_DIR / "mt_allvsall.blastn"
CLUSTER_TSV  = SCRIPT_DIR / "mt_clusters.tsv"
STATS_TSV    = SCRIPT_DIR / "mt_cluster_stats.tsv"
SUMMARY_FILE = SCRIPT_DIR / "mt_summary.txt"

MIN_ALN_LEN   = 500    # bp
MIN_COV_FRAC  = 0.50   # fraction of shorter sequence
N_SAMPLES     = 23
CORE_THRESH   = 0.80   # fraction of samples for "core" cluster
TOP_N_ALIGN   = 20     # align top N clusters by n_samples then mean_length
MAFFT_BIN     = "mafft"

# Samples whose mt assembly is a single large repeat-path sequence (not individual scaffolds).
# These are excluded from MSA alignment since they cannot be meaningfully globally aligned
# with individual scaffolds from other samples.
REPEAT_PATH_SAMPLES = {"axillaris_89", "frutescens_36", "frutescens_49"}
MAX_REP_LEN_FOR_ALIGN = 100_000  # bp; representatives larger than this are skipped in MSA

# ---------------------------------------------------------------------------
# Union-Find
# ---------------------------------------------------------------------------
class UnionFind:
    def __init__(self):
        self.parent = {}
        self.rank   = {}

    def find(self, x):
        self.parent.setdefault(x, x)
        self.rank.setdefault(x, 0)
        if self.parent[x] != x:
            self.parent[x] = self.find(self.parent[x])
        return self.parent[x]

    def union(self, x, y):
        rx, ry = self.find(x), self.find(y)
        if rx == ry:
            return
        if self.rank[rx] < self.rank[ry]:
            rx, ry = ry, rx
        self.parent[ry] = rx
        if self.rank[rx] == self.rank[ry]:
            self.rank[rx] += 1

# ---------------------------------------------------------------------------
# 1. Load sequences
# ---------------------------------------------------------------------------
def load_fasta(fasta_path):
    """Return dict: seq_id -> sequence string, and preserve order."""
    seqs  = {}
    order = []
    cur_id = None
    buf    = []
    with open(fasta_path) as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith(">"):
                if cur_id is not None:
                    seqs[cur_id] = "".join(buf)
                cur_id = line[1:].split()[0]  # first word after >
                order.append(cur_id)
                buf = []
            else:
                buf.append(line)
    if cur_id is not None:
        seqs[cur_id] = "".join(buf)
    return seqs, order

def sample_of(seq_id):
    """Extract sample prefix from 'SAMPLE|header'."""
    return seq_id.split("|")[0]

# ---------------------------------------------------------------------------
# 2. Parse BLAST and build Union-Find
# ---------------------------------------------------------------------------
def parse_blast_and_cluster(blast_path, seq_lengths):
    uf     = UnionFind()
    # Register all sequences as nodes
    for sid in seq_lengths:
        uf.find(sid)

    edge_count = 0
    with open(blast_path) as fh:
        for line in fh:
            cols = line.rstrip().split("\t")
            qid, sid2 = cols[0], cols[1]
            pident     = float(cols[2])
            aln_len    = int(cols[3])
            qlen       = int(cols[4])
            slen       = int(cols[5])

            # Skip self-hit
            if qid == sid2:
                continue
            # Skip same-sample
            if sample_of(qid) == sample_of(sid2):
                continue
            # Alignment length filter
            if aln_len < MIN_ALN_LEN:
                continue
            # Coverage filter: aln_len >= 50% of shorter sequence
            shorter = min(qlen, slen)
            if aln_len < MIN_COV_FRAC * shorter:
                continue

            uf.union(qid, sid2)
            edge_count += 1

    return uf, edge_count

# ---------------------------------------------------------------------------
# 3. Build cluster membership
# ---------------------------------------------------------------------------
def build_clusters(uf, seq_order):
    clusters = collections.defaultdict(list)  # root -> [seq_ids]
    for sid in seq_order:
        root = uf.find(sid)
        clusters[root].append(sid)
    return clusters

# ---------------------------------------------------------------------------
# 4. Cluster statistics
# ---------------------------------------------------------------------------
def cluster_stats(clusters, seq_lengths):
    stats = []
    for root, members in clusters.items():
        samples   = sorted(set(sample_of(m) for m in members))
        n_samp    = len(samples)
        lengths   = [seq_lengths[m] for m in members]
        stats.append({
            "root"       : root,
            "n_seqs"     : len(members),
            "n_samples"  : n_samp,
            "samples"    : samples,
            "members"    : members,
            "lengths"    : lengths,
            "min_len"    : min(lengths),
            "max_len"    : max(lengths),
            "mean_len"   : sum(lengths) / len(lengths),
            "total_bp"   : sum(lengths),
        })
    # Sort: n_samples DESC, mean_len DESC
    stats.sort(key=lambda x: (-x["n_samples"], -x["mean_len"]))
    # Assign cluster IDs
    for i, s in enumerate(stats, 1):
        s["cluster_id"] = f"cluster_{i:02d}"
    return stats

# ---------------------------------------------------------------------------
# 5. Extract representatives (longest per sample) and write FASTA
# ---------------------------------------------------------------------------
def extract_representatives(cluster_stat, seq_dict, out_dir):
    """Write one representative sequence per sample for this cluster.

    Returns:
        out_fasta  — path to FASTA with all reps (including large ones)
        reps       — all representative seq ids
        reps_align — representative seq ids suitable for MSA (size-filtered)
        skipped    — list of (sample, rep_len) excluded from MSA
    """
    out_dir.mkdir(parents=True, exist_ok=True)
    members  = cluster_stat["members"]
    by_sample = collections.defaultdict(list)
    for m in members:
        by_sample[sample_of(m)].append(m)

    cid       = cluster_stat["cluster_id"]
    out_fasta  = out_dir / f"{cid}_seqs.fasta"
    out_fasta_align = out_dir / f"{cid}_seqs_for_align.fasta"

    reps        = []
    reps_align  = []
    skipped     = []

    with open(out_fasta, "w") as fh_all, open(out_fasta_align, "w") as fh_aln:
        for sample, seqids in sorted(by_sample.items()):
            rep = max(seqids, key=lambda x: len(seq_dict[x]))
            reps.append(rep)
            rep_len = len(seq_dict[rep])
            seq = seq_dict[rep]
            # Write to all-reps FASTA always
            fh_all.write(f">{rep}\n")
            for i in range(0, len(seq), 60):
                fh_all.write(seq[i:i+60] + "\n")
            # Write to alignment FASTA only if not in repeat_path set and not too large
            if sample in REPEAT_PATH_SAMPLES or rep_len > MAX_REP_LEN_FOR_ALIGN:
                skipped.append((sample, rep_len))
            else:
                reps_align.append(rep)
                fh_aln.write(f">{rep}\n")
                for i in range(0, len(seq), 60):
                    fh_aln.write(seq[i:i+60] + "\n")

    # Remove the align FASTA if empty (so we don't accidentally run mafft on nothing)
    if not reps_align:
        out_fasta_align.unlink(missing_ok=True)
        out_fasta_align = None

    return out_fasta, reps, reps_align, skipped, out_fasta_align

# ---------------------------------------------------------------------------
# 6. MAFFT alignment
# ---------------------------------------------------------------------------
def run_mafft(in_fasta, out_fasta, n_threads=32):
    cmd = [MAFFT_BIN, "--auto", "--thread", str(n_threads), str(in_fasta)]
    try:
        with open(out_fasta, "w") as fh:
            result = subprocess.run(cmd, stdout=fh, stderr=subprocess.PIPE,
                                    check=True, text=True)
        return True, result.stderr
    except subprocess.CalledProcessError as e:
        return False, e.stderr
    except FileNotFoundError:
        return False, f"mafft not found at '{MAFFT_BIN}'"

# ---------------------------------------------------------------------------
# 7. Parsimony-informative sites
# ---------------------------------------------------------------------------
def count_pi_sites(aligned_fasta):
    """Count parsimony-informative sites in a MAFFT alignment."""
    seqs = {}
    order = []
    cur_id = None
    buf = []
    with open(aligned_fasta) as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith(">"):
                if cur_id is not None:
                    seqs[cur_id] = "".join(buf)
                cur_id = line[1:].split()[0]
                order.append(cur_id)
                buf = []
            else:
                buf.append(line)
    if cur_id is not None:
        seqs[cur_id] = "".join(buf)

    if len(seqs) < 2:
        return 0, 0, 0

    n_taxa = len(seqs)
    aln_len = len(next(iter(seqs.values())))
    pi_count = 0
    for col in range(aln_len):
        states = collections.Counter()
        for sid in order:
            ch = seqs[sid][col].upper()
            if ch not in ("-", "N", "?"):
                states[ch] += 1
        # PI site: >= 2 character states, each with >= 2 taxa
        informative = sum(1 for cnt in states.values() if cnt >= 2)
        if informative >= 2:
            pi_count += 1

    return aln_len, pi_count, round(100.0 * pi_count / aln_len, 2) if aln_len > 0 else 0

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    print("=== Mt Homolog Clustering Analysis ===\n")

    # Load sequences
    print(f"Loading sequences from {FASTA_FILE.name} ...")
    seq_dict, seq_order = load_fasta(FASTA_FILE)
    seq_lengths = {sid: len(seq) for sid, seq in seq_dict.items()}
    print(f"  {len(seq_dict)} sequences from {len(set(sample_of(s) for s in seq_dict))} samples")

    # Parse BLAST and cluster
    print(f"\nParsing BLAST hits from {BLAST_FILE.name} ...")
    print(f"  Filters: aln_len >= {MIN_ALN_LEN} bp, coverage >= {MIN_COV_FRAC*100:.0f}% of shorter seq")
    uf, n_edges = parse_blast_and_cluster(BLAST_FILE, seq_lengths)
    print(f"  Passed filter: {n_edges} edges used for clustering")

    # Build clusters
    clusters = build_clusters(uf, seq_order)
    print(f"  Connected components (clusters): {len(clusters)}")

    # Compute statistics
    stats = cluster_stats(clusters, seq_lengths)

    # Write per-scaffold cluster assignments
    with open(CLUSTER_TSV, "w") as fh:
        fh.write("seq_id\tsample\tcluster_id\tseq_len\tn_samples_in_cluster\n")
        for s in stats:
            for m in s["members"]:
                fh.write(f"{m}\t{sample_of(m)}\t{s['cluster_id']}\t{seq_lengths[m]}\t{s['n_samples']}\n")
    print(f"  Per-scaffold assignments written to {CLUSTER_TSV.name}")

    # Write cluster stats TSV
    core_min = math.ceil(N_SAMPLES * CORE_THRESH)
    with open(STATS_TSV, "w") as fh:
        fh.write("cluster_id\tn_seqs\tn_samples\tmin_len\tmax_len\tmean_len\ttotal_bp\tcore\tsamples\n")
        for s in stats:
            core = "yes" if s["n_samples"] >= core_min else "no"
            fh.write(
                f"{s['cluster_id']}\t{s['n_seqs']}\t{s['n_samples']}\t"
                f"{s['min_len']}\t{s['max_len']}\t{s['mean_len']:.0f}\t"
                f"{s['total_bp']}\t{core}\t{','.join(s['samples'])}\n"
            )
    print(f"  Cluster stats written to {STATS_TSV.name}")

    # Summary counts
    n_core = sum(1 for s in stats if s["n_samples"] >= core_min)
    n_singleton = sum(1 for s in stats if s["n_samples"] == 1)
    print(f"\n  Core clusters (>= {core_min}/{N_SAMPLES} samples): {n_core}")
    print(f"  Singleton clusters (1 sample):                    {n_singleton}")

    # Show top-20 clusters
    print(f"\n{'Cluster':<12} {'Seqs':>5} {'Samples':>8} {'MinLen':>8} {'MaxLen':>8} {'MeanLen':>8}")
    print("-" * 62)
    for s in stats[:20]:
        core_tag = " *" if s["n_samples"] >= core_min else ""
        print(f"{s['cluster_id']:<12} {s['n_seqs']:>5} {s['n_samples']:>8} "
              f"{s['min_len']:>8} {s['max_len']:>8} {s['mean_len']:>8.0f}{core_tag}")

    # Extract and align top clusters
    print(f"\n=== Extracting and aligning top {TOP_N_ALIGN} clusters ===\n")
    print(f"  Note: {sorted(REPEAT_PATH_SAMPLES)} excluded from MSA (repeat-path assemblies >100 kb)")
    print()
    alignment_results = []
    for s in stats[:TOP_N_ALIGN]:
        cid = s["cluster_id"]
        out_dir = SCRIPT_DIR / cid
        out_fasta, reps, reps_align, skipped, fasta_for_align = extract_representatives(
            s, seq_dict, out_dir)
        n_reps = len(reps)
        n_align = len(reps_align)

        skip_note = ""
        if skipped:
            skip_note = f" (skipped from MSA: {[sa for sa, _ in skipped]})"

        if n_align < 2:
            print(f"  {cid}: {n_reps} total reps, {n_align} alignable — skipping MSA{skip_note}")
            alignment_results.append((cid, s["n_samples"], n_reps, n_align, None, None, None, skipped))
            continue

        aligned_fasta = out_dir / f"{cid}_aligned.fasta"
        print(f"  {cid}: aligning {n_align}/{n_reps} sequences{skip_note} ...", end=" ", flush=True)
        ok, stderr = run_mafft(fasta_for_align, aligned_fasta)
        if ok:
            aln_len, pi_sites, pi_pct = count_pi_sites(aligned_fasta)
            print(f"OK  aln_len={aln_len:,}  PI={pi_sites:,} ({pi_pct}%)")
            alignment_results.append((cid, s["n_samples"], n_reps, n_align, aln_len, pi_sites, pi_pct, skipped))
        else:
            print(f"FAILED: {stderr[:120]}")
            alignment_results.append((cid, s["n_samples"], n_reps, n_align, None, None, None, skipped))

    # Write summary report
    lines = []
    lines.append("=== Mt Homolog Clustering Summary ===\n")
    lines.append(f"Total sequences  : {len(seq_dict)}")
    lines.append(f"Total samples    : {len(set(sample_of(s) for s in seq_dict))}")
    lines.append(f"Total clusters   : {len(stats)}")
    lines.append(f"Core clusters    : {n_core} (>= {core_min}/{N_SAMPLES} samples, {CORE_THRESH*100:.0f}% threshold)")
    lines.append(f"Singleton clusters: {n_singleton}")
    lines.append("")
    lines.append("--- Top 20 clusters ---")
    lines.append(f"{'Cluster':<12} {'Seqs':>5} {'Smpls':>6} {'MinLen':>8} {'MaxLen':>8} "
                 f"{'MeanLen':>8} {'Alignd':>7} {'Core':>5} {'AlnLen':>10} {'PI%':>6}")
    lines.append("-" * 80)
    for res, s in zip(alignment_results, stats[:TOP_N_ALIGN]):
        cid, n_samp, n_rep, n_align, aln_len, pi_sites, pi_pct, skipped = res
        core_str = "yes" if n_samp >= core_min else "no"
        aln_str  = f"{aln_len:,}" if aln_len else "—"
        pi_str   = f"{pi_pct}%" if pi_pct is not None else "—"
        lines.append(
            f"{cid:<12} {s['n_seqs']:>5} {n_samp:>6} {s['min_len']:>8} {s['max_len']:>8} "
            f"{s['mean_len']:>8.0f} {n_align:>7} {core_str:>5} {aln_str:>10} {pi_str:>6}"
        )
        if skipped:
            lines.append(f"  ^ excluded from MSA: " + ", ".join(f"{sa}({sl:,}bp)" for sa,sl in skipped))
    lines.append("")
    lines.append("--- Cluster coverage thresholds ---")
    for thresh in [0.5, 0.6, 0.7, 0.8, 0.9, 1.0]:
        n_min = math.ceil(N_SAMPLES * thresh)
        n_at  = sum(1 for s in stats if s["n_samples"] >= n_min)
        lines.append(f"  >= {thresh*100:.0f}% of samples ({n_min}/{N_SAMPLES}): {n_at} clusters")

    report = "\n".join(lines) + "\n"
    print("\n" + report)
    with open(SUMMARY_FILE, "w") as fh:
        fh.write(report)
    print(f"Summary written to {SUMMARY_FILE.name}")

if __name__ == "__main__":
    main()
