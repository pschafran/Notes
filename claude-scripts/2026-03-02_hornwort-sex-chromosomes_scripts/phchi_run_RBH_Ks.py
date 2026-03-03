#!/usr/bin/env python3
"""
Reciprocal-Best-BLAST hits + Ks analysis
Phaeomegaceros chiloensis female (PhchiF/14765-5) vs male (PhchiM/14765-4)

Pipeline:
  1. Build blastp databases for primary-transcript proteins
  2. Run blastp in both directions (F→M, M→F)
  3. Find reciprocal best hits (RBH)
  4. For each RBH pair: align proteins, back-translate, run yn00 for Ka/Ks
  5. Write results TSV + plain-text caption

Ks method: Yang & Nielsen 2000 (YN00) via external /usr/bin/yn00 binary.
NOTE: BioPython 1.85 cal_dn_ds(method='YN00') is broken in Python 3.10
      (RuntimeError: dictionary changed size during iteration). Use external yn00.
"""

import datetime
import math
import multiprocessing
import os
import re
import shutil
import subprocess
import tempfile
from collections import defaultdict

from Bio import SeqIO
from Bio.Align import PairwiseAligner, substitution_matrices
from Bio.Seq import Seq

# ── Paths ─────────────────────────────────────────────────────────────────────
BASE    = '/media/data/projects/hornwort_sex_chromosomes/analysis/Phaeomegaceros_fimbriatus'
WORKDIR = os.path.join(BASE, 'sex_chromosome_analyses')
OUTDIR  = os.path.join(WORKDIR, 'rbh_ks')

F_PROT = os.path.join(BASE, 'Phaeomegaceros_14765-5/final_genome_prep/braker_renamed_primary_transcripts.aa')
M_PROT = os.path.join(BASE, 'Phaeomegaceros_14765-4/final_genome_prep/braker_renamed_primary_transcripts.aa')
F_CDS  = os.path.join(BASE, 'Phaeomegaceros_14765-5/final_genome_prep/braker_renamed_primary_transcripts.codingseq')
M_CDS  = os.path.join(BASE, 'Phaeomegaceros_14765-4/final_genome_prep/braker_renamed_primary_transcripts.codingseq')

DB_F        = os.path.join(OUTDIR, 'db_PhchiF')
DB_M        = os.path.join(OUTDIR, 'db_PhchiM')
BLAST_FM    = os.path.join(OUTDIR, 'blast_F_vs_M.tsv')
BLAST_MF    = os.path.join(OUTDIR, 'blast_M_vs_F.tsv')
TSV_OUT     = os.path.join(OUTDIR, 'PhchiF_PhchiM_RBH_Ks.tsv')
CAPTION_OUT = os.path.join(OUTDIR, 'PhchiF_PhchiM_RBH_Ks.caption.txt')

YN00_BIN   = '/usr/bin/yn00'
MIN_CODONS = 50
EVALUE     = 1e-5
MAX_HITS   = 5


def extract_scaffold(gene_id):
    """
    PhchiF.S4G000100.t1  ->  PhchiF.S4
    PhchiM.C10G000100.t1 -> PhchiM.C10
    """
    m = re.match(r'^(.*?)G\d+\.t\d+$', gene_id)
    return m.group(1) if m else gene_id


def prepare_cds(raw_cds):
    cds = raw_cds.upper().replace('\n', '').replace(' ', '')
    if len(cds) < 3:
        raise ValueError(f'CDS too short ({len(cds)} bp)')
    if cds[-3:] in ('TAA', 'TAG', 'TGA'):
        cds = cds[:-3]
    if len(cds) == 0:
        raise ValueError('CDS is empty after stop codon removal')
    if len(cds) % 3 != 0:
        raise ValueError(f'CDS length {len(cds)} not divisible by 3')
    translated = str(Seq(cds).translate(to_stop=False))
    if '*' in translated:
        raise ValueError('Internal stop codon in CDS')
    return cds


def align_proteins(p1, p2):
    aligner = PairwiseAligner()
    aligner.substitution_matrix = substitution_matrices.load('BLOSUM62')
    aligner.open_gap_score    = -11
    aligner.extend_gap_score  = -1
    aligner.mode              = 'global'
    alignments = aligner.align(p1, p2)
    return next(iter(alignments))


def backtranslate(prot_aln, cds):
    codons = [cds[i:i+3] for i in range(0, len(cds), 3)]
    result = []
    ci = 0
    for aa in str(prot_aln):
        if aa == '-':
            result.append('---')
        else:
            result.append(codons[ci])
            ci += 1
    return ''.join(result)


def run_yn00(cds1, cds2):
    """Run external yn00 binary; return (Ka, Ks) or raise."""
    tmpdir = tempfile.mkdtemp(prefix='yn00_')
    try:
        n_bp = len(cds1)
        with open(os.path.join(tmpdir, 'aln.phy'), 'w') as fh:
            fh.write(f' 2 {n_bp}\n')
            fh.write(f'seq1      {cds1}\n')
            fh.write(f'seq2      {cds2}\n')

        with open(os.path.join(tmpdir, 'yn00.ctl'), 'w') as fh:
            fh.write('seqfile = aln.phy\n')
            fh.write('outfile = yn00.out\n')
            fh.write('verbose = 0\n')
            fh.write('icode = 0\n')
            fh.write('weighting = 0\n')
            fh.write('commonf3x4 = 0\n')

        proc = subprocess.run(
            [YN00_BIN],
            capture_output=True,
            cwd=tmpdir, timeout=60,
        )
        if proc.returncode != 0:
            raise RuntimeError(f'yn00 exit {proc.returncode}')

        out_file = os.path.join(tmpdir, 'yn00.out')
        with open(out_file) as fh:
            content = fh.read()
        if not content.strip():
            raise RuntimeError('yn00 output is empty')

        # Parse Yang & Nielsen (2000) table
        # Data line: seq_i seq_j  S  N  t  kappa  omega  dN +- SE  dS +- SE
        in_yn00 = False
        num = r'-?[\d.]+'
        for line in content.splitlines():
            if '(B) Yang & Nielsen' in line:
                in_yn00 = True
                continue
            if in_yn00 and line.strip().startswith('(C)'):
                break
            if in_yn00:
                m = re.match(
                    rf'^\s*\d+\s+\d+\s+'
                    rf'{num}\s+{num}\s+'
                    rf'{num}\s+{num}\s+{num}\s+'
                    rf'({num})\s+\+-\s+{num}\s+'
                    rf'({num})\s+\+-',
                    line,
                )
                if m:
                    return float(m.group(1)), float(m.group(2))

        raise RuntimeError('Could not parse yn00 output')
    finally:
        shutil.rmtree(tmpdir, ignore_errors=True)


def process_pair(args):
    idx, total, gid_F, gid_M, f_prot_seq, m_prot_seq, f_cds_seq, m_cds_seq = args
    result = {
        'gene_F': gid_F, 'gene_M': gid_M,
        'scaffold_F': extract_scaffold(gid_F),
        'scaffold_M': extract_scaffold(gid_M),
        'len_F': len(f_prot_seq), 'len_M': len(m_prot_seq),
        'Ka': float('nan'), 'Ks': float('nan'), 'Ka_Ks': float('nan'),
        'flag': 'error',
    }
    try:
        p_F = str(f_prot_seq.seq).rstrip('*').replace('*', 'X')
        p_M = str(m_prot_seq.seq).rstrip('*').replace('*', 'X')
        aln = align_proteins(p_F, p_M)
        f_aln = str(aln[0])
        m_aln = str(aln[1])

        f_cds_raw = str(f_cds_seq.seq)
        m_cds_raw = str(m_cds_seq.seq)
        f_cds = prepare_cds(f_cds_raw)
        m_cds = prepare_cds(m_cds_raw)

        f_bt = backtranslate(f_aln, f_cds)
        m_bt = backtranslate(m_aln, m_cds)

        if len(f_bt) != len(m_bt):
            result['flag'] = 'bt_mismatch'
            return result

        n_codons = len(f_bt) // 3
        if n_codons < MIN_CODONS:
            result['flag'] = 'short_aln'
            return result

        ka, ks = run_yn00(f_bt, m_bt)
        omega = ka / ks if ks > 0 else float('nan')

        result['Ka']    = ka
        result['Ks']    = ks
        result['Ka_Ks'] = omega

        if ks == 0 or not math.isfinite(ks):
            result['flag'] = 'zero_Ks'
        elif ks < 0 or ka < 0:
            result['flag'] = 'negative'
        else:
            result['flag'] = 'ok'

    except ValueError as e:
        result['flag'] = f'cds_err'
    except RuntimeError as e:
        result['flag'] = f'yn00_err'
    except Exception as e:
        result['flag'] = f'error'

    return result


def main():
    os.makedirs(OUTDIR, exist_ok=True)
    n_cpu = multiprocessing.cpu_count()
    print(f'CPUs available: {n_cpu}')

    # ── [1/5] BLAST databases ──────────────────────────────────────────────
    print('\n[1/5] Building BLAST databases...')
    for label, fasta, db in [('PhchiF', F_PROT, DB_F), ('PhchiM', M_PROT, DB_M)]:
        print(f'  makeblastdb: {os.path.basename(fasta)} → {os.path.basename(db)}')
        subprocess.run(
            ['makeblastdb', '-in', fasta, '-dbtype', 'prot', '-out', db],
            check=True, capture_output=True,
        )

    # ── [2/5] blastp ──────────────────────────────────────────────────────
    print('\n[2/5] Running blastp...')
    for label, query, db, out in [
        ('F→M', F_PROT, DB_M, BLAST_FM),
        ('M→F', M_PROT, DB_F, BLAST_MF),
    ]:
        print(f'  blastp: {os.path.basename(query)} vs {os.path.basename(db)} '
              f'({n_cpu} threads) → {os.path.basename(out)}')
        subprocess.run([
            'blastp', '-query', query, '-db', db, '-out', out,
            '-outfmt', '6 qseqid sseqid pident length evalue bitscore',
            '-evalue', str(EVALUE), '-max_target_seqs', str(MAX_HITS),
            '-num_threads', str(n_cpu),
        ], check=True, capture_output=True)

    # ── [3/5] Reciprocal best hits ────────────────────────────────────────
    print('\n[3/5] Finding reciprocal best hits...')

    def best_hits(tsv):
        bh = {}
        with open(tsv) as fh:
            for line in fh:
                p = line.rstrip().split('\t')
                if len(p) < 6: continue
                q, s = p[0], p[1]
                if q not in bh:
                    bh[q] = s
        return bh

    fm_best = best_hits(BLAST_FM)
    mf_best = best_hits(BLAST_MF)
    print(f'  F→M best hits: {len(fm_best):,}')
    print(f'  M→F best hits: {len(mf_best):,}')

    rbh_pairs = []
    for gid_F, gid_M in fm_best.items():
        if mf_best.get(gid_M) == gid_F:
            rbh_pairs.append((gid_F, gid_M))
    print(f'  Reciprocal best hits: {len(rbh_pairs):,}')

    # ── [4/5] Load sequences ──────────────────────────────────────────────
    print('\n[4/5] Loading sequences...')
    f_prot = {r.id: r for r in SeqIO.parse(F_PROT, 'fasta')}
    m_prot = {r.id: r for r in SeqIO.parse(M_PROT, 'fasta')}
    f_cds  = {r.id: r for r in SeqIO.parse(F_CDS,  'fasta')}
    m_cds  = {r.id: r for r in SeqIO.parse(M_CDS,  'fasta')}
    print(f'  F prot={len(f_prot):,}  M prot={len(m_prot):,}  '
          f'F CDS={len(f_cds):,}  M CDS={len(m_cds):,}')

    # ── [5/5] Ka/Ks ───────────────────────────────────────────────────────
    print(f'\n[5/5] Calculating Ka/Ks for {len(rbh_pairs):,} RBH pairs ({n_cpu} workers)...')

    tasks = []
    for i, (gid_F, gid_M) in enumerate(rbh_pairs):
        if gid_F not in f_prot or gid_M not in m_prot:
            continue
        if gid_F not in f_cds or gid_M not in m_cds:
            continue
        tasks.append((i, len(rbh_pairs), gid_F, gid_M,
                      f_prot[gid_F], m_prot[gid_M],
                      f_cds[gid_F],  m_cds[gid_M]))

    results = []
    with multiprocessing.Pool(n_cpu) as pool:
        for j, res in enumerate(pool.imap_unordered(process_pair, tasks, chunksize=8)):
            results.append(res)
            if (j + 1) % 500 == 0 or (j + 1) == len(tasks):
                n_ok = sum(1 for r in results if r['flag'] == 'ok')
                print(f'  {j+1:6d}/{len(tasks)}  ok={n_ok:,}')

    # ── Write output ──────────────────────────────────────────────────────
    print('\nWriting output...')
    header = ['gene_F', 'gene_M', 'scaffold_F', 'scaffold_M',
              'len_F', 'len_M', 'n_codons', 'aln_len',
              'Ks', 'Ka', 'Ka_Ks', 'flag']
    with open(TSV_OUT, 'w') as fh:
        fh.write('\t'.join(header) + '\n')
        for r in results:
            row = [
                r['gene_F'], r['gene_M'],
                r['scaffold_F'], r['scaffold_M'],
                str(r['len_F']), str(r['len_M']),
                '', '',
                f"{r['Ks']:.6f}" if math.isfinite(r['Ks']) else 'NA',
                f"{r['Ka']:.6f}" if math.isfinite(r['Ka']) else 'NA',
                f"{r['Ka_Ks']:.6f}" if math.isfinite(r['Ka_Ks']) else 'NA',
                r['flag'],
            ]
            fh.write('\t'.join(row) + '\n')
    print(f'  Wrote {len(results):,} rows → {TSV_OUT}')

    # Flag summary
    flag_counts = defaultdict(int)
    for r in results:
        flag_counts[r['flag']] += 1
    print('\nFlag summary:')
    for flag, n in sorted(flag_counts.items(), key=lambda x: -x[1]):
        print(f'  {flag:<40} {n:>8,}')

    ok_results = [r for r in results if r['flag'] == 'ok']
    ks_vals = [r['Ks'] for r in ok_results if math.isfinite(r['Ks'])]
    print(f'\nKs statistics (n={len(ks_vals):,} ok pairs):')
    if ks_vals:
        import statistics
        print(f'  median = {statistics.median(ks_vals):.4f}')
        print(f'  mean   = {sum(ks_vals)/len(ks_vals):.4f}')
        print(f'  min    = {min(ks_vals):.4f}')
        print(f'  max    = {max(ks_vals):.4f}')

    # Caption
    import statistics as _stats
    ks_median_str = f'{_stats.median(ks_vals):.4f}' if ks_vals else 'N/A'
    caption = (
        f'Reciprocal-best-BLAST ortholog pairs and synonymous substitution rates (Ks) '
        f'for Phaeomegaceros chiloensis female (PhchiF, 14765-5) and male (PhchiM, 14765-4).\n'
        f'Generated: {datetime.datetime.now().strftime("%Y-%m-%d %H:%M")}\n'
        f'Total RBH pairs: {len(rbh_pairs):,}\n'
        f'Pairs with Ks estimate (flag=ok): {len(ok_results):,}\n'
        f'Ks median (ok pairs): {ks_median_str}\n'
        f'Output: {TSV_OUT}\n'
    )
    with open(CAPTION_OUT, 'w') as fh:
        fh.write(caption)
    print(f'  Wrote caption → {CAPTION_OUT}')

    print(f'\nDone. Results in {OUTDIR}/')


if __name__ == '__main__':
    main()
