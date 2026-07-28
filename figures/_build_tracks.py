"""
Build per-chromosome density tracks for the multi-track ideogram.

Outputs:
    data/tracks/{Hap1,Hap2}.tracks.tsv
    columns: chrom, start, end, gene_count, repeat_bp, repeat_frac

Usage:
    python scripts/_build_tracks.py Hap1
    python scripts/_build_tracks.py Hap2
"""
import sys
import re
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
from _config import DATA_DIR
from _utils import read_fai


BIN_SIZE = 1_000_000   # 1 Mb bins
CHR_LEVEL_RE = r'^chr([1-9]|[12][0-9]|3[0-8]|X|Y)$'


def build_bins(fai_rows, bin_size=BIN_SIZE):
    out = {}
    for name, length in fai_rows:
        n_bins = int(np.ceil(length / bin_size))
        starts = np.arange(n_bins) * bin_size
        ends   = np.minimum(starts + bin_size, length)
        out[name] = pd.DataFrame({'start': starts, 'end': ends})
    return out


def gene_density(gff_path, bins_per_chr):
    print(f"  reading GFF: {gff_path}")
    rows = []
    with open(gff_path) as fh:
        for line in fh:
            if line.startswith('#'):
                continue
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 9 or parts[2] != 'gene':
                continue
            try:
                rows.append((parts[0], int(parts[3])))
            except ValueError:
                continue
    df = pd.DataFrame(rows, columns=['chrom', 'pos'])
    print(f"    {len(df):,} genes loaded")

    out = {}
    for chrom, bin_df in bins_per_chr.items():
        sub = df[df['chrom'] == chrom]
        if len(sub) == 0:
            counts = np.zeros(len(bin_df), dtype=int)
        else:
            bidx = (sub['pos'].values // BIN_SIZE).astype(int)
            bidx = np.clip(bidx, 0, len(bin_df) - 1)
            counts = np.bincount(bidx, minlength=len(bin_df))
        result = bin_df.copy()
        result['gene_count'] = counts
        out[chrom] = result
    return out


def repeat_coverage(rm_out_path, bins_per_chr, bin_size=BIN_SIZE):
    print(f"  reading RepeatMasker .out: {rm_out_path}")
    rows = []
    with open(rm_out_path) as fh:
        for _ in range(3):
            fh.readline()  # skip 3-line header
        for line in fh:
            line = line.strip()
            if not line:
                continue
            parts = line.split()
            if len(parts) < 7:
                continue
            try:
                rows.append((parts[4], int(parts[5]), int(parts[6])))
            except (ValueError, IndexError):
                continue
    print(f"    {len(rows):,} repeat hits loaded")

    out = {}
    for chrom, bin_df in bins_per_chr.items():
        result = bin_df.copy()
        result['repeat_bp'] = 0
        out[chrom] = result

    # vectorize per-chromosome
    rm_df = pd.DataFrame(rows, columns=['chrom', 'begin', 'end'])
    for chrom, bin_df in out.items():
        sub = rm_df[rm_df['chrom'] == chrom]
        if len(sub) == 0:
            continue
        begins = sub['begin'].values
        ends   = sub['end'].values
        b_starts = (begins - 1) // bin_size
        b_ends   = (ends - 1)   // bin_size
        for i in range(len(sub)):
            bs, be = b_starts[i], b_ends[i]
            beg, en = begins[i], ends[i]
            if bs == be:
                if 0 <= bs < len(bin_df):
                    bin_df.loc[bs, 'repeat_bp'] += (en - beg + 1)
            else:
                for b in range(max(0, bs), min(len(bin_df), be + 1)):
                    bl = bin_df.loc[b, 'start'] + 1
                    br = bin_df.loc[b, 'end']
                    sl = max(beg, bl)
                    sr = min(en, br)
                    bin_df.loc[b, 'repeat_bp'] += max(0, sr - sl + 1)

    for chrom in out:
        bin_df = out[chrom]
        sz = bin_df['end'] - bin_df['start']
        bin_df['repeat_frac'] = (bin_df['repeat_bp'] /
                                  sz.replace(0, 1)).clip(upper=1.0)
    return out


def main():
    if len(sys.argv) != 2 or sys.argv[1] not in ('Hap1', 'Hap2'):
        print("Usage: python _build_tracks.py {Hap1|Hap2}")
        sys.exit(1)

    hap = sys.argv[1]
    out_dir = DATA_DIR / 'tracks'
    out_dir.mkdir(exist_ok=True, parents=True)

    print(f"=== {hap} ===")
    fai_path = DATA_DIR / 'assembly_fai' / f'Jindo{hap}.fai'
    gff_path = DATA_DIR / 'braker3' / f'{hap}.braker.gff3'
    rm_path  = DATA_DIR / 'repeatmasker' / f'{hap}.RM.out'

    fai_rows = read_fai(fai_path)
    fai_rows = [(n, L) for (n, L) in fai_rows if re.search(CHR_LEVEL_RE, n)]

    bins = build_bins(fai_rows)
    gene_bins   = gene_density(gff_path, bins)
    repeat_bins = repeat_coverage(rm_path, bins)

    all_chunks = []
    for chrom in sorted(bins.keys()):
        g = gene_bins[chrom]
        r = repeat_bins[chrom]
        merged = g.merge(r[['start', 'end', 'repeat_bp', 'repeat_frac']],
                         on=['start', 'end'])
        merged.insert(0, 'chrom', chrom)
        all_chunks.append(merged)
    full = pd.concat(all_chunks, ignore_index=True)

    out_tsv = out_dir / f'{hap}.tracks.tsv'
    full.to_csv(out_tsv, sep='\t', index=False)
    print(f"  wrote: {out_tsv}  ({len(full):,} bins)")


if __name__ == '__main__':
    main()
