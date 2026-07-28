"""
Shared utility functions for manuscript figure generation.
"""
import re
from pathlib import Path
from typing import Optional


def parse_busco_short_summary(path) -> dict:
    """
    Parse a BUSCO short_summary.*.txt file.

    Returns a dict with keys:
        C, S, D, F, M  (percentages, float)
        n              (total BUSCOs, int)
        E              (internal-stop fraction, float; None if absent)
        lineage        (e.g. 'carnivora_odb12')
        mode           ('genome' or 'proteins')
        n_complete, n_single, n_dup, n_frag, n_missing  (raw counts, int)
    """
    txt = Path(path).read_text()

    out = {}

    # Lineage
    m = re.search(r'The lineage dataset is:\s*(\S+)', txt)
    out['lineage'] = m.group(1) if m else None

    # Mode
    m = re.search(r'BUSCO was run in mode:\s*(\S+)', txt)
    raw_mode = m.group(1) if m else ''
    if 'genome' in raw_mode:
        out['mode'] = 'genome'
    elif 'protein' in raw_mode:
        out['mode'] = 'proteins'
    else:
        out['mode'] = raw_mode or None

    # Percent line — supports both with and without E (internal stops)
    # e.g.  C:99.3%[S:98.3%,D:1.0%],F:0.2%,M:0.5%,n:13727,E:2.8%
    # e.g.  C:75.1%[S:62.2%,D:12.9%],F:1.9%,M:23.0%,n:9226
    m = re.search(
        r'C:([\d.]+)%\[S:([\d.]+)%,D:([\d.]+)%\],F:([\d.]+)%,M:([\d.]+)%,n:(\d+)(?:,E:([\d.]+)%)?',
        txt
    )
    if not m:
        raise ValueError(f"Cannot parse percentage line in: {path}")
    out['C'] = float(m.group(1))
    out['S'] = float(m.group(2))
    out['D'] = float(m.group(3))
    out['F'] = float(m.group(4))
    out['M'] = float(m.group(5))
    out['n'] = int(m.group(6))
    out['E'] = float(m.group(7)) if m.group(7) else None

    # Raw counts
    counts = {
        'n_complete': r'(\d+)\s+Complete BUSCOs \(C\)',
        'n_single':   r'(\d+)\s+Complete and single-copy BUSCOs \(S\)',
        'n_dup':      r'(\d+)\s+Complete and duplicated BUSCOs \(D\)',
        'n_frag':     r'(\d+)\s+Fragmented BUSCOs \(F\)',
        'n_missing':  r'(\d+)\s+Missing BUSCOs \(M\)',
    }
    for k, pat in counts.items():
        mm = re.search(pat, txt)
        out[k] = int(mm.group(1)) if mm else None

    return out


def parse_compleasm_summary(path) -> dict:
    """
    Parse a Compleasm summary.txt file.
    Format:
        ## lineage: carnivora_odb12
        S:98.76%, 13557
        D:0.74%, 102
        F:0.16%, 22
        I:0.00%, 0
        M:0.34%, 46
        N:13727
    """
    out = {'S': 0, 'D': 0, 'F': 0, 'I': 0, 'M': 0, 'N': 0, 'lineage': None}
    for line in Path(path).read_text().splitlines():
        line = line.strip()
        if line.startswith('## lineage'):
            out['lineage'] = line.split(':', 1)[1].strip()
        elif ':' in line:
            key, rest = line.split(':', 1)
            key = key.strip()
            if key in ('S', 'D', 'F', 'I', 'M'):
                pct = rest.split('%')[0].strip()
                try:
                    out[key] = float(pct)
                except ValueError:
                    pass
            elif key == 'N':
                try:
                    out['N'] = int(rest.split(',')[0].strip())
                except ValueError:
                    pass
    out['C'] = out['S'] + out['D']  # complete = single + duplicated
    return out


def read_fai(path) -> list:
    """
    Read a samtools faidx .fai file.
    Returns list of (name, length) tuples.
    """
    rows = []
    for line in Path(path).read_text().splitlines():
        parts = line.split('\t')
        if len(parts) >= 2:
            rows.append((parts[0], int(parts[1])))
    return rows


def save_fig(fig, path, dpi: int = 600, also_png: bool = True):
    """
    Save a matplotlib figure as PDF (vector) and optionally PNG (raster preview).
    Auto-creates parent directory.
    """
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=dpi, bbox_inches='tight', pad_inches=0.05)
    print(f"  saved: {path}")
    if also_png and path.suffix.lower() == '.pdf':
        png_path = path.with_suffix('.png')
        fig.savefig(png_path, dpi=dpi, bbox_inches='tight', pad_inches=0.05)
        print(f"  saved: {png_path}")


def telomere_status_per_chromosome(
    tidk_tsv,
    fai_path,
    window_size: int = 5000,
    repeat_threshold: int = 100,
):
    """
    Determine telomere presence at the 5' and 3' ends of every chromosome.

    Reads a tidk `_telomeric_repeat_windows.tsv` (columns: id, window,
    forward_repeat_number, reverse_repeat_number, telomeric_repeat) and the
    matching .fai (for chromosome lengths), and returns:

        { chrom_name: {
              'length':     int,
              'left_count': int,    # repeat count summed in 5' window
              'right_count':int,    # repeat count summed in 3' window
              'left':       bool,   # True if left_count >= threshold
              'right':      bool,
          },
          ... }

    Definitions follow the manuscript: telomere "present" at an end if at
    least `repeat_threshold` (forward+reverse) repeats fall within
    `window_size` bp of that end.
    """
    # chromosome lengths
    lengths = {name: L for name, L in read_fai(fai_path)}

    # init
    out = {n: {'length': L, 'left_count': 0, 'right_count': 0,
               'left': False, 'right': False}
           for n, L in lengths.items()}

    # stream the tsv (large file, avoid loading all)
    with open(tidk_tsv) as fh:
        header = fh.readline().rstrip('\n').split('\t')
        # Expect: id  window  forward_repeat_number  reverse_repeat_number  telomeric_repeat
        try:
            i_chr = header.index('id')
            i_win = header.index('window')
            i_fwd = header.index('forward_repeat_number')
            i_rev = header.index('reverse_repeat_number')
        except ValueError:
            raise ValueError(f"Unexpected tidk header: {header}")

        for line in fh:
            parts = line.rstrip('\n').split('\t')
            if len(parts) < max(i_chr, i_win, i_fwd, i_rev) + 1:
                continue
            chrom = parts[i_chr]
            if chrom not in out:
                continue
            try:
                win = int(parts[i_win])
                fwd = int(parts[i_fwd])
                rev = int(parts[i_rev])
            except ValueError:
                continue
            total = fwd + rev
            chrlen = out[chrom]['length']

            # tidk window field is the END coordinate of the bin.
            # The bin spans (win - bin_size, win].
            # We don't know bin_size from the file alone, but the assumption
            # "within `window_size` bp of either end" is satisfied if
            # win <= window_size  (left end)  or  win >= chrlen - window_size  (right end).
            if win <= window_size:
                out[chrom]['left_count'] += total
            if win >= chrlen - window_size:
                out[chrom]['right_count'] += total

    # apply threshold
    for c in out.values():
        c['left']  = c['left_count']  >= repeat_threshold
        c['right'] = c['right_count'] >= repeat_threshold

    return out


def telomere_status_per_chromosome(
    tidk_tsv,
    fai_path,
    window_size: int = 5000,
    repeat_threshold: int = 100,
):
    """
    Determine telomere presence at the 5' and 3' ends of every chromosome.

    Reads a tidk `_telomeric_repeat_windows.tsv` (columns: id, window,
    forward_repeat_number, reverse_repeat_number, telomeric_repeat) and the
    matching .fai (for chromosome lengths), and returns:

        { chrom_name: {
              'length':     int,
              'left_count': int,    # repeat count summed in 5' window
              'right_count':int,    # repeat count summed in 3' window
              'left':       bool,   # True if left_count >= threshold
              'right':      bool,
          },
          ... }

    Definitions follow the manuscript: telomere "present" at an end if at
    least `repeat_threshold` (forward+reverse) repeats fall within
    `window_size` bp of that end.
    """
    # chromosome lengths
    lengths = {name: L for name, L in read_fai(fai_path)}

    # init
    out = {n: {'length': L, 'left_count': 0, 'right_count': 0,
               'left': False, 'right': False}
           for n, L in lengths.items()}

    # stream the tsv (large file, avoid loading all)
    with open(tidk_tsv) as fh:
        header = fh.readline().rstrip('\n').split('\t')
        # Expect: id  window  forward_repeat_number  reverse_repeat_number  telomeric_repeat
        try:
            i_chr = header.index('id')
            i_win = header.index('window')
            i_fwd = header.index('forward_repeat_number')
            i_rev = header.index('reverse_repeat_number')
        except ValueError:
            raise ValueError(f"Unexpected tidk header: {header}")

        for line in fh:
            parts = line.rstrip('\n').split('\t')
            if len(parts) < max(i_chr, i_win, i_fwd, i_rev) + 1:
                continue
            chrom = parts[i_chr]
            if chrom not in out:
                continue
            try:
                win = int(parts[i_win])
                fwd = int(parts[i_fwd])
                rev = int(parts[i_rev])
            except ValueError:
                continue
            total = fwd + rev
            chrlen = out[chrom]['length']

            # tidk window field is the END coordinate of the bin.
            # The bin spans (win - bin_size, win].
            # We don't know bin_size from the file alone, but the assumption
            # "within `window_size` bp of either end" is satisfied if
            # win <= window_size  (left end)  or  win >= chrlen - window_size  (right end).
            if win <= window_size:
                out[chrom]['left_count'] += total
            if win >= chrlen - window_size:
                out[chrom]['right_count'] += total

    # apply threshold
    for c in out.values():
        c['left']  = c['left_count']  >= repeat_threshold
        c['right'] = c['right_count'] >= repeat_threshold

    return out
