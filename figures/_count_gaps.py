"""
Count gap stretches (>=10 consecutive Ns) and total Ns in each fasta.

Usage:
    python _count_gaps.py <fasta_path> <output_label>

Output (stdout):
    label\ttotal_N_bp\tgap_count\tcontig_count\n50_contig_bp

Where contig_count = scaffold_count + (gaps that split scaffolds into contigs).
"""
import sys
import re

if len(sys.argv) != 3:
    print("Usage: python _count_gaps.py <fasta> <label>")
    sys.exit(1)

fasta_path = sys.argv[1]
label = sys.argv[2]

# Collect contig lengths (split each scaffold at >=10 N runs)
GAP_PATTERN = re.compile(r'N{10,}', re.IGNORECASE)
total_N = 0
gap_count = 0
contig_lengths = []

with open(fasta_path) as fh:
    cur_seq = []
    for line in fh:
        if line.startswith('>'):
            if cur_seq:
                seq = ''.join(cur_seq).upper()
                # contigs = split at gaps
                contigs = GAP_PATTERN.split(seq)
                for c in contigs:
                    if c:
                        contig_lengths.append(len(c))
                # count gaps
                gaps = GAP_PATTERN.findall(seq)
                gap_count += len(gaps)
                total_N += sum(len(g) for g in gaps)
            cur_seq = []
        else:
            cur_seq.append(line.strip())
    # last
    if cur_seq:
        seq = ''.join(cur_seq).upper()
        contigs = GAP_PATTERN.split(seq)
        for c in contigs:
            if c:
                contig_lengths.append(len(c))
        gaps = GAP_PATTERN.findall(seq)
        gap_count += len(gaps)
        total_N += sum(len(g) for g in gaps)

# Contig N50
sorted_c = sorted(contig_lengths, reverse=True)
total_contig = sum(sorted_c)
half = total_contig / 2
cum = 0
n50_contig = 0
for L in sorted_c:
    cum += L
    if cum >= half:
        n50_contig = L
        break

print(f"{label}\t{total_N}\t{gap_count}\t{len(contig_lengths)}\t{n50_contig}")
