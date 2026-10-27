#!/usr/bin/env python3
"""Flag chrY sites with mixed allele calls.

chrY is haploid in males, so a male whose reads support two alleles at a site indicates
reads from paralogous or X-linked sequence. A male is counted as mixed at a site when the
minor-allele read fraction exceeds 0.2 at a total allelic depth of at least 4 (Methods).

Input : the chrY call set (VCF/BCF, FORMAT/AD), read with `bcftools query`.
Output: TSV with POS, site MQ and the number of males with a mixed call (mixed_call_males),
        the column reported in Supplementary Data 5.
Usage : python3 chrY_mixed_calls_v1.py chrY_clean.vcf.gz > chrY_mixed_calls.tsv
"""
import subprocess, sys

MIN_DEPTH, MIN_MINOR = 4, 0.2

def mixed(ad):
    v = [int(x) for x in ad.split(',') if x not in ('.', '')]
    s = sum(v)
    return s >= MIN_DEPTH and (s - max(v)) / s > MIN_MINOR

def main(vcf):
    q = subprocess.run(['bcftools', 'query', '-f', '%POS\t%INFO/MQ[\t%AD]\n', vcf],
                       check=True, capture_output=True, text=True).stdout
    print('pos\tsite_MQ\tmixed_call_males')
    n = flagged = 0
    for line in q.splitlines():
        f = line.split('\t')
        k = sum(mixed(ad) for ad in f[2:])
        print(f'{f[0]}\t{f[1]}\t{k}')
        n += 1; flagged += k > 0
    sys.stderr.write(f'sites {n}; sites with no mixed call {n - flagged}\n')

if __name__ == '__main__':
    main(sys.argv[1])
