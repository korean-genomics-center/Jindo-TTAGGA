#!/usr/bin/env python3
"""Autosomal relatedness (PLINK PI_HAT) versus chrY lineage, as reported in the manuscript.

Input: the Supplementary Data workbook (sheets 'S-Data 13 autosomal PI_HAT' and 'S-Data 6 chrY distances').

Test 1  Within- versus between-lineage PI_HAT.
        Statistic: median(within-lineage pairs) - median(between-lineage pairs).
        Null: permutation of the lineage labels of the 20 males (10:10); 20,000 random
        permutations (seed 1), as in Methods, and, as a check, exact enumeration of all
        C(20,10) = 184,756 label assignments. Two-sided: |stat| >= |observed|.
Test 2  Mantel test, Pearson correlation between the chrY distance matrix and the PI_HAT matrix,
        20,000 permutations of individual labels (seed 1), two-sided: |r| >= |observed r|.

Usage:  python pihat_lineage_tests_v1.py Jindo1_TTAGGA_Supplementary_Data_1-13.xlsx
Expected output: medians 0.0235 / 0.0153, P = 0.209 (two-sided); Mantel r = -0.13, P ~ 0.04.
Note: analysis/c1_test.py is an earlier exploratory version (difference in means, pair shuffling)
and is superseded by this script.
"""
import sys, itertools
import numpy as np, openpyxl

wb = openpyxl.load_workbook(sys.argv[1], read_only=True)
rows = list(wb['S-Data 13 autosomal PI_HAT'].iter_rows(values_only=True))
hdr = list(rows[1]); iA, lA, iB, lB, ip = (hdr.index(k) for k in ('IID1', 'lineage1', 'IID2', 'lineage2', 'PI_HAT'))
lin = {}; pairs = []
for r in rows[2:]:
    if r[iA] is None: continue
    lin[r[iA]] = r[lA]; lin[r[iB]] = r[lB]; pairs.append((r[iA], r[iB], float(r[ip])))
ids = sorted(lin); n = len(ids); ix = {s: i for i, s in enumerate(ids)}
P = np.zeros((n, n))
for a, b, v in pairs: P[ix[a], ix[b]] = P[ix[b], ix[a]] = v
iu = np.triu_indices(n, 1); pv = P[iu]

def stat(labels):
    same = labels[iu[0]] == labels[iu[1]]
    return np.median(pv[same]) - np.median(pv[~same])

obs_lab = np.array([lin[s] == 'A' for s in ids])
obs = stat(obs_lab)
nA = int(obs_lab.sum()); cnt = tot = 0
for comb in itertools.combinations(range(n), nA):
    lab = np.zeros(n, bool); lab[list(comb)] = True
    tot += 1; cnt += abs(stat(lab)) >= abs(obs) - 1e-12
same = obs_lab[iu[0]] == obs_lab[iu[1]]
print(f'within median {np.median(pv[same]):.4f}  between median {np.median(pv[~same]):.4f}')
rng1 = np.random.default_rng(1); c1 = 0
for _ in range(20000):
    c1 += abs(stat(rng1.permutation(obs_lab))) >= abs(obs) - 1e-12
print(f'Test 1: difference in medians {obs:+.4f}; 20,000 random label permutations (seed 1), two-sided P = {(c1+1)/20001:.3f}; '
      f'exact enumeration of {tot} assignments, two-sided P = {cnt/tot:.4f}')

r2 = list(wb['S-Data 6 chrY distances'].iter_rows(values_only=True))
did = list(r2[1][1:n+1]); D0 = np.array([row[1:n+1] for row in r2[2:2+n]], float)
o = [did.index(s) for s in ids]; D = D0[np.ix_(o, o)]
r0 = np.corrcoef(D[iu], pv)[0, 1]
rng = np.random.default_rng(1); N = 20000; c = 0
for _ in range(N):
    p = rng.permutation(n)
    c += abs(np.corrcoef(D[np.ix_(p, p)][iu], pv)[0, 1]) >= abs(r0)
print(f'Test 2: Mantel (Pearson) r = {r0:.3f}; {N} permutations, two-sided P = {(c+1)/(N+1):.3f}')
