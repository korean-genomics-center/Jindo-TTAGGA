#!/usr/bin/env python3
"""Demo: recover the two Jindo paternal chrY lineages from raw chrY genotypes.

Reproduces the core result of Figure 5a and Supplementary Data 6-8 starting
from the haploid chrY genotype matrix of the 20-male cohort.

Usage:
    python3 demo/run_demo.py                    # all 20 males (main analysis)
    python3 demo/run_demo.py --exclude J73478   # sensitivity analysis

Runtime: about 15 seconds on a normal desktop computer.
"""
import argparse
import numpy as np
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import squareform

GT = "demo/data/chrY_genotypes.tsv"
IDS = "demo/data/sample_ids.txt"
MAX_MISSING = 1
N_BOOT = 100
rng = np.random.default_rng(0)

ap = argparse.ArgumentParser()
ap.add_argument("--exclude", nargs="*", default=[], help="sample IDs to exclude (sensitivity analysis)")
args = ap.parse_args()

ids = np.array([l.strip() for l in open(IDS) if l.strip()])
raw = [l.rstrip("\n").split("\t")[4:] for l in open(GT)]
G = np.array(raw, dtype=object)
assert G.shape[1] == len(ids), f"{G.shape[1]} columns vs {len(ids)} ids"
if args.exclude:
    keep_ind = ~np.isin(ids, args.exclude)
    ids, G = ids[keep_ind], G[:, keep_ind]
    print(f"Excluded: {' '.join(args.exclude)}")
print(f"Loaded {G.shape[0]:,} chrY sites x {len(ids)} males")

miss = (G == ".")
keep = miss.sum(axis=1) <= MAX_MISSING
G, miss = G[keep], miss[keep]
print(f"After missingness filter (<= {MAX_MISSING} missing): {G.shape[0]:,} sites")

n = len(ids)

def dmat(mat, mmask):
    D = np.zeros((n, n))
    for i in range(n):
        for j in range(i + 1, n):
            both = ~mmask[:, i] & ~mmask[:, j]
            D[i, j] = D[j, i] = np.sum(mat[both, i] != mat[both, j])
    return D

def partition(D):
    Z = linkage(squareform(D, checks=False), method="average")
    return fcluster(Z, t=2, criterion="maxclust")

D = dmat(G, miss)
labels = partition(D)
print("\nTwo-lineage partition")
for g in np.unique(labels):
    members = sorted(ids[labels == g])
    print(f"  lineage {g} (n = {len(members)}): {' '.join(members)}")

iu = np.triu_indices(n, k=1)
same = labels[iu[0]] == labels[iu[1]]
w, b = np.median(D[iu][same]), np.median(D[iu][~same])
print(f"\nMedian within-lineage distance : {w:.0f}")
print(f"Median between-lineage distance: {b:.0f}")
print(f"Separation ratio               : {b / w:.1f}")

fixed = 0
for s in range(G.shape[0]):
    alleles = []
    ok = True
    for g in np.unique(labels):
        pick = (labels == g) & ~miss[s]
        a = set(G[s, pick])
        if len(a) != 1:
            ok = False
            break
        alleles.append(a.pop())
    if ok and len(set(alleles)) == 2:
        fixed += 1
print(f"Fixed differences between lineages: {fixed:,}")

ref = {frozenset(ids[labels == g]) for g in np.unique(labels)}
hits = 0
for _ in range(N_BOOT):
    idx = rng.integers(0, G.shape[0], G.shape[0])
    lab = partition(dmat(G[idx], miss[idx]))
    hits += {frozenset(ids[lab == g]) for g in np.unique(lab)} == ref
print(f"\nPartition recovered in {hits}/{N_BOOT} bootstrap replicates")
