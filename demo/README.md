# Demo

Recovers the two Korean Jindo paternal chrY lineages, the central result of
the study, from the haploid chrY genotype matrix of the 20-male cohort.

## Run

```bash
python3 run_demo.py        # from the demo/ directory
```

or, from the repository root:

```bash
python3 demo/run_demo.py
```

Requires only NumPy and SciPy (`pip install -r ../env/requirements.txt`).

## Input

`data/chrY_genotypes.tsv` (642 kB) - 11,404 chrY PASS SNP sites x 20 males.
Tab-separated, no header: chromosome, position, reference allele, alternative
allele, then one column per individual in the order given in
`data/sample_ids.txt`. Genotypes are `0` (reference), `1`-`3` (alternative
alleles) and `.` (missing). Because chrY is haploid in males, each individual
carries a single called allele per site.

This is the unfiltered call set, so the demo reproduces the missingness
filtering step as well as the downstream analysis.

## Expected output

```
Loaded 11,404 chrY sites x 20 males
After missingness filter (<= 1 missing): 7,375 sites

Two-lineage partition
  lineage 1 (n = 10): J12952315 J13044625 J288508 J352388 J496128 J498404 J61783 J62428 J63253 J76290
  lineage 2 (n = 10): J288326 J433555 J495439 J543672 J61435 J62957 J63542 J73201 J73478 J90474

Median within-lineage distance : 832
Median between-lineage distance: 4425
Separation ratio               : 5.3
Fixed differences between lineages: 2,525

Partition recovered in 100/100 bootstrap replicates
```

Every printed value appears in the manuscript: 7,375 filtered sites, the
balanced 10:10 partition, the 5.3-fold separation ratio, 2,525 fixed
differences (Supplementary Data 8) and full bootstrap support.

The sensitivity analysis that excludes J73478 (Supplementary Table 11) is
reproduced with

```
python3 demo/run_demo.py --exclude J73478
```

which gives 7,460 sites, a 10:9 partition with identical membership, a
separation ratio of 5.9, 3,310 fixed differences and 100/100 bootstrap
support.

## Expected run time

About 15 seconds on a normal desktop computer, of which most is the 100
bootstrap replicates.

## Correspondence with the Methods

| Step | Methods | Demo |
|---|---|---|
| Missingness filter | `bcftools view -i 'F_MISSING<0.1'` | at most 1 missing genotype of 20, equivalent over this cohort |
| Pairwise distance | sites with a genotype in both individuals | same |
| Clustering | average-linkage on the distance matrix | `scipy.cluster.hierarchy` |
| Support | 100 bootstrap replicates over sites | same, seeded for reproducibility |
| Individuals | all 20 males retained; J73478 excluded only in the sensitivity analysis | `--exclude J73478` |

The demo begins from called genotypes. Read alignment and GATK variant calling
that produce this matrix are in `analysis/` and require the raw sequencing data
deposited under BioProject PRJNA1494175.
