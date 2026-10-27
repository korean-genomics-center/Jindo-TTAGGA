# Revision analyses (NCOMMS-26-061979, revision 1)

Scripts added in response to the reviewers. Paths inside the scripts point to the
project server and document where each input was read; the analyses they run are
described in the Methods and in Supplementary Tables 8-11.

## chrY lineage sensitivity (Supplementary Table 11)

| Script | Supports |
|---|---|
| `chrY_numbers_v1.py` | All chrY lineage statistics for the 20-male main analysis and the J73478-excluded sensitivity analysis; length-matched 4-Mb control; ROS_Cfam_1.0 comparison with bootstrap; spanning-deletion sensitivity |
| `ros_cluster_v1.py` | Unsupervised clustering of genotypes called against ROS_Cfam_1.0 (same partition as against TTAGGA-P) |
| `single_copy_v2.py`, `run_single_copy_v2.sh` | chrY self-alignment, segmental-duplication intervals (>=10 kb, >=95% identity) and the analysis restricted to sites outside them and the centromere |
| `single_copy_boot_v1.py` | Bootstrap support for the segmental-duplication-excluded analysis |

## Assembly validation (Supplementary Table 8)

| Script | Supports |
|---|---|
| `run_final_syri_v3.sh` | SyRI against ROS_Cfam_1.0 on the deposited sequences (Supplementary Data 4) |
| `run_depth_v3.sh`, `summ_depth_v3.py` | 1-Mb long-read depth uniformity |
| `run_bindepth_v1.sh` | Depth of haplotype-specific and unassigned reads in windows of elevated depth |
| `run_final_checks_v1.sh` | N count, telomere data, Merqury completeness, rDNA array depth |
| `provenance_v2.py` | Sequence-level check that analysis inputs match the deposited assembly |

## Reference-absent gene catalogue (Supplementary Data 1; Supplementary Tables 9-10)

| Script | Supports |
|---|---|
| `F2_characterize_v1.py`, `F2_numbers_v1.py` | Refined catalogue (2,371 models; 1,658 loci): length, family class, domains, expression, haplotype pairing |
| `summarize_reverse_2x2_v1.py` | Reciprocal comparison (Supplementary Table 10) |

## Figures

| Script | Figure |
|---|---|
| `fig2ros_a_v5_chr9_inversion.py`, `fig2ros_combined_v11.py` | Figure 2 (chr9 inversion read from SyRI output) |
| `fig4_domain_counts_v1.py`, `fig4_refabsent_combined_v20.py` | Figure 4 (domain categories computed from InterProScan) |
| `fig4_funnel_v12.py` | Supplementary Figure 3 |

## Pseudoautosomal region and mixed allele calls (added after internal review)

These three scripts contain no server paths; inputs are passed as arguments or through `JINDO_ROOT`.
Requirements: samtools, minimap2, bcftools, Python 3 with openpyxl and scipy.

| Script | Supports |
|---|---|
| `chrY_par_alignment_v1.sh` | Aligns TTAGGA-P chrY to TTAGGA-M chrX to locate the pseudoautosomal region (PAR, chrY 14.70–21.26 Mb) and checks that ROS_Cfam_1.0 carries the PAR on chrX |
| `chrY_mixed_calls_v1.py` | Counts, per chrY site, the males with a mixed allele call (minor-allele read fraction > 0.2 at allelic depth ≥ 4) from FORMAT/AD; Supplementary Data 5 column `mixed_call_males`. On the chrY call set it reports `sites 7375; sites with no mixed call 4465` |
| `chrY_lineage_sensitivity_v1.py` | Recomputes the Supplementary Table 11 rows for all sites, sites without mixed calls, the male-specific region and the PAR from the Supplementary Data workbook |

Expected output of `chrY_lineage_sensitivity_v1.py`:

```
site set                                      sites  fixed   within  between    ratio  boot   pedigree pairs
All sites (main analysis)                      7375   2525      832     4425      5.3   100/100  435-1098
No mixed allele calls                          4465   2423      144     3294     23.0   100/100  59-172
No mixed allele calls, male-specific region    3759   2419        6     3152    525.3   100/100  0-6
Male-specific region (<14.70 Mb)               4763   2518      195     3652     18.7   100/100  14-385
Pseudoautosomal region (>=14.70 Mb)            2612      7      706      774      1.1     0/100  421-770
```
