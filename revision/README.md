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
