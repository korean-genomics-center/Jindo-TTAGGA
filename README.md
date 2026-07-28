# Jindo-TTAGGA

Analysis and figure code for:

**A complete dog Y chromosome from a gapless telomere-to-telomere Jindo genome
resolves canine paternal lineages**

Choi *et al.* (2026), manuscript under review.

---

## Data availability

Raw sequencing reads, the assembled haplotypes and the population resequencing
data are deposited at NCBI under BioProject **PRJNA1494175**.

This repository contains analysis code only; no sequence data are redistributed here.

---

## Paths

Scripts refer to input locations through two placeholders:

- `${JINDO_ROOT}` — project data root
- `${WORK_DIR}` — local working directory (tool installations, scratch)

Shell scripts read these from the environment. Python scripts contain them as
literal strings and must be edited before running. The code is provided so that
each reported value can be traced to the computation that produced it; it is not
packaged for turnkey re-execution.

---

## Script-to-result mapping

### `analysis/`

| Script | Supports |
|---|---|
| `kin_job.sh` | PLINK IBD on 13,945,113 autosomal SNPs; max PI_HAT 0.246 |
| `c1_test.py` | Within/between-lineage PI_HAT (0.024 / 0.015), permutation *P* = 0.21, Mantel *r* = −0.13 |
| `fixdiff.py` | 3,301 fixed differences between paternal lineages; Supplementary Data 8 |
| `_d.sh`, `_d2.sh` | chrY read depth with and without MAPQ filtering; Supplementary Table 7 |
| `verify_job.sh` | Cross-checks of reported counts against the call sets |
| `snp_t2t_v2.sh` | GATK calling against Jindo1-G-TTAGGA |
| `snp_ros_call.sh` | GATK calling against ROS_Cfam_1.0 (comparison baseline) |
| `snp_chrY_hap2.sh` | Haploid chrY genotyping; 11,404 PASS SNPs |
| `orthofinder_v2_full.sh` | Orthogroup assignment; 2,203 family members / 960 complete-novel |
| `build_db.sh`, `run_tblastn_fast.sh`, `summarize_v2.py` | tblastn classification; 487 sequence-absent / 2,676 annotation-absent |
| `methylation_pbcpg_both.sh` | chrY CpG methylation (pb-CpG-tools, pileup_calling_model v1) |
| `stable3.py` | Per-chromosome QV and centromere intervals; Supplementary Table 3 |
| `sdata11.sh`, `sdata_extra.py`, `fix2.py` | Supplementary Data tables, including the ROS and 4-Mb control distance matrices |
| `restruct.py` | Supplementary Data workbook assembly |

### `figures/`

| Script | Produces |
|---|---|
| `fig1_combined_v18.py` | Figure 1 (panels from `fig1a_v13_final.py`, `fig1b_landscape_v15.py`) |
| `fig2ros_combined_v10.py` | Figure 2 (panels from `fig2ros_a_v4_chr9_inversion.py`, `fig2ros_b_v6_length.py`, `fig2ros_d_v8_benchmark.py`) |
| `fig3_combined_v14.py` | Figure 3 (panels from `fig3a_v7_chrY_map_methyl.py`, `fig3_tspy_v10.py`, `fig3d_v18_canid_landscape.py`) |
| `fig4_refabsent_combined_v19.py` | Figure 4 (panel d from `fig5a_callable_v10.py`; the file name predates the current figure numbering) |
| `fig5_chrY_haplogroup_v5.py` | Figure 5 |
| `Sfig7_completeness_v2.py` | Supplementary Figure 7 |
| `measure_chrY_all.sh` | chrY depth and read-density measurements used in Figure 5 and Supplementary Table 7 |
| `_config.py`, `_utils.py`, `_build_tracks.py`, `_count_gaps.py` | Shared modules imported by the above |

Assembly generation itself (HiFi/ONT assembly, polishing, gap closing) used
published tools with the parameters given in the Methods and is not duplicated here.

---

## Requirements

External tools and versions: `env/tool_versions.txt`.
Python dependencies: `env/requirements.txt`.

---

## Contact

Jong Bhak — jongbhak@genomics.org (corresponding author)
Hyoungjin Choi — ORCID 0009-0007-0504-8870

## Funding

Korea Heritage Service, project no. 2.250780.01.
