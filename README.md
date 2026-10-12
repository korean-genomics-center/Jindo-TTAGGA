# Jindo-TTAGGA

Analysis and figure code for:

**A complete dog Y chromosome from a gapless telomere-to-telomere Jindo genome
resolves canine paternal lineages**

Choi *et al.* (2026), manuscript under review.

Licensed under the MIT License (see `LICENSE`).

---

## Quick start

```bash
git clone https://github.com/korean-genomics-center/Jindo-TTAGGA.git
cd Jindo-TTAGGA
pip install -r env/requirements.txt
python3 demo/run_demo.py
```

The demo reproduces the central result of the paper in about 15 seconds from
data included in this repository. See `demo/README.md`.

---

## Data availability

Raw sequencing reads, the assembled haplotypes and the population resequencing
data are deposited at NCBI under BioProject **PRJNA1494175**.

This repository contains analysis code only; no sequence data are
redistributed here, apart from the 642 kB genotype matrix used by the demo.

---

## System requirements

**Operating system.** Linux (x86_64); developed and run on Ubuntu 22.04.
The demo and the figure scripts also run on macOS. The pipelines in
`analysis/` assume a Linux cluster with a scheduler.

**Non-standard hardware.** None. The demo and the figure scripts run on any
desktop computer. The variant-calling and annotation pipelines in `analysis/`
were run on an SGE cluster (384 slots and 3 TB RAM per node).

**Python.** Python 3.9 or later. Package versions are pinned in
`env/requirements.txt` and the code has been run on those versions:
matplotlib 3.10.1, numpy 1.26.4, scipy 1.16.0, pandas 2.2.3, openpyxl 3.1.5,
cartopy 0.25.0. The demo needs only numpy and scipy.

**External bioinformatics tools.** Names and versions are in
`env/tool_versions.txt`. These are required only for the upstream pipelines in
`analysis/`, not for the demo or the figure scripts.

---

## Installation guide

No compilation or installation of this repository is required; the scripts run
in place.

```bash
git clone https://github.com/korean-genomics-center/Jindo-TTAGGA.git
cd Jindo-TTAGGA
python3 -m venv venv && source venv/bin/activate
pip install -r env/requirements.txt
```

Typical install time on a normal desktop computer: **2-3 minutes**, almost all
of it downloading NumPy, SciPy, pandas, Matplotlib and Cartopy.

The external tools listed in `env/tool_versions.txt` are needed only to re-run
the upstream pipelines from raw reads and are installed separately; most are
available through Bioconda.

---

## Demo

See `demo/README.md` for input format, expected output and run time.

```bash
python3 demo/run_demo.py
```

Expected: 11,404 chrY sites filtered to 7,375, a 10:10 partition of the cohort
into two paternal lineages, a separation ratio of 5.3, 2,525 fixed differences
and 100/100 bootstrap support. About 15 seconds.

---

## Instructions for use

### Paths

Scripts refer to input locations through two placeholders:

- `${JINDO_ROOT}` - project data root
- `${WORK_DIR}` - local working directory (tool installations, scratch)

Shell scripts read these from the environment. Python scripts contain them as
literal strings and must be edited before running.

### Running on your own data

The figure scripts take assembly FASTA, GFF, VCF and BED inputs in standard
formats and can be pointed at another assembly by editing the paths at the top
of each script (`ROOT`, `SCRIPTS_DIR`, `OUT_DIR`). The pipelines in `analysis/`
assume the file layout produced by the upstream steps described in the Methods.

The full pipeline from raw reads is not packaged for turnkey re-execution: it
requires several terabytes of sequencing data and cluster-scale compute. The
code is provided so that each reported value can be traced to the computation
that produced it, and the demo allows the central result to be verified without
any of that.

### Reproduction instructions

Each script maps to specific reported values, as tabulated below.

---

## Script-to-result mapping

### `analysis/`

| Script | Supports |
|---|---|
| `kin_job.sh` | PLINK IBD on 13,945,113 autosomal SNPs; max PI_HAT 0.246 |
| `c1_test.py` | Within/between-lineage PI_HAT (0.024 / 0.015), permutation *P* = 0.21, Mantel *r* = -0.13, *P* = 0.04 |
| `fixdiff.py` | 2,525 fixed differences between paternal lineages across all 20 males; Supplementary Data 8 |
| `_d.sh`, `_d2.sh` | chrY read depth with and without MAPQ filtering; Supplementary Table 7 |
| `verify_job.sh` | Cross-checks of reported counts against the call sets |
| `snp_t2t_v2.sh` | GATK calling against Jindo1-G-TTAGGA |
| `snp_ros_call.sh` | GATK calling against ROS_Cfam_1.0 (comparison baseline) |
| `snp_chrY_hap2.sh` | Haploid chrY genotyping; 11,404 PASS SNPs |
| `orthofinder_v2_full.sh` | Orthogroup assignment (family members / complete-novel) |
| `build_db.sh`, `run_tblastn_fast.sh`, `summarize_v2.py` | tblastn search of reference-absent models against the four dog reference genomes (refined in `revision/`) |
| `methylation_pbcpg_both.sh` | chrY CpG methylation (pb-CpG-tools, pileup_calling_model v1) |
| `stable3.py` | Per-chromosome QV and centromere intervals; Supplementary Table 3 |
| `sdata11.sh`, `sdata_extra.py`, `fix2.py` | Supplementary Data tables, including the ROS and 4-Mb control distance matrices |
| `restruct.py` | Supplementary Data workbook assembly |

### `figures/` - main figures

| Script | Produces |
|---|---|
| `fig1_combined_v18.py` | Figure 1 (panels from `fig1a_v13_final.py`, `fig1b_landscape_v15.py`; the output file is named `fig1_combined_v17` because the save path was not renamed) |
| `fig2ros_combined_v10.py` | Figure 2 (panels from `fig2ros_a_v4_chr9_inversion.py`, `fig2ros_b_v6_length.py`, `fig2ros_d_v8_benchmark.py`) |
| `fig3_combined_v14.py` | Figure 3 (panels from `fig3a_v7_chrY_map_methyl.py`, `fig3_tspy_v10.py`, `fig3d_v18_canid_landscape.py`) |
| `fig4_refabsent_combined_v19.py` | Figure 4 (panel d from `fig5a_callable_v10.py`; the file name predates the current figure numbering) |
| `fig5_chrY_haplogroup_v5.py` | Figure 5 |
| `measure_chrY_all.sh` | chrY depth and read-density measurements used in Figure 5 and Supplementary Table 7 |
| `_config.py`, `_utils.py`, `_build_tracks.py`, `_count_gaps.py` | Shared modules imported by the above |

### `figures/supplementary/` - supplementary figures

Several file names predate the current supplementary figure numbering; the
numbering below is that of the submitted Supplementary Information.

| Script | Produces |
|---|---|
| `Sfig7_completeness_v2.py` | Supplementary Figure 1, assembly completeness and consensus quality |
| `fig3c_v7_chrY_synteny_notitle.py` | Supplementary Figure 2, chrY alignment to ROS_Cfam_1.0 |
| `fig4_funnel_v10.py` | Supplementary Figure 3, reference-absent gene funnel |
| `fig5b_chrY_population_v8.py` | Supplementary Figure 4, chrY population variation |
| `fig5a_v2_phylogeny_notitle.py` | Supplementary Figure 5, canid phylogenetic placement |
| `mito_supp_s7_v6_notitle_22trna.py` | Supplementary Figure 6, mitochondrial genome |
| `fig2ros_c_v6_haps.py` | Supplementary Figure 7, whole-genome Hap1-versus-Hap2 synteny |

Assembly generation itself (HiFi/ONT assembly, polishing, gap closing) used
published tools with the parameters given in the Methods and is not duplicated
here.

---

## Contact

- Jong Bhak - jongbhak@genomics.org (corresponding author)
- Hyoungjin Choi - hjchoi349@unist.ac.kr, ORCID 0009-0007-0504-8870 (analysis code)

## Funding

Korea Heritage Service, project no. 2.250780.01.
