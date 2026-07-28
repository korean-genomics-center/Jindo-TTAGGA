#!/bin/bash
set -eu
source ${WORK_DIR}/miniconda3/etc/profile.d/conda.sh
set +u; conda activate TTT; set -u

REF=${JINDO_ROOT}/Resources/Reference
OUT=${JINDO_ROOT}/Results/Gene_Annotation/RefAbsent_Functional/tblastn_genome

GENOMES=(
 "ros:$REF/GCF_014441545.1/GCF_014441545.1_ROS_Cfam_1.0_genomic.fna"
 "gsd:$REF/GCF_011100685.1/GCF_011100685.1_UU_Cfam_GSD_1.0_genomic.fna"
 "boxer:$REF/GCF_000002285.5/GCF_000002285.5_Dog10K_Boxer_Tasha_genomic.fna"
 "zoey:$REF/GCF_005444595.1/GCF_005444595.1_UMICH_Zoey_3.1_genomic.fna"
)
IFS=':' read -r TAG FA <<< "${GENOMES[$((SGE_TASK_ID-1))]}"
echo "=== $TAG : $FA ==="
makeblastdb -in "$FA" -dbtype nucl -out $OUT/db/$TAG -title $TAG
echo "DONE $TAG"
