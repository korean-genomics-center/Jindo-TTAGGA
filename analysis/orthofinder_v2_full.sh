#!/bin/bash

set -e
export PATH=${WORK_DIR}/miniconda3/envs/TTT/bin:$PATH

NEW_RUN_DIR=${JINDO_ROOT}/Results/Gene_Annotation/OrthoFinder/Run_v2_full_2026May27
INPUT_DIR=${NEW_RUN_DIR}/input
# Output directory must NOT exist when calling -o; OrthoFinder will create it
OUTPUT_DIR=${NEW_RUN_DIR}/output_v2_full

echo "==================================================="
echo "OrthoFinder v2 Full (Stage 1)"
echo "Started: $(date) on $(hostname)"
echo "==================================================="

orthofinder --version 2>&1 | head -2
echo ""

# Output 디렉토리가 이미 있으면 삭제 (재실행 안전)
if [ -d "${OUTPUT_DIR}" ]; then
    echo "Removing existing output directory: ${OUTPUT_DIR}"
    rm -rf ${OUTPUT_DIR}
fi

echo "=== Input species (18) ==="
ls -la ${INPUT_DIR}/*.fa
echo ""

# OrthoFinder (DIAMOND search, MSA mode)
echo "=== Running OrthoFinder ==="
orthofinder \
    -f ${INPUT_DIR} \
    -t 64 \
    -a 16 \
    -M msa \
    -S diamond \
    -A mafft \
    -T fasttree \
    -o ${OUTPUT_DIR} \
    -n v2_full

echo ""
echo "=== OrthoFinder Finished: $(date) ==="
echo ""

# Result location
RESULT_DIR=$(ls -d ${OUTPUT_DIR}/Results_v2_full* 2>/dev/null | head -1)
echo "Result directory: ${RESULT_DIR}"
ls -la ${RESULT_DIR}/ 2>/dev/null

# Quick stats
if [ -f "${RESULT_DIR}/Comparative_Genomics_Statistics/Statistics_Overall.tsv" ]; then
    echo ""
    echo "=== Statistics Overall ==="
    cat ${RESULT_DIR}/Comparative_Genomics_Statistics/Statistics_Overall.tsv | head -50
fi

if [ -f "${RESULT_DIR}/Orthogroups/Orthogroups.GeneCount.tsv" ]; then
    echo ""
    echo "=== Orthogroup counts ==="
    head -1 ${RESULT_DIR}/Orthogroups/Orthogroups.GeneCount.tsv
    echo "Total: $(tail -n +2 ${RESULT_DIR}/Orthogroups/Orthogroups.GeneCount.tsv | wc -l)"
fi

echo ""
echo "OrthoFinder v2 Full DONE: $(date)"
