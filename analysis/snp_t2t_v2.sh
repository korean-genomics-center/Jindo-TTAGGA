#!/bin/bash

set -e

export PATH=${WORK_DIR}/miniconda3/envs/TTT/bin:$PATH

REF=${JINDO_ROOT}/Results/Repeat_Annotation/RepeatMasker/Hap1/FinalASM.J495799_Child.Trio_Specific.hap1.fasta.masked
BAM_LIST=${JINDO_ROOT}/Analysis/Population_20Jindo/VCF/bam_list_t2t.txt
VCF_DIR=${JINDO_ROOT}/Analysis/Population_20Jindo/VCF
STATS=${JINDO_ROOT}/Analysis/Population_20Jindo/Stats
mkdir -p ${VCF_DIR} ${VCF_DIR}/per_chr_t2t ${STATS}

echo "================================================="
echo "snp_t2t_v2: $(date) on $(hostname)"
echo "================================================="

# BAM list
ls ${JINDO_ROOT}/Analysis/Population_20Jindo/BAM/jindo_*_hap1.sorted.bam | sort > ${BAM_LIST}
echo "BAM count: $(wc -l < ${BAM_LIST})"

# Reference fai
if [ ! -f "${REF}.fai" ]; then
    samtools faidx ${REF}
fi

# Chromosome list (autosome + chrX)
CHR_LIST=${VCF_DIR}/chr_list_t2t.txt
awk '{print $1}' ${REF}.fai | grep -E "^chr[0-9X]+$" | sort -V > ${CHR_LIST}
echo "Chromosomes ($(wc -l < ${CHR_LIST})):"
cat ${CHR_LIST}

# 함수 정의 (직접 호출, subshell 없이)
call_chr() {
    local CHR=$1
    local OUT=${VCF_DIR}/per_chr_t2t/cohort_${CHR}.vcf.gz
    local LOG=${VCF_DIR}/per_chr_t2t/cohort_${CHR}.log
    
    # 이미 있고 정상이면 skip
    if [ -f "${OUT}" ] && bcftools view -h ${OUT} > /dev/null 2>&1; then
        echo "[$(date)] ${CHR}: already done, skipping"
        return 0
    fi
    
    echo "[$(date)] ${CHR}: starting"
    
    # mpileup + call (직접 파일에 쓰기)
    bcftools mpileup \
        -f ${REF} \
        -b ${BAM_LIST} \
        -r ${CHR} \
        -a "FORMAT/AD,FORMAT/DP,FORMAT/SP,INFO/AD" \
        -q 20 -Q 20 \
        --threads 4 2>${LOG} | \
    bcftools call \
        -mv \
        -P 0.001 \
        --threads 4 \
        -Oz -o ${OUT} 2>>${LOG}
    
    if [ ! -s "${OUT}" ]; then
        echo "[$(date)] ${CHR}: FAILED (empty output)"
        return 1
    fi
    
    bcftools index --threads 2 ${OUT}
    
    local NVAR=$(bcftools view -H ${OUT} | wc -l)
    echo "[$(date)] ${CHR}: DONE (${NVAR} variants)"
}

# Background parallel execution with semaphore (max 16 jobs)
echo ""
echo "=== Starting parallel calling (max 16 background jobs × 4t = 64t) ==="

MAX_JOBS=16
JOB_COUNT=0

while read CHR; do
    # 16 jobs 차면 일부 끝날 때까지 대기
    while [ "$(jobs -r | wc -l)" -ge ${MAX_JOBS} ]; do
        sleep 5
    done
    
    # background에서 실행
    call_chr ${CHR} &
    JOB_COUNT=$((JOB_COUNT+1))
    echo "[$(date)] Submitted ${CHR} (job ${JOB_COUNT}, running: $(jobs -r | wc -l))"
done < ${CHR_LIST}

# 모든 background job 끝날 때까지 대기
echo ""
echo "=== Waiting for all background jobs to finish ==="
wait
echo "[$(date)] All chromosomes processed"

# Verification
echo ""
echo "=== Verification ==="
N_DONE=0
N_FAIL=0
FAILED_CHRS=""
for chr in $(cat ${CHR_LIST}); do
    VCF=${VCF_DIR}/per_chr_t2t/cohort_${chr}.vcf.gz
    if [ -f "${VCF}" ] && bcftools view -h ${VCF} > /dev/null 2>&1; then
        N_DONE=$((N_DONE+1))
    else
        N_FAIL=$((N_FAIL+1))
        FAILED_CHRS="${FAILED_CHRS} ${chr}"
    fi
done
echo "Done: ${N_DONE}, Failed: ${N_FAIL}"
if [ ${N_FAIL} -gt 0 ]; then
    echo "FAILED chromosomes: ${FAILED_CHRS}"
    exit 1
fi

# Concat
echo ""
echo "=== Concat ==="
RAW_VCF=${VCF_DIR}/cohort_20jindo.t2t.raw.vcf.gz
bcftools concat --threads 16 -Oz -o ${RAW_VCF} \
    $(cat ${CHR_LIST} | sed "s|^|${VCF_DIR}/per_chr_t2t/cohort_|;s|$|.vcf.gz|")
bcftools index --threads 4 ${RAW_VCF}
echo "Raw: $(bcftools view -H ${RAW_VCF} | wc -l) variants"

# Norm
echo ""
echo "=== Norm ==="
NORM_VCF=${VCF_DIR}/cohort_20jindo.t2t.norm.vcf.gz
bcftools norm -m -any -f ${REF} --threads 16 -Oz -o ${NORM_VCF} ${RAW_VCF}
bcftools index --threads 4 ${NORM_VCF}

# Filter
echo ""
echo "=== Filter ==="
FILT_VCF=${VCF_DIR}/cohort_20jindo.t2t.filtered.vcf.gz
bcftools view --threads 16 -v snps -m 2 -M 2 \
    -e 'QUAL<20 || INFO/DP<100 || F_MISSING>0.1 || MAF<0.05' \
    -Oz -o ${FILT_VCF} ${NORM_VCF}
bcftools index --threads 4 ${FILT_VCF}

# Indels
INDEL_VCF=${VCF_DIR}/cohort_20jindo.t2t.indels.vcf.gz
bcftools view --threads 8 -v indels \
    -e 'QUAL<20 || INFO/DP<100' \
    -Oz -o ${INDEL_VCF} ${NORM_VCF}
bcftools index --threads 4 ${INDEL_VCF}

bcftools stats ${FILT_VCF} > ${STATS}/t2t_filtered_stats.txt

echo ""
echo "=== Summary ==="
echo "Filtered SNPs: $(bcftools view -H ${FILT_VCF} | wc -l)"
echo "Filtered Indels: $(bcftools view -H ${INDEL_VCF} | wc -l)"
ls -la ${VCF_DIR}/*.vcf.gz

echo ""
echo "snp_t2t_v2 Finished: $(date)"
