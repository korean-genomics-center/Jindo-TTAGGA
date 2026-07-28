#!/bin/bash

set -e
set -x

export PATH=${WORK_DIR}/miniconda3/envs/TTT/bin:$PATH

REF=${JINDO_ROOT}/Analysis/Population_20Jindo/ROS_Reference/ROS_Cfam_1.0.renamed.fa
BAM_DIR=${JINDO_ROOT}/Analysis/Population_20Jindo/BAM_ROS
VCF_DIR=${JINDO_ROOT}/Analysis/Population_20Jindo/VCF_ROS
STATS=${JINDO_ROOT}/Analysis/Population_20Jindo/Stats_ROS
mkdir -p ${VCF_DIR} ${STATS}

echo "==================================================="
echo "ROS cohort SNP calling | Started: $(date) on $(hostname)"
echo "==================================================="

# Step 0: Reference fai
if [ ! -f "${REF}.fai" ]; then
    samtools faidx ${REF}
fi

# Step 1: BAM list
BAM_LIST=${VCF_DIR}/bam_list_ros.txt
ls ${BAM_DIR}/jindo_*_ros.sorted.bam | sort > ${BAM_LIST}
N_BAMS=$(wc -l < ${BAM_LIST})
echo "Total BAMs: ${N_BAMS}"

# Step 2: chromosome list (ROS has chr1, chr2, ... format from .renamed.fa)
CHR_LIST=${VCF_DIR}/chr_list_ros.txt
awk '{print $1}' ${REF}.fai | grep -E "^chr[0-9XY]+$" | sort -V > ${CHR_LIST}
echo "Chromosomes to call:"
cat ${CHR_LIST}

# Step 3: bcftools mpileup + call per chromosome
echo ""
echo "=== Step 3: bcftools mpileup + call ==="
mkdir -p ${VCF_DIR}/per_chr_ros
cat ${CHR_LIST} | xargs -n 1 -P 16 -I {} bash -c '
    CHR=$1
    REF='${REF}'
    BAM_LIST='${BAM_LIST}'
    OUT='${VCF_DIR}'/per_chr_ros/cohort_${CHR}.vcf.gz
    
    if [ ! -f ${OUT} ]; then
        echo "[$(date)] Starting ${CHR}"
        bcftools mpileup \
            -f ${REF} \
            -b ${BAM_LIST} \
            -r ${CHR} \
            -a "FORMAT/AD,FORMAT/DP,FORMAT/SP,INFO/AD" \
            -q 20 -Q 20 \
            --threads 4 2>/dev/null | \
        bcftools call \
            -mv \
            -P 0.001 \
            --threads 4 \
            -Oz -o ${OUT} 2>/dev/null
        
        bcftools index --threads 2 ${OUT}
        echo "[$(date)] ${CHR} done"
    else
        echo "${CHR} already done, skipping"
    fi
' _ {}

# Step 4: Concat
RAW_VCF=${VCF_DIR}/cohort_20jindo.ros.raw.vcf.gz
bcftools concat \
    --threads 16 \
    -Oz -o ${RAW_VCF} \
    $(cat ${CHR_LIST} | sed "s|^|${VCF_DIR}/per_chr_ros/cohort_|;s|$|.vcf.gz|")
bcftools index --threads 4 ${RAW_VCF}

# Step 5: Norm
NORM_VCF=${VCF_DIR}/cohort_20jindo.ros.norm.vcf.gz
bcftools norm \
    -m -any \
    -f ${REF} \
    --threads 16 \
    -Oz -o ${NORM_VCF} \
    ${RAW_VCF}
bcftools index --threads 4 ${NORM_VCF}

# Step 6: Filter
FILT_VCF=${VCF_DIR}/cohort_20jindo.ros.filtered.vcf.gz
bcftools view \
    --threads 16 \
    -v snps -m 2 -M 2 \
    -e 'QUAL<20 || INFO/DP<100 || F_MISSING>0.1 || MAF<0.05' \
    -Oz -o ${FILT_VCF} \
    ${NORM_VCF}
bcftools index --threads 4 ${FILT_VCF}

# Step 7: Indels
INDEL_VCF=${VCF_DIR}/cohort_20jindo.ros.indels.vcf.gz
bcftools view \
    --threads 8 \
    -v indels \
    -e 'QUAL<20 || INFO/DP<100' \
    -Oz -o ${INDEL_VCF} \
    ${NORM_VCF}
bcftools index --threads 4 ${INDEL_VCF}

# Step 8: Stats
bcftools stats ${FILT_VCF} > ${STATS}/ros_filtered_stats.txt
bcftools stats ${INDEL_VCF} > ${STATS}/ros_indels_stats.txt

echo ""
echo "=== Summary ==="
echo "Filtered SNPs: $(bcftools view -H ${FILT_VCF} | wc -l)"
echo "Filtered Indels: $(bcftools view -H ${INDEL_VCF} | wc -l)"

ls -la ${VCF_DIR}/

echo ""
echo "ROS cohort SNP calling | Finished: $(date)"
