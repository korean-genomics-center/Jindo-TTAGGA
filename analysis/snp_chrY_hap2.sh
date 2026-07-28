#!/bin/bash

set -e
set -x

export PATH=${WORK_DIR}/miniconda3/envs/TTT/bin:$PATH

REF=${JINDO_ROOT}/Results/Repeat_Annotation/RepeatMasker/Hap2/FinalASM.J495799_Child.Trio_Specific.hap2.fasta.masked
BAM_DIR=${JINDO_ROOT}/Analysis/Population_20Jindo/BAM_Hap2
VCF_DIR=${JINDO_ROOT}/Analysis/Population_20Jindo/VCF_chrY
STATS=${JINDO_ROOT}/Analysis/Population_20Jindo/Stats_chrY
mkdir -p ${VCF_DIR} ${STATS}

echo "==================================================="
echo "chrY SNP calling (haploid mode, 20 male Jindo)"
echo "Started: $(date) on $(hostname)"
echo "==================================================="

# BAM list 생성 (Hap2 BAM 20개)
BAM_LIST=${VCF_DIR}/bam_list_chrY.txt
ls ${BAM_DIR}/jindo_*_hap2.sorted.bam | sort > ${BAM_LIST}
N_BAMS=$(wc -l < ${BAM_LIST})
echo "Total BAMs: ${N_BAMS} / 20 expected"
cat ${BAM_LIST}

if [ ${N_BAMS} -ne 20 ]; then
    echo "ERROR: Expected 20 BAMs, got ${N_BAMS}"
    exit 1
fi

# chrY SNP calling (haploid mode, 모든 sample이 male)
echo ""
echo "=== bcftools mpileup + call (chrY only, haploid) ==="
RAW_VCF=${VCF_DIR}/cohort_20jindo.chrY.raw.vcf.gz

bcftools mpileup \
    -f ${REF} \
    -b ${BAM_LIST} \
    -r chrY \
    -a "FORMAT/AD,FORMAT/DP,FORMAT/SP,INFO/AD" \
    -q 20 -Q 20 \
    --threads 16 | \
bcftools call \
    -mv \
    -P 0.001 \
    --ploidy 1 \
    --threads 16 \
    -Oz -o ${RAW_VCF}

bcftools index --threads 8 ${RAW_VCF}
echo "Raw chrY variants: $(bcftools view -H ${RAW_VCF} | wc -l)"

# Normalize
NORM_VCF=${VCF_DIR}/cohort_20jindo.chrY.norm.vcf.gz
bcftools norm \
    -m -any \
    -f ${REF} \
    --threads 16 \
    -Oz -o ${NORM_VCF} \
    ${RAW_VCF}
bcftools index --threads 8 ${NORM_VCF}

# Filter biallelic SNPs (haploid)
FILT_VCF=${VCF_DIR}/cohort_20jindo.chrY.filtered.vcf.gz
bcftools view \
    --threads 16 \
    -v snps -m 2 -M 2 \
    -e 'QUAL<20 || INFO/DP<50 || F_MISSING>0.2' \
    -Oz -o ${FILT_VCF} \
    ${NORM_VCF}
bcftools index --threads 8 ${FILT_VCF}

# Indels
INDEL_VCF=${VCF_DIR}/cohort_20jindo.chrY.indels.vcf.gz
bcftools view \
    --threads 8 \
    -v indels \
    -e 'QUAL<20 || INFO/DP<50' \
    -Oz -o ${INDEL_VCF} \
    ${NORM_VCF}
bcftools index --threads 8 ${INDEL_VCF}

# Stats
bcftools stats ${FILT_VCF} > ${STATS}/chrY_filtered_stats.txt
bcftools stats ${INDEL_VCF} > ${STATS}/chrY_indels_stats.txt

# Summary
echo ""
echo "==================================================="
echo "chrY SNP calling COMPLETE: $(date)"
echo "==================================================="
echo "Raw chrY variants: $(bcftools view -H ${RAW_VCF} | wc -l)"
echo "Filtered chrY SNPs: $(bcftools view -H ${FILT_VCF} | wc -l)"
echo "chrY Indels: $(bcftools view -H ${INDEL_VCF} | wc -l)"
echo ""
ls -la ${VCF_DIR}/*.vcf.gz
