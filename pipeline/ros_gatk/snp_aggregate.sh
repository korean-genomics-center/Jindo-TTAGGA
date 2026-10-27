#!/bin/bash
set -e
conda activate TTT
BASE=${JINDO_ROOT}/Analysis/GATK_Recall
cd ${BASE}
VJ=${BASE}/VCF_joint
OUT=${BASE}/VCF_merged
mkdir -p ${OUT}

# .fai 순서대로 상염색체 concat 리스트 만들기 (chr 순서 보장)
make_list () {
  local dir=$1; local out=$2
  for c in $(cut -f1 ${BASE}/intervals/${3} ); do
    echo "${dir}/${c}.vcf.gz"
  done > ${out}
}

# 상염색체 concat (염색체 순서대로)
echo "[$(date)] concat hap1_auto"
make_list ${VJ}/perchr_hap1_auto ${OUT}/list_hap1_auto.txt autosomes_hap1.list
bcftools concat -f ${OUT}/list_hap1_auto.txt -Oz -o ${OUT}/hap1_auto.vcf.gz --threads 8
bcftools index -t ${OUT}/hap1_auto.vcf.gz

echo "[$(date)] concat ros_auto"
make_list ${VJ}/perchr_ros_auto ${OUT}/list_ros_auto.txt autosomes_ros.list
bcftools concat -f ${OUT}/list_ros_auto.txt -Oz -o ${OUT}/ros_auto.vcf.gz --threads 8
bcftools index -t ${OUT}/ros_auto.vcf.gz

# T2T 전체 = hap1_auto + hap1_chrX + hap2_chrY
echo "[$(date)] T2T 통합 (auto + chrX + chrY)"
bcftools concat --allow-overlaps \
  ${OUT}/hap1_auto.vcf.gz ${VJ}/cohort_hap1_chrX.vcf.gz ${VJ}/cohort_hap2_chrY.vcf.gz \
  -Oz -o ${OUT}/T2T_all.vcf.gz --threads 8 2>/dev/null || \
bcftools concat \
  ${OUT}/hap1_auto.vcf.gz ${VJ}/cohort_hap1_chrX.vcf.gz ${VJ}/cohort_hap2_chrY.vcf.gz \
  -Oz -o ${OUT}/T2T_all.vcf.gz --threads 8
bcftools index -t ${OUT}/T2T_all.vcf.gz

# ROS 전체 = ros_auto + ros_chrX + ros_chrY
echo "[$(date)] ROS 통합 (auto + chrX + chrY)"
bcftools concat \
  ${OUT}/ros_auto.vcf.gz ${VJ}/cohort_ros_chrX.vcf.gz ${VJ}/cohort_ros_chrY.vcf.gz \
  -Oz -o ${OUT}/ROS_all.vcf.gz --threads 8
bcftools index -t ${OUT}/ROS_all.vcf.gz

# SNP만 추출 + 카운트
echo "[$(date)] SNP 카운트"
for set in T2T ROS; do
  TOTAL=$(bcftools view -v snps ${OUT}/${set}_all.vcf.gz | grep -vc '^#')
  echo "${set} 전체 SNP (필터 전): ${TOTAL}"
done
echo "[$(date)] DONE"
