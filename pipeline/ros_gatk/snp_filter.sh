#!/bin/bash
set -e
conda activate gatk4
BASE=${JINDO_ROOT}/Analysis/GATK_Recall
OUT=${BASE}/VCF_merged
cd ${OUT}

for set in T2T ROS; do
  REF=${BASE}/ref/$( [ "$set" = "T2T" ] && echo hap1 || echo ros ).fasta
  echo "[$(date)] ${set}: SNP 추출"
  gatk SelectVariants -V ${set}_all.vcf.gz --select-type-to-include SNP \
    -O ${set}_snp.vcf.gz 2>/dev/null

  echo "[$(date)] ${set}: 하드필터 표시"
  gatk VariantFiltration -V ${set}_snp.vcf.gz \
    --filter-expression "QD < 2.0" --filter-name "QD2" \
    --filter-expression "FS > 60.0" --filter-name "FS60" \
    --filter-expression "MQ < 40.0" --filter-name "MQ40" \
    --filter-expression "QUAL < 30.0" --filter-name "QUAL30" \
    -O ${set}_snp_flagged.vcf.gz 2>/dev/null

  RAW=$(zcat ${set}_snp.vcf.gz | grep -vc '^#')
  PASS=$(zcat ${set}_snp_flagged.vcf.gz | grep -v '^#' | awk '$7=="PASS"' | wc -l)
  echo "=== ${set}: raw SNP ${RAW}, PASS ${PASS} ==="
done
echo "[$(date)] DONE"
