#!/bin/bash
set -e
conda activate gatk4
# 인자: $1=tag $2=reftag $3=chr
BASE=${JINDO_ROOT}/Analysis/GATK_Recall
TAG=$1; REFTAG=$2; CHR=$3
REF=${BASE}/ref/${REFTAG}.fasta
GVCFDIR=${BASE}/gvcf_${TAG}
OUTDIR=${BASE}/VCF_joint/perchr_${TAG}
TMP=${BASE}/tmp/comb_${TAG}_${CHR}
mkdir -p ${OUTDIR} ${TMP}

VARGS=""
for g in ${GVCFDIR}/*.g.vcf.gz; do VARGS="${VARGS} -V ${g}"; done
N=$(echo ${VARGS} | grep -o '\-V' | wc -l)
[ ${N} -eq 20 ] || { echo "ERROR: expected 20, got ${N}"; exit 1; }

echo "[$(date)] CombineGVCFs start ${TAG} ${CHR}"
gatk --java-options "-Xmx8g -Djava.io.tmpdir=${TMP}" CombineGVCFs \
  -R ${REF} ${VARGS} -L ${CHR} \
  -O ${TMP}/combined_${CHR}.g.vcf.gz

echo "[$(date)] GenotypeGVCFs start"
gatk --java-options "-Xmx8g -Djava.io.tmpdir=${TMP}" GenotypeGVCFs \
  -R ${REF} -V ${TMP}/combined_${CHR}.g.vcf.gz -L ${CHR} \
  -O ${OUTDIR}/${CHR}.vcf.gz

gatk IndexFeatureFile -I ${OUTDIR}/${CHR}.vcf.gz
echo "[$(date)] DONE ${TAG} ${CHR}"
rm -rf ${TMP}
