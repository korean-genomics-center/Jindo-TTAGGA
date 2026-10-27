#!/bin/bash
set -e
conda activate gatk4
# 인자: $1=sample $2=ref_tag(hap1|hap2|ros) $3=ploidy $4=interval_list $5=tag(출력구분)
BASE=${JINDO_ROOT}/Analysis/GATK_Recall
SAMPLE=$1; REFTAG=$2; PLOIDY=$3; INTERVALS=$4; TAG=$5
REF=${BASE}/ref/${REFTAG}.fasta
TMP=${BASE}/tmp/${SAMPLE}_${TAG}
mkdir -p ${TMP}

case ${REFTAG} in
  hap1) BAM=${BASE}/../Population_20Jindo/BAM/jindo_${SAMPLE}_hap1.sorted.bam ;;
  hap2) BAM=${BASE}/../Population_20Jindo/BAM_Hap2/jindo_${SAMPLE}_hap2.sorted.bam ;;
  ros)  BAM=${BASE}/../Population_20Jindo/BAM_ROS/jindo_${SAMPLE}_ros.sorted.bam ;;
esac

OUTDIR=${BASE}/gvcf_${TAG}
mkdir -p ${OUTDIR}
OUT=${OUTDIR}/${SAMPLE}.g.vcf.gz

AVAIL=$(df -P ${OUTDIR} | awk 'NR==2{print $4}')
echo "[$(date)] ${SAMPLE} tag=${TAG} ref=${REFTAG} ploidy=${PLOIDY} avail_KB=${AVAIL} host=$(hostname)"
if [ ${AVAIL} -lt 52428800 ]; then echo "ERROR: <50GB free, abort"; exit 1; fi
echo "BAM=${BAM}"; echo "REF=${REF}"; echo "INTERVALS=${INTERVALS}"
[ -f ${BAM} ] || { echo "ERROR: BAM not found"; exit 1; }

gatk --java-options "-Xmx12g -Djava.io.tmpdir=${TMP}" HaplotypeCaller \
  -R ${REF} -I ${BAM} -O ${OUT} \
  -L ${INTERVALS} --sample-ploidy ${PLOIDY} \
  -ERC GVCF --native-pair-hmm-threads 8

gatk IndexFeatureFile -I ${OUT}
echo "[$(date)] DONE ${SAMPLE} tag=${TAG}"
rm -rf ${TMP}
