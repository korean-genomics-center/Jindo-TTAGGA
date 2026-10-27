#!/bin/bash
set -uo pipefail
conda activate gatk4
JINDO=${JINDO_ROOT}
BASE=$JINDO/Analysis/Population_20Jindo/Reanalysis_Combined
WG=$BASE/WGS_Recall_Combined
REF=$BASE/REF/jindo_male_combined.fasta
BAM_DIR=$BASE/BAM
TMPBASE=$WG/tmp
mkdir -p $TMPBASE

SAMPLES=( J288326 J352388 J433555 J495439 J496128 J498404 J543672 J61435 J61783 J63542 J12952315 J13044625 J288508 J62428 J62957 J63253 J73201 J73478 J76290 J90474 )
# task 1-20 = autosome, task 21-40 = chrX
IDX=$((SGE_TASK_ID-1))
if [ $IDX -lt 20 ]; then
  REGION=auto; PLOIDY=2; INTERVALS=$WG/ref/autosomes.list
  SAMPLE=${SAMPLES[$IDX]}
else
  REGION=chrX; PLOIDY=1; INTERVALS=$WG/ref/chrX.list
  SAMPLE=${SAMPLES[$((IDX-20))]}
fi
GVCF_DIR=$WG/gvcf_${REGION}
mkdir -p $GVCF_DIR
BAM=$BAM_DIR/${SAMPLE}_combined.sorted.bam
OUT=$GVCF_DIR/${SAMPLE}.g.vcf.gz
DONE=$GVCF_DIR/${SAMPLE}.DONE
TMP=$TMPBASE/${SAMPLE}_${REGION}; mkdir -p $TMP

echo "===== wg02 task $SGE_TASK_ID/40 : $SAMPLE [$REGION ploidy=$PLOIDY] : $(date) on $(hostname) ====="
if [ -f "$DONE" ] && [ -s "$OUT" ]; then echo "SKIP done"; exit 0; fi
[ -s "$BAM" ] || { echo "ERROR BAM missing $BAM"; exit 1; }

gatk --java-options "-Xmx16g -Djava.io.tmpdir=${TMP}" HaplotypeCaller \
  -R "$REF" -I "$BAM" -O "$OUT" \
  -L "$INTERVALS" --sample-ploidy $PLOIDY \
  -ERC GVCF --native-pair-hmm-threads 16
[ -s "$OUT" ] || { echo "ERROR HC failed $SAMPLE $REGION"; exit 1; }
gatk IndexFeatureFile -I "$OUT"
rm -rf $TMP; touch "$DONE"
echo "===== wg02 DONE $SAMPLE [$REGION] : $(date) ====="
