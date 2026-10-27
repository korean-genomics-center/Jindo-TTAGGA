#!/bin/bash
set -uo pipefail
conda activate gatk4
JINDO=${JINDO_ROOT}
BASE=$JINDO/Analysis/Population_20Jindo/Reanalysis_Combined
GC=$BASE/GATK_Recall_Combined
REF=$BASE/REF/jindo_male_combined.fasta
BAM_DIR=$BASE/BAM
GVCF_DIR=$GC/gvcf_chrY
INTERVALS=$GC/ref/chrY.list
TMPBASE=$GC/tmp
mkdir -p $GVCF_DIR $TMPBASE
SAMPLES=( J288326 J352388 J433555 J495439 J496128 J498404 J543672 J61435 J61783 J63542 J12952315 J13044625 J288508 J62428 J62957 J63253 J73201 J73478 J76290 J90474 )
SAMPLE=${SAMPLES[$((SGE_TASK_ID-1))]}
BAM=$BAM_DIR/${SAMPLE}_combined.sorted.bam
OUT=$GVCF_DIR/${SAMPLE}.g.vcf.gz
DONE=$GVCF_DIR/${SAMPLE}.DONE
TMP=$TMPBASE/${SAMPLE}; mkdir -p $TMP
echo "===== g02 HC $SGE_TASK_ID/20 : $SAMPLE : $(date) on $(hostname) ====="
if [ -f "$DONE" ] && [ -s "$OUT" ]; then echo "SKIP done"; exit 0; fi
[ -s "$BAM" ] || { echo "ERROR BAM missing"; exit 1; }
gatk --java-options "-Xmx7g -Djava.io.tmpdir=${TMP}" HaplotypeCaller \
  -R "$REF" -I "$BAM" -O "$OUT" -L "$INTERVALS" --sample-ploidy 1 \
  -ERC GVCF --native-pair-hmm-threads 16
[ -s "$OUT" ] || { echo "ERROR HC failed"; exit 1; }
gatk IndexFeatureFile -I "$OUT"; rm -rf $TMP; touch "$DONE"
echo "===== g02 DONE $SAMPLE : $(date) ====="
