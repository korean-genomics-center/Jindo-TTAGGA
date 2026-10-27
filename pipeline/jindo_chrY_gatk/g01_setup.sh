#!/bin/bash
set -e
conda activate gatk4
JINDO=${JINDO_ROOT}
BASE=$JINDO/Analysis/Population_20Jindo/Reanalysis_Combined
GC=$BASE/GATK_Recall_Combined
REF=$BASE/REF/jindo_male_combined.fasta
echo "===== g01 setup : $(date) on $(hostname) ====="
if [ ! -s "${REF%.fasta}.dict" ]; then
  echo "[1] CreateSequenceDictionary"; gatk CreateSequenceDictionary -R "$REF"
else echo "[1] dict exists, skip"; fi
echo "chrY" > $GC/ref/chrY.list
echo "[2] interval:"; cat $GC/ref/chrY.list
echo "[3] ref 점검:"; ls -la ${REF} ${REF}.fai ${REF%.fasta}.dict
echo "===== g01 DONE : $(date) ====="
