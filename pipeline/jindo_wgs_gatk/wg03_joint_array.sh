#!/bin/bash
set -uo pipefail
conda activate gatk4
JINDO=${JINDO_ROOT}
BASE=$JINDO/Analysis/Population_20Jindo/Reanalysis_Combined
WG=$BASE/WGS_Recall_Combined
REF=$BASE/REF/jindo_male_combined.fasta

# chr 목록: autosome 38 + chrX (task index -> chr)
CHRS=( $(cat $WG/ref/autosomes.list | tr ' ' '\n') chrX )
CHR=${CHRS[$((SGE_TASK_ID-1))]}
# chrX는 gvcf_chrX, 나머지는 gvcf_auto
if [ "$CHR" = "chrX" ]; then GVCF_DIR=$WG/gvcf_chrX; else GVCF_DIR=$WG/gvcf_auto; fi

DB=$WG/gdb/${CHR}
JOINT_DIR=$WG/VCF_joint
TMP=$WG/tmp/joint_${CHR}
mkdir -p $JOINT_DIR $TMP $WG/gdb

echo "===== wg03 joint task $SGE_TASK_ID/39 : $CHR : $(date) on $(hostname) ====="
OUT=$JOINT_DIR/${CHR}.vcf.gz
DONE=$JOINT_DIR/${CHR}.DONE
if [ -f "$DONE" ] && [ -s "$OUT" ]; then echo "SKIP $CHR done"; exit 0; fi

VARGS=""
for g in $GVCF_DIR/*.g.vcf.gz; do VARGS="$VARGS -V $g"; done
N=$(echo $VARGS | grep -o '\-V' | wc -l)
[ "$N" -eq 20 ] || { echo "ERROR $CHR: expected 20 gvcf, got $N"; exit 1; }

rm -rf $DB
gatk --java-options "-Xmx16g -Djava.io.tmpdir=$TMP" GenomicsDBImport \
  $VARGS --genomicsdb-workspace-path $DB -L $CHR --tmp-dir $TMP --reader-threads 4 --batch-size 5
gatk --java-options "-Xmx16g -Djava.io.tmpdir=$TMP" GenotypeGVCFs \
  -R "$REF" -V gendb://$DB -O "$OUT" --tmp-dir $TMP
gatk IndexFeatureFile -I "$OUT"
rm -rf $TMP $DB; touch "$DONE"
echo "===== wg03 DONE $CHR : $(date) ====="
