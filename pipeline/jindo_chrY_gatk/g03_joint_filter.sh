#!/bin/bash
set -uo pipefail
conda activate gatk4
JINDO=${JINDO_ROOT}
BASE=$JINDO/Analysis/Population_20Jindo/Reanalysis_Combined
GC=$BASE/GATK_Recall_Combined
OLDGATK=$JINDO/Analysis/GATK_Recall
REF=$BASE/REF/jindo_male_combined.fasta
GVCF_DIR=$GC/gvcf_chrY
DB=$GC/genomicsdb_chrY
JOINT_DIR=$GC/VCF_joint
FILT_DIR=$GC/VCF_filtered
STATS=$GC/Stats
TMP=$GC/tmp/joint
mkdir -p $JOINT_DIR $FILT_DIR $STATS $TMP
echo "===== g03 joint+filter : $(date) on $(hostname) ====="
SAMPLES=( J288326 J352388 J433555 J495439 J496128 J498404 J543672 J61435 J61783 J63542 J12952315 J13044625 J288508 J62428 J62957 J63253 J73201 J73478 J76290 J90474 )
MISS=0
for s in "${SAMPLES[@]}"; do [ -s "$GVCF_DIR/${s}.g.vcf.gz" ] || { echo "MISSING gvcf $s"; MISS=$((MISS+1)); }; done
[ "$MISS" -ne 0 ] && { echo "ERROR $MISS gvcf missing. 중단."; exit 1; }
echo ">>> OK 20/20 gvcf"
VARGS=""; for g in $GVCF_DIR/*.g.vcf.gz; do VARGS="$VARGS -V $g"; done
N=$(echo $VARGS | grep -o '\-V' | wc -l); [ "$N" -eq 20 ] || { echo "ERROR got $N"; exit 1; }
rm -rf $DB
gatk --java-options "-Xmx8g -Djava.io.tmpdir=$TMP" GenomicsDBImport $VARGS --genomicsdb-workspace-path $DB -L chrY --tmp-dir $TMP --reader-threads 4
JOINT=$JOINT_DIR/chrY.combined.vcf.gz
gatk --java-options "-Xmx8g -Djava.io.tmpdir=$TMP" GenotypeGVCFs -R "$REF" -V gendb://$DB -O "$JOINT" --tmp-dir $TMP
gatk IndexFeatureFile -I "$JOINT"
SNP=$FILT_DIR/chrY.combined_snp.vcf.gz; FLAG=$FILT_DIR/chrY.combined_snp_flagged.vcf.gz
gatk SelectVariants -V "$JOINT" --select-type-to-include SNP -O "$SNP"
gatk VariantFiltration -V "$SNP" \
  --filter-expression "QD < 2.0" --filter-name "QD2" \
  --filter-expression "FS > 60.0" --filter-name "FS60" \
  --filter-expression "MQ < 40.0" --filter-name "MQ40" \
  --filter-expression "QUAL < 30.0" --filter-name "QUAL30" -O "$FLAG"
NEW_RAW=$(zcat "$SNP" | grep -vc '^#')
NEW_PASS=$(zcat "$FLAG" | grep -v '^#' | awk '$7=="PASS"' | wc -l)
OLD_PASS=$(zcat $OLDGATK/VCF_merged/T2T_snp_flagged.vcf.gz 2>/dev/null | awk '$1=="chrY" && $7=="PASS"' | wc -l)
ROS_PASS=$(zcat $OLDGATK/VCF_merged/ROS_snp_flagged.vcf.gz 2>/dev/null | awk '($1=="chrY"||$1=="Y") && $7=="PASS"' | wc -l)
echo "" | tee $STATS/comparison.txt
echo "########## SAME-CALLER 비교 (GATK, chrY PASS SNP) ##########" | tee -a $STATS/comparison.txt
echo "  OLD (chrX 없는 ref)   T2T : $OLD_PASS" | tee -a $STATS/comparison.txt
echo "  NEW (combined ref)    T2T : $NEW_PASS (raw SNP $NEW_RAW)" | tee -a $STATS/comparison.txt
echo "  ROS (기존 참고)           : $ROS_PASS" | tee -a $STATS/comparison.txt
echo "" | tee -a $STATS/comparison.txt
echo "########## chrY 1Mb bin (OLD vs NEW PASS) ##########" | tee -a $STATS/comparison.txt
echo -e "Mb\tOLD\tNEW" | tee -a $STATS/comparison.txt
paste \
  <(zcat $OLDGATK/VCF_merged/T2T_snp_flagged.vcf.gz | awk '$1=="chrY"&&$7=="PASS"{print int($2/1e6)}' | sort -n | uniq -c | awk '{print $2"\t"$1}') \
  <(zcat "$FLAG" | awk '$1=="chrY"&&$7=="PASS"{print int($2/1e6)}' | sort -n | uniq -c | awk '{print $2"\t"$1}') \
  2>/dev/null | awk '{printf "%s\t%s\t%s\n",$1,$2,$4}' | tee -a $STATS/comparison.txt
echo ""
echo "===== g03 ALL DONE : NEW=$NEW_PASS old=$OLD_PASS : $(date) =====" | tee -a $STATS/comparison.txt
