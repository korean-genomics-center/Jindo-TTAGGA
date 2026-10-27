#!/bin/bash
set -uo pipefail
conda activate gatk4
export JAVA_HOME=${CONDA_PREFIX}/lib/jvm
export PATH=$JAVA_HOME/bin:$PATH
JINDO=${JINDO_ROOT}
BASE=$JINDO/Analysis/Population_20Jindo/Reanalysis_Combined
WG=$BASE/WGS_Recall_Combined
REF=$BASE/REF/jindo_male_combined.fasta
JOINT_DIR=$WG/VCF_joint
CHRY_JOINT=$BASE/GATK_Recall_Combined/VCF_joint/chrY.combined.vcf.gz   # 파일럿 재활용
OUT=$WG/VCF_merged
STATS=$WG/Stats
mkdir -p $OUT $STATS

echo "===== wg04 aggregate : $(date) on $(hostname) ====="
# [1] 완료 점검 (autosome 38 + chrX = 39)
CHRS=( $(cat $WG/ref/autosomes.list | tr ' ' '\n') chrX )
MISS=0
for c in "${CHRS[@]}"; do [ -s "$JOINT_DIR/${c}.vcf.gz" ] || { echo "MISSING joint $c"; MISS=$((MISS+1)); }; done
[ -s "$CHRY_JOINT" ] || { echo "MISSING chrY joint (파일럿)"; MISS=$((MISS+1)); }
[ "$MISS" -ne 0 ] && { echo "ERROR $MISS missing. 중단."; exit 1; }
echo ">>> OK 39 + chrY present"

# [2] concat: .fai 순서대로 (auto + chrX + chrY)
LIST=$OUT/concat_list.txt; : > $LIST
for c in $(cat $WG/ref/autosomes.list | tr ' ' '\n'); do echo "$JOINT_DIR/${c}.vcf.gz" >> $LIST; done
echo "$JOINT_DIR/chrX.vcf.gz" >> $LIST
echo "$CHRY_JOINT" >> $LIST

ALL=$OUT/T2T_combined_all.vcf.gz
bcftools concat -f $LIST --allow-overlaps -Oz -o $ALL --threads 8 2>/dev/null || bcftools concat -f $LIST -Oz -o $ALL --threads 8
bcftools index -t $ALL

# [3] SNP 추출 + 필터 (기존과 동일 기준)
SNP=$OUT/T2T_combined_snp.vcf.gz
FLAG=$OUT/T2T_combined_snp_flagged.vcf.gz
gatk SelectVariants -R "$REF" -V $ALL --select-type-to-include SNP -O $SNP
gatk VariantFiltration -R "$REF" -V $SNP \
  --filter-expression "QD < 2.0" --filter-name "QD2" \
  --filter-expression "FS > 60.0" --filter-name "FS60" \
  --filter-expression "MQ < 40.0" --filter-name "MQ40" \
  --filter-expression "QUAL < 30.0" --filter-name "QUAL30" -O $FLAG

# [4] 집계 + old vs new
NEW_PASS=$(zcat $FLAG | grep -v '^#' | awk '$7=="PASS"' | wc -l)
OLD_PASS=13523155  # 기존 T2T 전체 PASS (SNP_SUMMARY 확정값)
echo "" | tee $STATS/wg_comparison.txt
echo "########## 전체 T2T PASS SNP: OLD vs NEW (combined/unmasked) ##########" | tee -a $STATS/wg_comparison.txt
echo "  OLD (masked, chrX 없음) : $OLD_PASS" | tee -a $STATS/wg_comparison.txt
echo "  NEW (combined, unmasked): $NEW_PASS" | tee -a $STATS/wg_comparison.txt
echo "" | tee -a $STATS/wg_comparison.txt

# [5] 염색체별 old vs new
echo "########## 염색체별 NEW PASS SNP ##########" | tee -a $STATS/wg_comparison.txt
zcat $FLAG | awk '$7=="PASS"{c[$1]++} END{for(k in c) print k"\t"c[k]}' | sort -V | tee -a $STATS/wg_comparison.txt
echo ""
echo "===== wg04 DONE : NEW=$NEW_PASS (old=$OLD_PASS) : $(date) =====" | tee -a $STATS/wg_comparison.txt
