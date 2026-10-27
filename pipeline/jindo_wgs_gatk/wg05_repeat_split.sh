#!/bin/bash
set -uo pipefail
set +u
conda activate TTT
set -u
command -v bedtools >/dev/null || { echo "FATAL: bedtools 없음"; exit 1; }
JINDO=${JINDO_ROOT}
BASE=$JINDO/Analysis/Population_20Jindo/Reanalysis_Combined
WG=$BASE/WGS_Recall_Combined
OUT=$WG/VCF_merged
STATS=$WG/Stats
RMBED=$WG/ref/combined_repeat.bed
FLAG=$OUT/T2T_combined_snp_flagged.vcf.gz

echo "===== wg05 repeat/unique 분리 : $(date) on $(hostname) ====="
# repeat BED merge (wg01에서 안했으면)
MERGED=${RMBED%.bed}.merged.bed
[ -s "$MERGED" ] || bedtools merge -i "$RMBED" > "$MERGED"

# PASS SNP만 BED로
PASSBED=$OUT/T2T_pass_snp.bed
zcat $FLAG | awk '$7=="PASS" && $1!~/^#/{print $1"\t"($2-1)"\t"$2}' > $PASSBED
TOTAL=$(wc -l < $PASSBED)

# repeat 영역과 교차 (repeat 내 / unique)
IN_REPEAT=$(bedtools intersect -a $PASSBED -b $MERGED -u | wc -l)
UNIQUE=$((TOTAL - IN_REPEAT))

echo "" | tee $STATS/wg_repeat_split.txt
echo "########## T2T PASS SNP: repeat vs unique ##########" | tee -a $STATS/wg_repeat_split.txt
echo "  전체 PASS SNP : $TOTAL" | tee -a $STATS/wg_repeat_split.txt
echo "  repeat 영역   : $IN_REPEAT ($(awk "BEGIN{printf \"%.1f\", 100*$IN_REPEAT/$TOTAL}")%)" | tee -a $STATS/wg_repeat_split.txt
echo "  unique 영역   : $UNIQUE ($(awk "BEGIN{printf \"%.1f\", 100*$UNIQUE/$TOTAL}")%)  <- 신뢰 변이" | tee -a $STATS/wg_repeat_split.txt
echo "" | tee -a $STATS/wg_repeat_split.txt
echo "[다음 수동 과제] ROS도 동일하게 unique 영역 산출 후 T2T vs ROS unique 차집합 = 방어가능한 +N" | tee -a $STATS/wg_repeat_split.txt
echo "===== wg05 DONE : unique=$UNIQUE / total=$TOTAL : $(date) =====" | tee -a $STATS/wg_repeat_split.txt
