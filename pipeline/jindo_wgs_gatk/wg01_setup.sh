#!/bin/bash
set -e
conda activate gatk4
JINDO=${JINDO_ROOT}
BASE=$JINDO/Analysis/Population_20Jindo/Reanalysis_Combined
WG=$BASE/WGS_Recall_Combined
REF=$BASE/REF/jindo_male_combined.fasta
GATK=$JINDO/Analysis/GATK_Recall

echo "===== wg01 setup : $(date) on $(hostname) ====="
# dict 확인 (chrY 파일럿에서 이미 생성됨)
[ -s "${REF%.fasta}.dict" ] || gatk CreateSequenceDictionary -R "$REF"

# interval: 상염색체 38개(기존 순서 계승), chrX
cp $GATK/intervals/autosomes_hap1.list $WG/ref/autosomes.list
echo "chrX" > $WG/ref/chrX.list
echo "[1] autosome interval: $(cat $WG/ref/autosomes.list | tr ' ' '\n' | wc -l) chr"
echo "    chrX interval 준비"

# RepeatMasker .out -> BED (repeat/unique 분리용, 최종 단계에서 사용)
RMOUT=$JINDO/Results/Repeat_Annotation/RepeatMasker/Hap1/FinalASM.J495799_Child.Trio_Specific.hap1.fasta.out
RMBED=$WG/ref/hap1_repeat.bed
echo "[2] RepeatMasker .out -> BED"
awk 'NR>3 && $5!="" {print $5"\t"($6-1)"\t"$7}' "$RMOUT" | sort -k1,1 -k2,2n > "$RMBED"
echo "    repeat 구간: $(wc -l < $RMBED)"
# 병합 (bedtools 있으면)
if command -v bedtools >/dev/null 2>&1; then
  bedtools merge -i "$RMBED" > "${RMBED%.bed}.merged.bed"
  echo "    merged: $(wc -l < ${RMBED%.bed}.merged.bed) 구간, $(awk '{s+=$3-$2}END{print s}' ${RMBED%.bed}.merged.bed) bp"
else
  echo "    (bedtools 없음 — merge는 wg05에서)"
fi
echo "===== wg01 DONE : $(date) ====="
