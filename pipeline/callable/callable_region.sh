#!/bin/bash
set -eu
set +u; conda activate TTT; set -u

JINDO=${JINDO_ROOT}
POP=$JINDO/Analysis/Population_20Jindo
OUT=$POP/Reanalysis_Combined/Callable
mkdir -p $OUT/perbam $OUT/Logs

SAMPLES=($(ls $POP/Reanalysis_Combined/BAM/*_combined.sorted.bam | xargs -n1 basename | sed 's/_combined.sorted.bam//' | sort))
N=${#SAMPLES[@]}
TID=$((SGE_TASK_ID-1))
SIDX=$((TID % N)); REFIDX=$((TID / N))
SAMPLE=${SAMPLES[$SIDX]}
MINDP=10

if [ $REFIDX -eq 0 ]; then
  BAM=$POP/Reanalysis_Combined/BAM/${SAMPLE}_combined.sorted.bam; TAG=t2t
else
  BAM=$POP/BAM_ROS/jindo_${SAMPLE}_ros.sorted.bam; TAG=ros
fi
OUTBED=$OUT/perbam/${TAG}_${SAMPLE}.callable.bed.gz
[ -s "$OUTBED" ] && { echo "skip $TAG $SAMPLE"; exit 0; }

# pipefail 없이도 안전: 중간 파일 없이 파이프. 실패 시 빈 출력은 아래 크기체크로 감지.
samtools view -b -q 20 "$BAM" \
  | bedtools genomecov -ibam - -bga \
  | awk -v m=$MINDP '$4>=m {print $1"\t"$2"\t"$3}' \
  | bedtools merge -i - \
  | gzip > "$OUTBED"

# 결과 검증: 비었으면 실패로 처리
if [ "$(zcat "$OUTBED" | head -1 | wc -l)" -eq 0 ]; then
  echo "FAIL $TAG $SAMPLE : empty output"; rm -f "$OUTBED"; exit 1
fi
echo "DONE $TAG $SAMPLE : $(zcat "$OUTBED" | wc -l) intervals"
