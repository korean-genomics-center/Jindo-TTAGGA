#!/bin/bash
set -uo pipefail
export PATH=${CONDA_PREFIX}/bin:$PATH
JINDO=${JINDO_ROOT}
BASE=$JINDO/Analysis/Population_20Jindo/Reanalysis_Combined
REF=$BASE/REF/jindo_male_combined.fasta
BAM_DIR=$BASE/BAM; STATS=$BASE/Stats; LOG_DIR=$BASE/Logs; TMP_DIR=$BASE/tmp
mkdir -p "$BAM_DIR" "$STATS" "$LOG_DIR" "$TMP_DIR"
THREADS=80
F001=$JINDO/Resources/SequencingData/WGS_Other_20/Jindo_WGS_20/HDD1/X209SC25127458-Z01-F001/01.RawData
F002=$JINDO/Resources/SequencingData/WGS_Other_20/Jindo_WGS_20/HDD1/X209SC25127458-Z01-F002/01.RawData
SAMPLES=( J288326 J352388 J433555 J495439 J496128 J498404 J543672 J61435 J61783 J63542 J12952315 J13044625 J288508 J62428 J62957 J63253 J73201 J73478 J76290 J90474 )
SAMPLE=${SAMPLES[$((SGE_TASK_ID-1))]}
echo "===== Realign $SGE_TASK_ID/20 : $SAMPLE : $(date) on $(hostname) t=$THREADS ====="
OUT_BAM=$BAM_DIR/${SAMPLE}_combined.sorted.bam; DONE_FLAG=$BAM_DIR/${SAMPLE}.DONE
if [ -f "$DONE_FLAG" ] && [ -s "$OUT_BAM" ]; then echo "SKIP $SAMPLE done"; exit 0; fi
[ -s "$REF" ] && [ -s "${REF}.bwt" ] || { echo "ERROR ref/index missing"; exit 1; }
R1_LIST=(); R2_LIST=()
for DIR in "$F001/$SAMPLE" "$F002/$SAMPLE"; do
  [ -d "$DIR" ] || continue
  for r1 in "$DIR"/*_1.fq.gz; do [ -e "$r1" ] || continue
    r2=${r1/_1.fq.gz/_2.fq.gz}; [ -s "$r2" ] || { echo "ERROR R2 missing $r1"; exit 1; }
    R1_LIST+=("$r1"); R2_LIST+=("$r2"); done
done
NLANE=${#R1_LIST[@]}; [ "$NLANE" -eq 0 ] && { echo "ERROR no fastq $SAMPLE"; exit 1; }
echo "  lanes: $NLANE"
LANE_BAMS=(); i=0
for idx in "${!R1_LIST[@]}"; do i=$((i+1)); R1=${R1_LIST[$idx]}; R2=${R2_LIST[$idx]}
  LANE_TAG=$(basename "$R1" | sed -E 's/.*_([0-9A-Za-z]+)_L([0-9]+)_1.fq.gz/\1_L\2/')
  RG="@RG\tID:${SAMPLE}.${LANE_TAG}\tSM:${SAMPLE}\tLB:${SAMPLE}\tPL:ILLUMINA\tPU:${LANE_TAG}"
  LANE_BAM=$TMP_DIR/${SAMPLE}.lane${i}.bam
  echo "  [lane $i/$NLANE] $LANE_TAG"
  # bwa -> name-sort -> fixmate(-m) -> coord-sort  (markdup 전제 조건)
  bwa mem -t $THREADS -R "$RG" "$REF" "$R1" "$R2" 2>> "$LOG_DIR/${SAMPLE}.bwa.err" \
    | samtools sort -n -@ 8 -m 2G -T "$TMP_DIR/${SAMPLE}.lane${i}.nst" -O bam - \
    | samtools fixmate -m -@ 8 - - \
    | samtools sort -@ 8 -m 2G -T "$TMP_DIR/${SAMPLE}.lane${i}.cst" -o "$LANE_BAM" -
  [ -s "$LANE_BAM" ] || { echo "ERROR bwa/fixmate/sort $SAMPLE lane $i"; exit 1; }
  LANE_BAMS+=("$LANE_BAM")
done
MERGED=$TMP_DIR/${SAMPLE}.merged.bam
if [ "${#LANE_BAMS[@]}" -eq 1 ]; then mv "${LANE_BAMS[0]}" "$MERGED"
else echo "  merge ${#LANE_BAMS[@]} lanes"; samtools merge -@ 16 -f "$MERGED" "${LANE_BAMS[@]}" || { echo "ERROR merge"; exit 1; }; rm -f "${LANE_BAMS[@]}"; fi
echo "  markdup"; samtools markdup -@ 16 "$MERGED" "$OUT_BAM" || { echo "ERROR markdup"; exit 1; }; rm -f "$MERGED"
samtools index -@ 8 "$OUT_BAM"; samtools flagstat "$OUT_BAM" > "$STATS/${SAMPLE}_combined_flagstat.txt"
touch "$DONE_FLAG"; echo "===== DONE $SAMPLE : $(date) ====="
