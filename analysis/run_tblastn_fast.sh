#!/bin/bash
set -eu
source ${WORK_DIR}/miniconda3/etc/profile.d/conda.sh
set +u; conda activate TTT; set -u

OUT=${JINDO_ROOT}/Results/Gene_Annotation/RefAbsent_Functional/tblastn_genome
TAGS=(ros gsd boxer zoey)
TID=$((SGE_TASK_ID-1))
GI=$((TID / 20)); CI=$((TID % 20 + 1))
TAG=${TAGS[$GI]}
Q=$(printf "%s/chunks/q%02d.faa" $OUT $CI)
O=$(printf "%s/results_chunked/%s_q%02d.tsv" $OUT $TAG $CI)

[ -s "$O" ] && { echo "skip $TAG q$CI"; exit 0; }

tblastn -query $Q -db $OUT/db/$TAG \
        -evalue 1e-5 \
        -outfmt "6 qseqid sseqid pident length qlen slen qstart qend evalue bitscore" \
        -max_target_seqs 5 -num_threads 8 -out $O
echo "DONE $TAG chunk$CI : $(wc -l < $O) hits"
