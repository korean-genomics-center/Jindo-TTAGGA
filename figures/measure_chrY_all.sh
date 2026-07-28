#!/bin/bash
# chrY length / N / gap measurement from FASTA
# Single definition applied to every assembly: gap = run of >= 10 consecutive N
set -u
REFDIR=${JINDO_ROOT}/Resources/Reference
OUTDIR=${JINDO_ROOT}/Results/Manuscript_Figures/data/stats
OUT="$OUTDIR/chrY_measured.tsv"
mkdir -p "$OUTDIR"

echo -e "assembly\tseq_id\tclass\tlength_bp\tN_bp\tgaps_ge10bp" > "$OUT"

measure() {  # $1=fasta $2=seqid $3=assembly $4=class
  local seq
  seq=$(samtools faidx "$1" "$2" 2>/dev/null | grep -v ">" | tr -d '\n')
  [ -z "$seq" ] && return
  local L N G
  L=${#seq}
  N=$(printf '%s' "$seq" | tr -cd 'Nn' | wc -c)
  G=$(printf '%s' "$seq" | grep -o '[Nn]\+' | awk 'length($0)>=10' | wc -l)
  printf "%s\t%s\t%s\t%d\t%d\t%d\n" "$3" "$2" "$4" "$L" "$N" "$G" >> "$OUT"
}

for d in "$REFDIR"/GC*; do
  FA=$(ls "$d"/*_genomic.fna 2>/dev/null | head -1)
  [ -z "$FA" ] && continue
  ASM=$(basename "$FA" _genomic.fna | sed 's/^[^_]*_[^_]*_//')
  [ ! -f "${FA}.fai" ] && samtools faidx "$FA"
  grep -i "chromosome Y," "$FA" | grep -vi "unlocalized" | sed 's/^>//' | awk '{print $1}' \
    | while read -r s; do measure "$FA" "$s" "$ASM" "chr-scale"; done
  grep -i "chromosome Y unlocalized" "$FA" | sed 's/^>//' | awk '{print $1}' \
    | while read -r s; do measure "$FA" "$s" "$ASM" "unlocalized"; done
done
column -t "$OUT"
