#!/usr/bin/env bash
# Locate the pseudoautosomal region (PAR) of the Jindo1-G-TTAGGA chrY by aligning
# TTAGGA-P chrY to TTAGGA-M chrX, and check that ROS_Cfam_1.0 carries the PAR on chrX.
# Usage: JINDO_ROOT=/path/to/data THREADS=8 bash chrY_par_alignment_v1.sh
# Requires: samtools, minimap2 (v2.28-r1209 used in the paper)
set -euo pipefail
: "${JINDO_ROOT:?set JINDO_ROOT}"; THREADS=${THREADS:-8}
M_FA=${M_FA:-$JINDO_ROOT/assembly/Jindo1-G-TTAGGA.hap1.fasta}   # TTAGGA-M (carries chrX)
P_FA=${P_FA:-$JINDO_ROOT/assembly/Jindo1-G-TTAGGA.hap2.fasta}   # TTAGGA-P (carries chrY)
ROS_FA=${ROS_FA:-$JINDO_ROOT/reference/GCF_014441545.1_ROS_Cfam_1.0_genomic.fna}
OUT=${OUT:-par_alignment}; mkdir -p "$OUT"

samtools faidx "$M_FA" chrX > "$OUT/M_chrX.fa"
samtools faidx "$P_FA" chrY > "$OUT/P_chrY.fa"

# 1. TTAGGA-P chrY vs TTAGGA-M chrX
minimap2 -x asm20 -c -t "$THREADS" "$OUT/M_chrX.fa" "$OUT/P_chrY.fa" > "$OUT/P_chrY_vs_M_chrX.paf"
echo "== alignment blocks >=10 kb (chrY Mb, chrX Mb, strand, aligned kb, identity)"
awk '$11>=10000{for(i=13;i<=NF;i++) if($i~/^de:f:/) d=substr($i,6);
  printf "Y %.2f-%.2f\tX %.2f-%.2f\t%s\t%d kb\tid %.3f\n",$3/1e6,$4/1e6,$8/1e6,$9/1e6,$5,$11/1000,1-d}' \
  "$OUT/P_chrY_vs_M_chrX.paf" | sort -k2,2n
echo "== chrY length aligned to chrX per 1-Mb bin (blocks >=10 kb, identity >=0.95)"
awk '$11>=10000{for(i=13;i<=NF;i++) if($i~/^de:f:/) d=substr($i,6);
  if(1-d>=0.95) for(p=$3;p<$4;p+=1000) c[int(p/1e6)]+=1}
  END{for(i=0;i<=21;i++) printf "%d-%d Mb\t%d kb\n",i,i+1,c[i]}' "$OUT/P_chrY_vs_M_chrX.paf"

# 2. Jindo PAR (chrY 14.70 Mb to end) vs ROS_Cfam_1.0 chrX (first 10 Mb)
YLEN=$(awk '$1=="chrY"{print $2}' "$P_FA.fai")
samtools faidx "$OUT/P_chrY.fa" "chrY:14700001-$YLEN" > "$OUT/P_chrY_PAR.fa"
[ -f "$ROS_FA.fai" ] || samtools faidx "$ROS_FA"
samtools faidx "$ROS_FA" NC_051843.1:1-10000000 > "$OUT/ROS_chrX_10Mb.fa"
minimap2 -x asm20 -c -t "$THREADS" "$OUT/ROS_chrX_10Mb.fa" "$OUT/P_chrY_PAR.fa" > "$OUT/P_PAR_vs_ROS_chrX.paf"
echo "== Jindo PAR vs ROS_Cfam_1.0 chrX"
awk '$11>=10000{s+=$11; for(i=13;i<=NF;i++) if($i~/^de:f:/) d=substr($i,6);
  printf "PAR %.2f-%.2f Mb -> ROS chrX %.2f-%.2f Mb\t%s\t%d kb\tid %.3f\n",$3/1e6,$4/1e6,$8/1e6,$9/1e6,$5,$11/1000,1-d}
  END{printf "total aligned %.2f Mb\n",s/1e6}' "$OUT/P_PAR_vs_ROS_chrX.paf"
