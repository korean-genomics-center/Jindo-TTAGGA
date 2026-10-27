#!/bin/bash
set -uo pipefail
export PATH=${CONDA_PREFIX}/bin:$PATH
JINDO=${JINDO_ROOT}
POP=$JINDO/Analysis/Population_20Jindo
BASE=$POP/Reanalysis_Combined
REF=$BASE/REF/jindo_male_combined.fasta
BAM_DIR=$BASE/BAM; OLD_BAM_DIR=$POP/BAM_Hap2; VCF_DIR=$BASE/VCF_chrY; STATS=$BASE/Stats
mkdir -p "$VCF_DIR" "$STATS"
SAMPLES=( J288326 J352388 J433555 J495439 J496128 J498404 J543672 J61435 J61783 J63542 J12952315 J13044625 J288508 J62428 J62957 J63253 J73201 J73478 J76290 J90474 )
echo "===== chrY stats + variant recall : $(date) on $(hostname) ====="
echo "########## [1] 완료 점검 ##########"
MISSING=0
for s in "${SAMPLES[@]}"; do b=$BAM_DIR/${s}_combined.sorted.bam
  { [ -s "$b" ] && [ -s "${b}.bai" ]; } || { echo "  MISSING $s"; MISSING=$((MISSING+1)); }; done
[ "$MISSING" -ne 0 ] && { echo ">>> ERROR $MISSING missing. 중단."; exit 1; }
echo ">>> OK 20/20 present"
echo "########## [2] chrY MAPQ0 before vs after ##########"
OUT=$STATS/chrY_MAPQ0_before_after.tsv
printf "sample\tbefore_MAPQ0%%\tbefore_total\tafter_MAPQ0%%\tafter_total\n" > "$OUT"
for s in "${SAMPLES[@]}"; do
  OLDB=$OLD_BAM_DIR/jindo_${s}_hap2.sorted.bam; NEWB=$BAM_DIR/${s}_combined.sorted.bam
  if [ -s "$OLDB" ]; then bstat=$(samtools view "$OLDB" chrY 2>/dev/null | awk '{if($5==0)z++;t++} END{if(t>0)printf "%.1f\t%d",100*z/t,t; else printf "NA\t0"}'); else bstat="NA\t0"; fi
  nstat=$(samtools view "$NEWB" chrY 2>/dev/null | awk '{if($5==0)z++;t++} END{if(t>0)printf "%.1f\t%d",100*z/t,t; else printf "NA\t0"}')
  printf "%s\t%b\t%b\n" "$s" "$bstat" "$nstat" >> "$OUT"; done
echo "  -> $OUT"; column -t "$OUT"
echo "########## [3] chrX/chrY depth + X:Y ratio ##########"
OUT2=$STATS/chrX_chrY_depth_ratio.tsv
printf "sample\tchrX_dp\tchrY_dp\tX/Y_ratio\n" > "$OUT2"
for s in "${SAMPLES[@]}"; do b=$BAM_DIR/${s}_combined.sorted.bam
  x=$(samtools depth -r chrX "$b" 2>/dev/null | awk '{sum+=$3;n++} END{if(n>0)printf "%.1f",sum/n; else print 0}')
  y=$(samtools depth -r chrY "$b" 2>/dev/null | awk '{sum+=$3;n++} END{if(n>0)printf "%.1f",sum/n; else print 0}')
  r=$(echo "$x $y" | awk '{if($2>0)printf "%.2f",$1/$2; else print "NA"}')
  printf "%s\t%s\t%s\t%s\n" "$s" "$x" "$y" "$r" >> "$OUT2"; done
echo "  -> $OUT2"; column -t "$OUT2"
echo "########## [4] haploid chrY 변이 재호출 ##########"
BAM_LIST=$VCF_DIR/bam_list_chrY.txt; : > "$BAM_LIST"
for s in "${SAMPLES[@]}"; do echo "$BAM_DIR/${s}_combined.sorted.bam" >> "$BAM_LIST"; done
N=$(wc -l < "$BAM_LIST"); echo "  BAMs: $N/20"; [ "$N" -ne 20 ] && { echo "ERROR"; exit 1; }
RAW=$VCF_DIR/cohort_20jindo.chrY.combined.raw.vcf.gz
FILT=$VCF_DIR/cohort_20jindo.chrY.combined.filtered.vcf.gz
bcftools mpileup -f "$REF" -b "$BAM_LIST" -r chrY -a "FORMAT/AD,FORMAT/DP,FORMAT/SP,INFO/AD" -q 20 -Q 20 --threads 32 -Ou | bcftools call -m --ploidy 1 -v --threads 32 -Oz -o "$RAW"
bcftools index -t "$RAW"
bcftools filter -e 'QUAL<30 || INFO/DP<40' -Oz -o "$FILT" "$RAW"; bcftools index -t "$FILT"
NVAR_RAW=$(bcftools view -H "$RAW" 2>/dev/null | wc -l); NVAR_FILT=$(bcftools view -H "$FILT" 2>/dev/null | wc -l)
echo "  chrY variants raw=$NVAR_RAW filtered=$NVAR_FILT"
bcftools stats "$FILT" > "$STATS/chrY_combined_vcf_stats.txt"
echo "===== ALL DONE : MAPQ0=$OUT | ratio=$OUT2 | chrY_filt=$NVAR_FILT : $(date) ====="
