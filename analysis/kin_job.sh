#!/bin/bash
set -uo pipefail
ROOT=${JINDO_ROOT}
CB=${WORK_DIR}/miniconda3
BCF=$CB/envs/TTT/bin/bcftools
PLK=$CB/envs/TTT/bin/plink
NP=32
cd ${WORK_DIR}
echo "=== $(date) | $(hostname) ==="

PC=$ROOT/Analysis/Population_20Jindo/VCF/per_chr_t2t
ls $PC/cohort_chr[0-9]*.vcf.gz 2>/dev/null | grep -viE "chr(X|Y|M)" | sort -V > vcflist.txt
echo "상염색체 VCF $(wc -l < vcflist.txt) 개"
echo ">> 샘플 확인 (첫 파일)"
$BCF query -l $(head -1 vcflist.txt) | tr '\n' ' '; echo

$BCF concat -f vcflist.txt --threads $NP -Oz -o autosomes.vcf.gz && $BCF index -f autosomes.vcf.gz
echo ">> 사이트 $($BCF view -H autosomes.vcf.gz | wc -l)"
echo ">> 샘플 $($BCF query -l autosomes.vcf.gz | wc -l) 마리"

$PLK --vcf autosomes.vcf.gz --double-id --dog --maf 0.05 --geno 0.1 \
     --make-bed --out auto --threads $NP 2>&1 | tail -4
if [ -f auto.bed ]; then
  $PLK --bfile auto --dog --genome --out kinship --threads $NP 2>&1 | tail -3
  if [ -f kinship.genome ]; then
    echo ">> PI_HAT 상위 25쌍"
    awk 'NR>1{printf "   %-11s %-11s Z0=%.3f Z1=%.3f Z2=%.3f PI_HAT=%.4f\n",$2,$4,$7,$8,$9,$10}' kinship.genome \
      | sort -t= -k5,5gr | head -25
    echo ">> 분포"
    awk 'NR>1{p=$10; if(p>0.354)a++; else if(p>0.177)b++; else if(p>0.0884)c++; else d++}
      END{printf "   1촌급(>0.354)  %d\n   2촌급(>0.177)  %d\n   3촌급(>0.0884) %d\n   무관(<=0.0884) %d\n   총 %d쌍\n",a,b,c,d,a+b+c+d}' kinship.genome
    awk 'NR>1{s+=$10;n++; if($10>m)m=$10}END{printf ">> 평균 PI_HAT %.4f / 최대 %.4f\n",s/n,m}' kinship.genome
    echo ">> lineage별 (S-Data 5c 배정 필요시 별도)"
  fi
fi
echo "=== 완료 $(date) ==="
