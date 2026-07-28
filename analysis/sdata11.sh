#!/bin/bash
set -uo pipefail
J=${JINDO_ROOT}
BCF=${WORK_DIR}/miniconda3/envs/TTT/bin/bcftools
H=$J/Analysis/Population_20Jindo/chrY_Analysis/Haplogroup
cd ${WORK_DIR}

for TAG in ros trunc4Mb; do
  case $TAG in
    ros)      V=$H/ros_chrY_pass.vcf.gz ;;
    trunc4Mb) V=$H/chrY_trunc4Mb.vcf.gz ;;
  esac
  [ -f "$V" ] || { echo "없음: $V"; continue; }
  echo "=== $TAG : $V ==="
  $BCF query -l "$V" > ${TAG}_samples.txt
  echo "  샘플 $(wc -l < ${TAG}_samples.txt) / 사이트 $($BCF view -H "$V" | wc -l)"
  $BCF query -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n' "$V" > ${TAG}_gt.tsv
  echo "  ${TAG}_gt.tsv 생성"
done
