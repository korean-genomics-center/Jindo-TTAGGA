#!/bin/bash
set -eu
R=${JINDO_ROOT}/Results/Genome_Assembly/Assembly_Evaluation/GCI_Evaluation/GCI_Result
J=(Hap1:hifi Hap1:nano Hap1:two_type Hap2:hifi Hap2:nano Hap2:two_type)
IFS=':' read -r H T <<< "${J[$((SGE_TASK_ID-1))]}"
F=$(ls $R/$H/*_${T}.depth.gz)
O=depth_out/${H}_${T}_win1Mb_v3.tsv
zcat $F | awk -v W=1000000 '/^>/{c=substr($1,2); p=0; next} {w=int(p/W); k=c"\t"w; s[k]+=$1; n[k]++; p++} END{for(k in s) printf "%s\t%.3f\t%d\n",k,s[k]/n[k],n[k]}' | sort -k1,1V -k2,2n > $O.tmp
mv $O.tmp $O
echo "DONE $H $T $(wc -l < $O) windows"
