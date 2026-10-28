#!/bin/bash
set -eu
source ${WORK_DIR}/miniconda3/etc/profile.d/conda.sh
set +u; conda activate TTT; set -u
J=${JINDO_ROOT}
python3 -c "
s=''.join(l.strip() for l in open('$J/Results/Manuscript_Figures/data/centromere/fasta/Hap2_chrY.fa') if not l.startswith('>'))
a,b=7269076,10826466; s=s[:a-1]+'N'*(b-a+1)+s[b:]
open('Hap2_chrY_cenmasked_v1.fa','w').write('>chrY\n'+'\n'.join(s[i:i+80] for i in range(0,len(s),80))+'\n'); print('masked chrY', len(s))"
Y=$PWD/Hap2_chrY_cenmasked_v1.fa
echo "start minimap2 $(date)"
$J/Resources/Tools/minimap2/minimap2_v2.30_r1287/minimap2/minimap2 -t 16 -c -X -x asm5 -N 20 -p 0.1 $Y $Y > chrY_self_v2.paf 2> chrY_self_v2.log
echo "end minimap2 $(date) | alignments $(wc -l < chrY_self_v2.paf)"
awk '$11>=10000 && $10/$11>=0.95{print $1"\t"$3"\t"$4; print $6"\t"$8"\t"$9}' chrY_self_v2.paf | sort -k2,2n | awk 'BEGIN{OFS="\t"} {if(NR>1 && $2<=e){if($3>e)e=$3} else {if(NR>1)print c,s,e; c=$1;s=$2;e=$3}} END{if(NR>0)print c,s,e}' > chrY_SD_v2.bed
awk '{printf "  SD %.2f-%.2f Mb\n",$2/1e6,$3/1e6}' chrY_SD_v2.bed | head -40
python3 single_copy_v2.py
echo DONE
