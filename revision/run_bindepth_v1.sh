#!/bin/bash
#$ -S /bin/bash
#$ -N bindep
#$ -cwd
#$ -q all.q
#$ -pe smp 4
#$ -t 1-3
#$ -o logs/bindep.$TASK_ID.log
#$ -j y
set -eu
source ${WORK_DIR}/miniconda3/etc/profile.d/conda.sh
set +u; conda activate TTT; set -u
A=${JINDO_ROOT}/Results/Genome_Assembly/Assembly_Evaluation/GCI_Evaluation/Align.Winnowmap
BINS=(p 0 a); b=${BINS[$((SGE_TASK_ID-1))]}
B=$A/Winnowmap_Align.Ref.FinalASM.J495799_Child.hap1.QRY.PacBio_HiFi_ccs_converted.Filter_L1000_QV20.triobin_split.triobin_${b}.bam
samtools view -@ 4 -F 0x904 $B | awk -v b=$b '
BEGIN{n=split("chr10:4000001:7000000 chr17:15000001:20000000 chr1:80000001:83000000",R," "); for(i=1;i<=n;i++){split(R[i],x,":"); C[i]=x[1]; S[i]=x[2]; E[i]=x[3]}}
{ c=$3; p=$4; cg=$6; L=0
  while(match(cg,/[0-9]+[MDN=X]/)){ L+=substr(cg,RSTART,RLENGTH-1); cg=substr(cg,RSTART+RLENGTH) }
  e=p+L-1
  for(i=1;i<=n;i++) if(c==C[i] && e>=S[i] && p<=E[i]){ a=(p>S[i]?p:S[i]); z=(e<E[i]?e:E[i]); T[i]+=z-a+1 }
}
END{for(i=1;i<=n;i++) printf "bin %s\t%s:%d-%d\tmeandepth %.1f\n",b,C[i],S[i],E[i],T[i]/(E[i]-S[i]+1)}' > bindepth_${b}_v1.txt
echo "DONE $b"
