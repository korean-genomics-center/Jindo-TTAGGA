#!/bin/bash
set -eu
J=${JINDO_ROOT}
MM=$J/Resources/Tools/minimap2/minimap2_v2.30_r1287/minimap2/minimap2
SY=syri
ROS=$J/Results/Genome_Assembly/Assembly_Scaffolding/Hifiasm/J495799_Child.ont_L50kb_QV15.hifi_1kb_QV20.Trio_Specific.default_option/Reference_Fasta.Major_only/GCA_014441545.1/GCA_014441545.1_ROS_Cfam_1.0_genomic.major_only.fna
H=(hap1 hap2); L=(TTAGGA-M TTAGGA-P)
h=${H[$((SGE_TASK_ID-1))]}; l=${L[$((SGE_TASK_ID-1))]}
REF=$J/Results/Genome_Assembly/Final_Assembly/FinalASM.J495799_Child.Trio_Specific.$h.fasta
D=$PWD/$l.R7orient.ROS_Cfam_1.0; mkdir -p $D
P=$D/Minimap2.Ref_${l}.Qry_ROS_Cfam_1.0.asm5c_eqx_R7orient.paf
[ -s "$P" ] || { $MM -t 32 -c --eqx -x asm5 $REF $ROS > $P.tmp && mv $P.tmp $P; }
$SY -c $P -r $REF -q $ROS -F P --dir $D --prefix syri_v3.ref_${l}.qry_ROS_Cfam_1.0 --nc 16 --no-chrmatch
S=$D/syri_v3.ref_${l}.qry_ROS_Cfam_1.0syri.summary
if [ -s "$S" ]; then echo "DONE $l"; cat $S; else echo "FAILED $l: no summary"; tail -5 $D/syri_v3.*log; exit 1; fi
