#!/bin/bash
#$ -S /bin/bash
#$ -N fchk
#$ -cwd
#$ -q all.q
#$ -pe smp 4
#$ -o logs/final_checks_v1.log
#$ -j y
set -u
source ${WORK_DIR}/miniconda3/etc/profile.d/conda.sh
set +u; conda activate TTT; set -u
J=${JINDO_ROOT}
echo "########## [1] N count per sequence (SUBMISSION_FINAL)"
for h in hap1 hap2; do
  awk -v h=$h '/^>/{if(n!="")printf "%s\t%s\t%d\t%d\n",h,n,len,N; n=substr($1,2); len=0; N=0; next} {len+=length($0); N+=gsub(/[Nn]/,"")} END{printf "%s\t%s\t%d\t%d\n",h,n,len,N}' $J/Resources/SUBMISSION_FINAL/assembly/Jindo1-G-TTAGGA.$h.fasta > nfinal_$h.tsv
  echo "  $h: sequences $(wc -l < nfinal_$h.tsv), total N $(awk '{s+=$4}END{print s}' nfinal_$h.tsv), sequences with N>0: $(awk '$4>0' nfinal_$h.tsv | wc -l)"
  awk '$4>0{print "    "$0}' nfinal_$h.tsv | head
done
echo "########## [2] telomere data"
T=$J/Results/Manuscript_Figures/data/telomere/final
ls -l $T
for f in $T/*; do echo "== $(basename $f)"; head -4 "$f"; echo "   lines: $(wc -l < "$f")"; done 2>/dev/null | head -60
echo "########## [3] Merqury"
M=$J/Results/Genome_Assembly/Assembly_Evaluation/Merqury_Evaluation/FinalASM.J495799_Child.merqury.k21
ls $M | head -40
for f in $M/*.completeness.stats $M/*.qv; do [ -f "$f" ] && { echo "== $(basename $f)"; grep -v "^chr" "$f" | head -6; }; done
for f in $M/*phased_block.stats $M/*switch* $M/*hamming*; do [ -f "$f" ] && { echo "== $(basename $f)"; head -6 "$f"; }; done
find $J/Results/Genome_Assembly/Assembly_Evaluation -maxdepth 3 \( -iname "*trioeval*" -o -iname "*hamming*" -o -iname "*switch*" \) 2>/dev/null | head
echo "########## [4] rDNA array depth (HiFi, all MAPQ; Round7 coords)"
A=$J/Results/Genome_Assembly/Assembly_Evaluation/GCI_Evaluation/Align.Winnowmap
B1=$(ls $A/Winnowmap_Align.Ref.FinalASM.J495799_Child.hap1.QRY.PacBio_HiFi*triobin_p0a.sorted.bam | head -1)
B2=$(ls $A/Winnowmap_Align.Ref.FinalASM.J495799_Child.hap2.QRY.PacBio_HiFi*triobin_m0a.sorted.bam | head -1)
md(){ samtools depth -a -r $2 $1 | awk '{s+=$3;n++}END{printf "%.1f",s/n}'; }
echo "  hap1 chr32 rDNA 0.15-0.45 Mb: $(md $B1 chr32:150001-450000)x | flank chr32 1-2 Mb: $(md $B1 chr32:1000001-2000000)x"
echo "  hap2 chr17 rDNA 66.40-66.70 Mb: $(md $B2 chr17:66400001-66700000)x | flank chr17 60-61 Mb: $(md $B2 chr17:60000001-61000000)x"
echo "  hap2 chrY  rDNA 0.15-0.45 Mb: $(md $B2 chrY:150001-450000)x | flank chrY 1-2 Mb: $(md $B2 chrY:1000001-2000000)x"
echo DONE
