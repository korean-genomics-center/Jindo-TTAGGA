#!/bin/bash

set -uo pipefail
ROOT=${JINDO_ROOT}
CB=${WORK_DIR}/miniconda3
NP=32
cd ${WORK_DIR}

echo "=== $(date) | $(hostname) | slots=${NSLOTS:-?} ==="

ft(){ for d in $CB/bin $CB/envs/*/bin /usr/bin /usr/local/bin; do
        [ -x "$d/$1" ] && { echo "$d/$1"; return; }; done; }
BCF=$(ft bcftools); SAM=$(ft samtools)
PLK=$(ft plink2); [ -z "$PLK" ] && PLK=$(ft plink)
echo "bcftools=${BCF:-NONE}"
echo "samtools=${SAM:-NONE}"
echo "plink   =${PLK:-NONE}"
echo

VCF=$ROOT/Analysis/Population_20Jindo/chrY_Analysis/Haplogroup/chrY_pass.vcf.gz
[ -f "$VCF" ] || VCF=$ROOT/Analysis/Population_20Jindo/chrY_Analysis/Haplogroup/chrY_clean.vcf.gz
echo "VCF=$VCF"; echo

echo "########## [13] chrY VCF 결측 ##########"
if [ -n "$BCF" ]; then
  $BCF query -l "$VCF" > samples.txt
  NS=$(wc -l < samples.txt); NV=$($BCF view -H "$VCF" | wc -l)
  echo "샘플 $NS / 사이트 $NV / genotype $((NS*NV))"
  $BCF query -f '[%GT\t]\n' "$VCF" | tr '\t' '\n' | grep -v '^$' \
    | sort | uniq -c | sort -rn > gt_dist.txt
  cat gt_dist.txt
  M=$(awk '$2=="./."||$2=="."{s+=$1}END{print s+0}' gt_dist.txt)
  awk -v m=$M -v t=$((NS*NV)) 'BEGIN{printf ">> 결측 %d (%.4f%%)  [0이 아니면 코멘트13은 오탐]\n",m,100*m/t}'
  $BCF query -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n' "$VCF" > gt_matrix.tsv
  echo ">> gt_matrix.tsv $(wc -l < gt_matrix.tsv) 행 생성"
else echo "SKIP: bcftools 없음"; fi
echo

echo "########## [6] chrY depth ##########"
if [ -n "$SAM" ]; then
  ls $ROOT/Analysis/Population_20Jindo/BAM/jindo_*_hap1.sorted.bam > bams.txt 2>/dev/null
  echo "BAM $(wc -l < bams.txt) 개, ${NP}병렬"
  cat > _d.sh <<'DEOF'
#!/bin/bash
B="$1"; SAM="$2"
S=$(basename "$B" | sed 's/^jindo_//; s/_hap1.sorted.bam$//')
f(){ $SAM depth -a $1 -r "$2" "$B" | awk '{s+=$3;n++}END{if(n)printf "%.3f",s/n; else printf "NA"}'; }
a=$(f "-Q 20" chrY:1-4000000); b=$(f "" chrY:1-4000000)
c=$(f "-Q 20" chrY);           d=$(f "" chrY)
printf "%s\t%s\t%s\t%s\t%s\n" "$S" "$a" "$b" "$c" "$d"
DEOF
  chmod +x _d.sh
  printf "sample\tQ20_0_4Mb\tnoQ_0_4Mb\tQ20_full\tnoQ_full\n" > depth_mapq.tsv
  xargs -a bams.txt -I{} -P $NP ./_d.sh {} "$SAM" >> depth_mapq.tsv
  echo ">> 코호트 평균"
  awk -F'\t' 'NR>1&&$2!="NA"{a+=$2;b+=$3;c+=$4;d+=$5;n++}END{
    printf "   0-4Mb : Q20 %.2fx | noQ %.2fx | noQ/Q20 %.3f\n",a/n,b/n,b/a;
    printf "   전체  : Q20 %.2fx | noQ %.2fx | noQ/Q20 %.3f\n",c/n,d/n,d/c;
    printf "   (n=%d)\n",n}' depth_mapq.tsv
else echo "SKIP: samtools 없음"; fi
echo

echo "########## [2] chrY QV 61.4 원인 ##########"
MQD=$ROOT/Results/Genome_Assembly/Assembly_Evaluation/Merqury_Evaluation/FinalASM.J495799_Child.merqury.k21
find "$MQD" -maxdepth 1 -name "*.qv" 2>/dev/null | while read f; do
  awk -v F="$(basename $f)" '$1=="chrY"{print "   "F": "$0}' "$f" 2>/dev/null; done | head -5
for H in hap1 hap2; do
  BED=$(find "$MQD" -maxdepth 1 -name "*${H}_only.bed" 2>/dev/null | head -1)
  [ -n "$BED" ] || continue
  N=$(awk '$1=="chrY"' "$BED" | wc -l)
  echo ">> $H : chrY 오류 k-mer $N 개"
  [ "$N" -gt 0 ] && awk '$1=="chrY"{print int($2/1000000)}' "$BED" | sort -n | uniq -c \
    | awk '{printf "     %2d-%2d Mb : %6d\n",$2,$2+1,$1}'
done
echo

echo "########## [3] kinship ##########"
PC=$ROOT/Results/Manuscript_Figures/data/jindo_roh/plink/per_chrom
if [ -n "$BCF" ] && [ -d "$PC" ]; then
  ls $PC/joint_chr*.vcf.gz 2>/dev/null | sort -V > vcflist.txt
  echo "VCF $(wc -l < vcflist.txt) 개 concat 중"
  $BCF concat -f vcflist.txt --threads $NP -Oz -o autosomes.vcf.gz && $BCF index -f autosomes.vcf.gz
  echo ">> $($BCF view -H autosomes.vcf.gz | wc -l) 사이트"
  if [ -n "$PLK" ]; then
    $PLK --vcf autosomes.vcf.gz --dog --maf 0.05 --geno 0.1 \
         --make-bed --out auto --threads $NP 2>&1 | tail -3
    $PLK --bfile auto --dog --genome --min 0.0 --out kinship --threads $NP 2>&1 | tail -3
    if [ -f kinship.genome ]; then
      echo ">> PLINK IBD 상위 30쌍 (PI_HAT = 2*kinship)"
      head -1 kinship.genome
      tail -n +2 kinship.genome | sort -k10,10gr | head -30
      tail -n +2 kinship.genome | awk '{p=$10; if(p>0.354)a++; else if(p>0.177)b++; else if(p>0.0884)c++; else d++}
        END{printf ">> 1촌급 %d / 2촌급 %d / 3촌급 %d / 무관 %d  (총 %d쌍)\n",a,b,c,d,a+b+c+d}'
    fi
  else
    echo "plink 없음 — python KING-robust 대체"
    $BCF query -l autosomes.vcf.gz > auto_samples.txt
    $BCF query -f '[%GT\t]\n' autosomes.vcf.gz | head -300000 > gt_auto.tsv
    python3 - <<'PY'
import itertools, collections
S=[l.strip() for l in open('auto_samples.txt')]
def c(g):
    g=g.replace('|','/')
    return {'0/0':0,'0/1':1,'1/0':1,'1/1':2}.get(g,-1)
rows=[]
for line in open('gt_auto.tsv'):
    v=[c(x) for x in line.rstrip('\n').split('\t') if x!='']
    if len(v)==len(S): rows.append(v)
print(f"   사용 사이트 {len(rows)}, 샘플 {len(S)}")
res=[]
for i,j in itertools.combinations(range(len(S)),2):
    hi=hj=nd=nb=0
    for v in rows:
        a,b=v[i],v[j]
        if a<0 or b<0: continue
        if a==1: hi+=1
        if b==1: hj+=1
        if (a==0 and b==2) or (a==2 and b==0): nd+=1
        if a==1 and b==1: nb+=1
    d=min(hi,hj)
    if d: res.append(((nb-2*nd)/(2*d), S[i], S[j]))
res.sort(reverse=True)
print(">> KING-robust 상위 30쌍")
for k,a,b in res[:30]: print(f"     {a:14s} {b:14s} {k:+.4f}")
cc=collections.Counter('1촌' if k>0.177 else '2촌' if k>0.0884 else '3촌' if k>0.0442 else '무관' for k,_,_ in res)
print(">>", dict(cc), f"/ 총 {len(res)}쌍")
PY
  fi
else echo "SKIP"; fi

echo "=== 완료 $(date) ==="
