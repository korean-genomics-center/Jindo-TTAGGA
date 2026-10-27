#!/bin/bash
set -uo pipefail
export PATH=${CONDA_PREFIX}/bin:$PATH
JINDO=${JINDO_ROOT}
BASE=$JINDO/Analysis/Population_20Jindo/Reanalysis_Combined
REF_DIR=$BASE/REF; LOG_DIR=$BASE/Logs
mkdir -p "$REF_DIR" "$LOG_DIR"
HAP1=$JINDO/Results/Genome_Assembly/Polishing/Reference_Fasta/Round7.SV_SNV_Polished_Fasta/PolishedASM.J495799_Child.Trio_Specific.Polish_R7.SV_SNV.hap1.fasta
HAP2=$JINDO/Results/Genome_Assembly/Polishing/Reference_Fasta/Round7.SV_SNV_Polished_Fasta/PolishedASM.J495799_Child.Trio_Specific.Polish_R7.SV_SNV.hap2.fasta
MITO=$JINDO/Results/Mitochondrial_Genome_Assembly/potential_contigs/atg000005l/atg000005l.mitogenome.rotated.fa
REF=$REF_DIR/jindo_male_combined.fasta
echo "===== Build combined MALE reference : $(date) on $(hostname) ====="
for f in "$HAP1" "$HAP2" "$MITO"; do [ -s "$f" ] || { echo "ERROR missing $f"; exit 1; }; done
echo "[1] Hap1 (autosomes + chrX)"; samtools faidx "$HAP1"
grep ">" "$HAP1" | sed 's/>//; s/ .*//' | grep -v -i -x "chrY" > "$REF_DIR/hap1_keep.list"
samtools faidx "$HAP1" $(cat "$REF_DIR/hap1_keep.list" | tr '\n' ' ') > "$REF_DIR/part_hap1.fa"
echo "    Hap1 seqs kept: $(grep -c '>' "$REF_DIR/part_hap1.fa")"
echo "[2] Hap2 chrY"; samtools faidx "$HAP2"; samtools faidx "$HAP2" chrY > "$REF_DIR/part_chrY.fa"
echo "[3] mito -> chrM"; echo ">chrM" > "$REF_DIR/part_chrM.fa"; grep -v ">" "$MITO" >> "$REF_DIR/part_chrM.fa"
echo "[4] concatenate"; cat "$REF_DIR/part_hap1.fa" "$REF_DIR/part_chrY.fa" "$REF_DIR/part_chrM.fa" > "$REF"
DUP=$(grep ">" "$REF" | sort | uniq -d); [ -n "$DUP" ] && { echo "ERROR dup: $DUP"; exit 1; }
echo "    total seqs: $(grep -c '>' "$REF") (expect 41)"; grep ">" "$REF" | tr '\n' ' '; echo ""
echo "[5] faidx"; samtools faidx "$REF"
echo "[6] bwa index (~30-60min)"; bwa index "$REF"
for c in chrX chrY chrM; do grep -q ">$c$" "$REF" || grep -q ">$c " "$REF" && echo "    OK $c" || { echo "ERROR $c missing"; exit 1; }; done
rm -f "$REF_DIR"/part_*.fa "$REF_DIR"/part_*.fai "$REF_DIR/hap1_keep.list"
echo "===== REFERENCE BUILD DONE: $REF : $(date) ====="
