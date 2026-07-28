#!/bin/bash

set -e
set -x

export PATH=${WORK_DIR}/miniconda3/envs/TTT/bin:$PATH

# Inputs
REF_HAP1=${JINDO_ROOT}/Results/Genome_Assembly/Final_Assembly/FinalASM.J495799_Child.Trio_Specific.hap1.fasta
REF_HAP2=${JINDO_ROOT}/Results/Genome_Assembly/Final_Assembly/FinalASM.J495799_Child.Trio_Specific.hap2.fasta
HIFI_DIR=${JINDO_ROOT}/Resources/SequencingData/HiFi/Data-X209SC25077815-Z02-F002
PBCPG=${WORK_DIR}/Tools/pb-CpG-tools/bin/aligned_bam_to_cpg_scores
PBCPG_MODEL=${WORK_DIR}/Tools/pb-CpG-tools/models/pileup_calling_model.v1.tflite

# Output
OUT_DIR=${JINDO_ROOT}/Analysis/Population_20Jindo/Methylation
TMP_DIR=${OUT_DIR}/tmp
mkdir -p ${OUT_DIR} ${TMP_DIR}

echo "================================================="
echo "5mC Methylation Calling (Hap1 + Hap2) | Started: $(date)"
echo "Running on: $(hostname)"
echo "================================================="

# Reference check
echo ""
echo "=== Reference files ==="
ls -la ${REF_HAP1} ${REF_HAP2}

# Step 1: HiFi BAM list
echo ""
echo "=== Step 1: HiFi BAM 파일 list ==="
HIFI_BAMS=$(find ${HIFI_DIR} -name "*.hifi_reads.bam" | sort)
echo "${HIFI_BAMS}" | nl
N_BAMS=$(echo "${HIFI_BAMS}" | wc -l)
echo "Total HiFi BAM 파일: ${N_BAMS}"

# Step 2: Merge HiFi BAMs
echo ""
echo "=== Step 2: Merge HiFi BAMs (preserving MM/ML tags) ==="
MERGED_BAM=${TMP_DIR}/J495799_hifi_merged.bam

if [ ! -f "${MERGED_BAM}" ]; then
    echo "Merging ${N_BAMS} HiFi BAMs..."
    samtools merge -@ 32 -f ${MERGED_BAM} ${HIFI_BAMS}
    echo "Merged BAM size: $(ls -la ${MERGED_BAM} | awk '{print $5}') bytes"
else
    echo "Merged BAM already exists, skipping"
fi

# Verify MM/ML tags preserved
echo ""
echo "=== Verify MM/ML tags in merged BAM ==="
samtools view ${MERGED_BAM} | head -1 | tr '\t' '\n' | grep -E "^(MM|ML):" | head -2

# Step 3: pbmm2 align to BOTH haplotypes
for HAP in hap1 hap2; do
    echo ""
    echo "============================================"
    echo "=== Step 3 [${HAP}]: pbmm2 alignment ==="
    echo "============================================"
    
    if [ "${HAP}" == "hap1" ]; then
        REF=${REF_HAP1}
    else
        REF=${REF_HAP2}
    fi
    
    ALIGNED_BAM=${TMP_DIR}/J495799_hifi_aligned_${HAP}.bam
    
    if [ ! -f "${ALIGNED_BAM}.bai" ]; then
        pbmm2 align \
            ${REF} \
            ${MERGED_BAM} \
            ${ALIGNED_BAM} \
            --preset CCS \
            --sort \
            -j 32 \
            --rg "@RG\tID:J495799\tSM:J495799\tLB:HiFi\tPL:PACBIO"
        samtools index -@ 8 ${ALIGNED_BAM}
    else
        echo "${HAP} aligned BAM already exists, skipping"
    fi
    
    ls -la ${ALIGNED_BAM}
    samtools flagstat ${ALIGNED_BAM} > ${OUT_DIR}/J495799_hifi_aligned_${HAP}_flagstat.txt
    cat ${OUT_DIR}/J495799_hifi_aligned_${HAP}_flagstat.txt
    
    # Verify MM/ML tags preserved
    echo ""
    echo "=== Verify MM/ML tags in ${HAP} aligned BAM ==="
    samtools view -F 4 ${ALIGNED_BAM} | head -1 | tr '\t' '\n' | grep -E "^(MM|ML):" | head -2
done

# Step 4: pb-CpG-tools for BOTH haplotypes
for HAP in hap1 hap2; do
    echo ""
    echo "============================================"
    echo "=== Step 4 [${HAP}]: 5mC CpG calling ==="
    echo "============================================"
    
    ALIGNED_BAM=${TMP_DIR}/J495799_hifi_aligned_${HAP}.bam
    OUT_PREFIX=${OUT_DIR}/J495799_${HAP}_cpg
    
    ${PBCPG} \
        --bam ${ALIGNED_BAM} \
        --output-prefix ${OUT_PREFIX} \
        --model ${PBCPG_MODEL} \
        --threads 32
    
    echo ""
    echo "=== Output files [${HAP}] ==="
    ls -la ${OUT_PREFIX}*
done

# Step 5: Summary statistics for BOTH haplotypes
echo ""
echo "============================================"
echo "=== Step 5: Methylation summary ==="
echo "============================================"

for HAP in hap1 hap2; do
    BED=${OUT_DIR}/J495799_${HAP}_cpg.combined.bed
    
    if [ -f "${BED}" ]; then
        echo ""
        echo "=== ${HAP} CpG methylation ==="
        echo "Total CpG sites: $(wc -l < ${BED})"
        
        awk '{print $4}' ${BED} | \
            awk 'BEGIN{n=0; sum=0; high=0; mid=0; low=0} 
                 {n++; sum+=$1; if($1>=80) high++; else if($1<=20) low++; else mid++} 
                 END {printf "Mean methylation: %.2f%%\nHigh (>=80%%): %d (%.1f%%)\nLow (<=20%%): %d (%.1f%%)\nIntermediate (20-80%%): %d (%.1f%%)\n", 
                      sum/n, high, high/n*100, low, low/n*100, mid, mid/n*100}'
    fi
done

echo ""
echo "================================================="
echo "5mC Methylation (Hap1+Hap2) | Finished: $(date)"
echo "================================================="
