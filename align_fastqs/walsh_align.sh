#!/usr/bin/env bash

#set -euo pipefail

FASTQ_DIR=/mnt/sda/scz_meta_analysis/walsh_studer_neuron_2025/fastq
OUT_BASE=/mnt/sdb/scz_meta_analysis_processed/walsh_studer
GENOME_DIR=/home/deepak/datasets/annotations/CELLRANGER_GENOME_DIR
WHITELIST=/home/deepak/datasets/annotations/3M-february-2018_TRU.txt

mkdir -p ${OUT_BASE}/star_outputs
mkdir -p ${OUT_BASE}/excluded
mkdir -p ${OUT_BASE}/logs

cd "$FASTQ_DIR"

for fq1 in *_3.fastq.gz; do

    sample=${fq1%_3.fastq.gz}
    fq2=${sample}_4.fastq.gz

    # sanity check
    [[ -f "$fq2" ]] || { echo "Missing $fq2, skipping"; continue; }
     
    # detect R1 length from first read
    R1_LEN=$(zcat "$fq1" | awk 'NR==2 { print length($0); exit }' || true)

    echo "[$sample] R1 length = $R1_LEN"

    echo "running STARsolo"

    STAR \
        --runThreadN 25 \
        --genomeDir "$GENOME_DIR" \
        --readFilesIn "$fq2" "$fq1" \
        --readFilesCommand zcat \
        --soloType CB_UMI_Simple \
        --soloCBwhitelist "$WHITELIST" \
	--soloCBstart 1 \
        --soloCBlen 16 \
	--soloUMIstart 17 \
        --soloUMIlen 12 \
        --soloFeatures Gene Velocyto \
	--soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts \
	--soloUMIdedup 1MM_CR \
        --clipAdapterType CellRanger4 \
        --outFilterScoreMin 30 \
	--clip3pNbases 0 28 \
	--soloBarcodeReadLength 29 \
        --outSAMtype BAM SortedByCoordinate \
        --outSAMattributes CR UR CY UY CB UB \
        --outFileNamePrefix ${OUT_BASE}/star_outputs/${sample}_ \
        --outStd Log Progress \
	> ${OUT_BASE}/logs/${sample}.log 2>&1

    echo "STARsolo finished"

done

shopt -s nullglob

for fq1 in *_1.fastq.gz; do
    sample=${fq1%_1.fastq.gz}

    fq2=${sample}_2.fastq.gz

    # must have R2
    [[ -f "$fq2" ]] || continue

    # skip 4-read libraries (already handled elsewhere)
    [[ -f "${sample}_3.fastq.gz" ]] && continue

    # sanity check: barcode read length (expect ~16+12 = 28–29)
    R1_LEN=$(zcat "$fq1" | awk 'NR==2 { print length($0); exit }')

    echo "[$sample] 2-read library detected (R1=$R1_LEN bp)"
    echo "running STARsolo"

    STAR \
        --runThreadN 25 \
        --genomeDir "$GENOME_DIR" \
        --readFilesIn "$fq2" "$fq1" \
        --readFilesCommand zcat \
        --soloType CB_UMI_Simple \
        --soloCBwhitelist "$WHITELIST" \
        --soloCBstart 1 \
        --soloCBlen 16 \
        --soloUMIstart 17 \
        --soloUMIlen 12 \
        --soloFeatures Gene Velocyto \
        --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts \
        --soloUMIdedup 1MM_CR \
        --clipAdapterType CellRanger4 \
        --outFilterScoreMin 30 \
	--clip3pNbases 0 28 \
	--soloBarcodeReadLength 29 \
        --outSAMtype BAM SortedByCoordinate \
        --outSAMattributes CR UR CY UY CB UB \
        --outFileNamePrefix ${OUT_BASE}/star_outputs/${sample}_ \
        --outStd Log Progress \
        > ${OUT_BASE}/logs/${sample}.log 2>&1

    echo "STARsolo finished for $sample"
done


