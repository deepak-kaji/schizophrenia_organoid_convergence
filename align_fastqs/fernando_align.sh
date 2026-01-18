#!/usr/bin/env bash

#set -euo pipefail

FASTQ_DIR=/mnt/lacie/scz_meta_analysis/fernando_brennand_nature_2025/fastq
OUT_BASE=/mnt/sdb/scz_meta_analysis_processed/fernando_brennand
GENOME_DIR=/home/deepak/datasets/annotations/CELLRANGER_GENOME_DIR
WHITELIST=/home/deepak/datasets/annotations/3M-february-2018_TRU.txt

mkdir -p ${OUT_BASE}/star_outputs
mkdir -p ${OUT_BASE}/excluded
mkdir -p ${OUT_BASE}/logs

cd "$FASTQ_DIR"
pwd 

for fqI in *_1.fastq.gz; do

    sample=${fqI%_1.fastq.gz}
    fqCDNA=${sample}_3.fastq.gz
    fqCB=${sample}_2.fastq.gz

    [[ -f "$fqCDNA" && -f "$fqCB" ]] || {
        echo "Missing FASTQs for $sample, skipping"
        continue
    }

    echo "[$sample] Assuming 10x v3 chemistry"

    STAR \
        --runThreadN 25 \
        --genomeDir "$GENOME_DIR" \
        --readFilesIn "$fqCDNA" "$fqCB" \
        --readFilesCommand zcat \
        --soloType CB_UMI_Simple \
        --soloCBwhitelist "$WHITELIST" \
	--soloCBstart 1 \
	--soloCBlen 16 \
	--soloUMIstart 17 \
        --soloUMIlen 12 \
        --soloFeatures Gene Velocyto \
        --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts \
        --soloUMIfiltering MultiGeneUMI_CR \
        --soloUMIdedup 1MM_CR \
        --clipAdapterType CellRanger4 \
        --outFilterScoreMin 30 \
        --soloBarcodeReadLength 101 \
	--clip3pNbases 0 28 \
	--outSAMtype BAM SortedByCoordinate \
        --outSAMattributes CR UR CY UY CB UB \
        --outFileNamePrefix ${OUT_BASE}/star_outputs/${sample}_ \
        --outStd Log Progress \
        > ${OUT_BASE}/logs/${sample}.log 2>&1

    echo "[$sample] STARsolo finished"

done

