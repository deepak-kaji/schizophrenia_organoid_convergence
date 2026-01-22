#!/usr/bin/env bash

#set -euo pipefail

FASTQ_DIR=/mnt/sda/scz_meta_analysis/shin_nowakowski_cellstemcell_2025/fastq
OUT_BASE=/mnt/sdb/scz_meta_analysis_processed/shin_nowakowski
GENOME_DIR=/home/deepak/datasets/annotations/CELLRANGER_GENOME_DIR
WHITELIST=/home/deepak/datasets/annotations/3M-february-2018_TRU.txt

mkdir -p ${OUT_BASE}/star_outputs
mkdir -p ${OUT_BASE}/excluded
mkdir -p ${OUT_BASE}/logs

cd "$FASTQ_DIR"

for fq1 in *_1.fastq.gz; do

    sample=${fq1%_1.fastq.gz}
    fq2=${sample}_2.fastq.gz

    # sanity check
    [[ -f "$fq2" ]] || { echo "Missing $fq2, skipping"; continue; }
     
    echo $fq1
    # detect R1 length from first read
    R1_LEN=$(zcat "$fq1" | awk 'NR==2 { print length($0); exit }' || true)

    echo "[$sample] R1 length = $R1_LEN"

    # decide chemistry
    if [[ "$R1_LEN" -eq 28 ]]; then
        CHEM="10xv3"
        CB_LEN=16
        UMI_LEN=12
    elif [[ "$R1_LEN" -eq 26 ]]; then
        CHEM="10xv2"
        CB_LEN=16
        UMI_LEN=10
    else
        echo "[$sample] Non-10x or unsupported chemistry (R1=$R1_LEN). Excluding."
        mv "$fq1" "$fq2" ${OUT_BASE}/excluded/
        continue
    fi

    echo "[$sample] Detected ${CHEM}, running STARsolo"

    STAR \
        --runThreadN 25 \
        --genomeDir "$GENOME_DIR" \
        --readFilesIn "$fq2" "$fq1" \
        --readFilesCommand zcat \
        --soloType CB_UMI_Simple \
        --soloCBwhitelist "$WHITELIST" \
        --soloCBlen "$CB_LEN" \
        --soloUMIlen "$UMI_LEN" \
        --soloFeatures Gene Velocyto \
        --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts \
        --soloUMIfiltering MultiGeneUMI_CR \
        --soloUMIdedup 1MM_CR \
        --clipAdapterType CellRanger4 \
        --outFilterScoreMin 30 \
        --outSAMtype BAM SortedByCoordinate \
        --outSAMattributes CR UR CY UY CB UB \
        --outFileNamePrefix ${OUT_BASE}/star_outputs/${sample}_ \
        --outStd Log Progress \
	> ${OUT_BASE}/logs/${sample}.log 2>&1

    echo "STARsolo finished"

done

