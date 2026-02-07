#!/usr/bin/env bash

#set -euo pipefail

FASTQ_DIR=/mnt/sda/scz_meta_analysis/sebastian_pak_nat_comm_2023/fastq
OUT_BASE=/mnt/sdb/scz_meta_analysis_processed/sebastian_pak
GENOME_DIR=/home/deepak/datasets/annotations/CELLRANGER_GENOME_DIR
WHITELIST=/home/deepak/datasets/annotations/3M-february-2018_TRU.txt

mkdir -p ${OUT_BASE}/star_outputs
mkdir -p ${OUT_BASE}/excluded
mkdir -p ${OUT_BASE}/logs

cd "$FASTQ_DIR"

for fq1 in *_2.fastq.gz; do

    sample=${fq1%_2.fastq.gz}
    fq2=${sample}_3.fastq.gz

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
        --readFilesCommand "pigz -dc -p 8" \
        --soloType CB_UMI_Simple \
        --soloCBwhitelist "$WHITELIST" \
        --soloCBlen "$CB_LEN" \
        --soloUMIlen "$UMI_LEN" \
        --soloFeatures Gene Velocyto \
	--soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts \
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

SECOND_FASTQ_DIR=/mnt/sda/scz_meta_analysis/sebastian_pak_nat_comm_2023/second_round_fastqs

cd "SECOND_FASTQ_DIR"

for sample in SRR23984432 SRR23984433 SRR23984442 SRR23984443 SRR23984444 SRR23984445 SRR23984446 SRR23984447 SRR23984494 SRR23984495 SRR23984496 SRR23984497
do
    fq1=${sample}_3.fastq.gz
    fq2=${sample}_4.fastq.gz

    # sanity check
    [[ -f "$fq1" ]] || { echo "Missing $fq1, skipping"; continue; }
    [[ -f "$fq2" ]] || { echo "Missing $fq2, skipping"; continue; }
     
    echo "$sample"
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

