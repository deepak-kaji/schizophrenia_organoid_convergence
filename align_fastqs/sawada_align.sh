#!/usr/bin/env bash

set -euo pipefail

FASTQ_DIR=/mnt/sda/scz_meta_analysis/sawada_kato_molecular_psych_2020/fastq
OUT_BASE=/mnt/sdb/scz_meta_analysis_processed/sawada_kato
GENOME_DIR=/home/deepak/datasets/annotations/CELLRANGER_GENOME_DIR

mkdir -p ${OUT_BASE}/star_outputs
mkdir -p ${OUT_BASE}/logs

cd "$FASTQ_DIR"

for r1 in *_1.fastq.gz; do
    sample=${r1%_1.fastq.gz}
    r2=${sample}_2.fastq.gz

    [[ -f "$r2" ]] || {
        echo "Missing R2 for $sample, skipping"
        continue
    }

    echo "[$sample] Quartz-Seq2 alignment"

    STAR \
        --runThreadN 32 \
        --genomeDir "$GENOME_DIR" \
        --readFilesIn "$r1" "$r2" \
        --readFilesCommand zcat \
        --sjdbGTFfile /home/deepak/datasets/annotations/gencode.v43.nochr.gtf \
        --outFileNamePrefix "${OUT_BASE}/star_outputs/${sample}_" \
        --outTmpDir "${OUT_BASE}/star_tmp/${sample}" \
        --outTmpKeep None \
        --outSAMtype BAM SortedByCoordinate \
        --outFilterMultimapNmax 20 \
        --outFilterMismatchNoverReadLmax 0.04 \
        --alignSJoverhangMin 8 \
        --alignSJDBoverhangMin 1 \
        --outSAMattributes NH HI AS nM MD \
        --limitBAMsortRAM 300000000000

done

