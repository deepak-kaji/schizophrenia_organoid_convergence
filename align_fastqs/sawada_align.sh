#!/usr/bin/env bash

set -euo pipefail

FASTQ_DIR=/mnt/sda/scz_meta_analysis/sawada_kato_molecular_psych_2020/fastq
OUT_BASE=/mnt/sdb/scz_meta_analysis_processed/sawada_kato
GENOME_DIR=/home/deepak/datasets/annotations/CELLRANGER_GENOME_DIR

mkdir -p ${OUT_BASE}/star_outputs
mkdir -p ${OUT_BASE}/logs

cd "$FASTQ_DIR"

for sample in SRR7878531 SRR7878532 SRR7878533 SRR7878534 SRR7878535 SRR7878536 SRR7878537 SRR7878538; do

    r1="${sample}_1.fastq.gz"
    r2="${sample}_2.fastq.gz" 

    [[ -f "$r2" ]] || {
        echo "Missing R2 for $sample, skipping"
        continue
    }

    echo "[$sample] Quartz-Seq2 alignment"

    STAR \
        --runThreadN 25 \
        --genomeDir "$GENOME_DIR" \
        --readFilesIn "$r2" "$r1" \
        --readFilesCommand zcat \
	--soloType CB_UMI_Simple \
	--soloCBwhitelist /mnt/sda/scz_meta_analysis/sawada_kato_molecular_psych_2020/quartzseq2_barcodes.txt \
	--soloCBlen 14 \
	--soloUMIlen 8 \
	--soloFeatures Gene Velocyto \
	--soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts \
	--soloUMIdedup 1MM_CR \
        --outFilterScoreMin 30 \
        --outSAMtype BAM SortedByCoordinate \
        --outSAMattributes CR UR CY UY CB UB \
	--outTmpDir "${OUT_BASE}/tmp/${sample}" \
        --outTmpKeep None \
        --outFileNamePrefix "${OUT_BASE}/star_outputs/${sample}_" \
        --outStd Log Progress \ 
        > ${OUT_BASE}/logs/${sample}.log 2>&1

    echo "STARsolo finished"

done

