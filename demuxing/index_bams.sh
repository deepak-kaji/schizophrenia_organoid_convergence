#!/bin/bash

BAM_DIR=/mnt/sdb/scz_meta_analysis_processed/shin_nowakowski/star_outputs
THREADS=25

for BAM in ${BAM_DIR}/*_Aligned.sortedByCoord.out.bam
do
    echo "Checking $BAM"

    # Check sort order
    SORT_ORDER=$(samtools view -H "$BAM" | grep '^@HD' | grep -o 'SO:[^[:space:]]*')

    if [[ "$SORT_ORDER" != "SO:coordinate" ]]; then
        echo "  ❌ BAM is not coordinate sorted. Skipping."
        continue
    fi

    # Check if index exists
    if [[ -f "${BAM}.bai" || -f "${BAM%.bam}.bai" ]]; then
        echo "  Index exists."
    else
        echo "  Index missing. Creating index..."
        samtools index -@ $THREADS "$BAM"
    fi

done

