#pixi run --manifest-path /home/deepak/pixi_envs/cellbender/pixi.toml   plink   --file /mnt/sda/scz_meta_analysis/shin_nowakowski_cellstemcell_2025/genotypes/NowakowskiLabOnly_Genotyping2021QB3/20210817_GSA08_NowakowskiLabOnly   --recode vcf bgz   --out /mnt/sda/scz_meta_analysis/shin_nowakowski_cellstemcell_2025/genotypes/NowakowskiLabOnly/vcf

#!/bin/bash

# Paths
BAM_DIR=/mnt/sdb/scz_meta_analysis_processed/shin_nowakowski/star_outputs
VCF=/mnt/sda/scz_meta_analysis/shin_nowakowski_cellstemcell_2025/genotypes/vcf/NowakowskiLabOnly.vcf.gz
CELLSNP_BASE=/mnt/sdb/scz_meta_analysis_processed/shin_nowakowski/vireo_cellsnp
VIREO_BASE=/mnt/sdb/scz_meta_analysis_processed/shin_nowakowski/vireo_out

# Number of donors
NDONORS=44
THREADS=25

# List of SRR runs
SRRS=(SRR26424568 SRR26424569 SRR26424570 SRR26424571 SRR26424572 SRR26424573 \
      SRR26424574 SRR26424575 SRR26424576 SRR26424577 SRR26424578 SRR26424579 \
      SRR26424580 SRR26424581 SRR26424582 SRR26424583 SRR26424584 SRR26424585 \
      SRR26424586 SRR26424587 SRR26424588 SRR26424589 SRR26424590)

for SRR in "${SRRS[@]}"
do
    echo "=== Processing $SRR ==="

    BAM=${BAM_DIR}/${SRR}_Aligned.sortedByCoord.out.bam
    BARCODE=${BAM_DIR}/${SRR}_Solo.out/Gene/raw/barcodes.tsv
    CELLSNP_OUT=${CELLSNP_BASE}/${SRR}
    VIREO_OUT=${VIREO_BASE}/${SRR}

    mkdir -p $CELLSNP_OUT $VIREO_OUT

    # Run cellSNP-lite if not done
    if [ ! -f ${CELLSNP_OUT}/cellSNP.raw.count.h5 ]; then
        echo "Running cellSNP-lite for $SRR..."
        pixi run --manifest-path /home/deepak/pixi_envs/cellbender/pixi.toml cellsnp-lite \
            -s $BAM \
            -b $BARCODE \
            -R $VCF \
            -O $CELLSNP_OUT \
            --minMAF 0.1 \
            --minCOUNT 20 \
            -p $THREADS
    else
        echo "cellSNP output exists for $SRR, skipping..."
    fi

    # Run Vireo if not done
#    if [ ! -f ${VIREO_OUT}/donor_assignments.tsv ]; then
#        echo "Running Vireo for $SRR..."
#        pixi run --manifest-path /home/deepak/pixi_envs/cellbender/pixi.toml vireo -c $CELLSNP_OUT \
#              -d $VCF \
#              -N $NDONORS \
#              -o $VIREO_OUT
#    else
#        echo "Vireo output exists for $SRR, skipping..."
#    fi

    echo "=== Done with $SRR ==="
done

