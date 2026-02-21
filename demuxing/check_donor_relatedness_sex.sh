# Step 1: Convert VCF to PLINK BED, handling underscores in sample IDs
pixi run --manifest-path ~/pixi_envs/cellbender/ plink \
    --vcf /mnt/sda/scz_meta_analysis/shin_nowakowski_cellstemcell_2025/genotypes/gtc2vcf/vcfs/from_idat_chr_mapped.vcf.gz \
    --make-bed \
    --double-id \
    --out /mnt/sda/scz_meta_analysis/shin_nowakowski_cellstemcell_2025/genotypes/gtc2vcf/vcfs/donors

# Step 2: Predict sex from X chromosome heterozygosity
pixi run --manifest-path ~/pixi_envs/cellbender/ plink \
    --bfile /mnt/sda/scz_meta_analysis/shin_nowakowski_cellstemcell_2025/genotypes/gtc2vcf/vcfs/donors \
    --check-sex \
    --out /mnt/sda/scz_meta_analysis/shin_nowakowski_cellstemcell_2025/genotypes/gtc2vcf/vcfs/donors_sex

# Step 3: Compute pairwise relatedness (IBD / PI_HAT)
pixi run --manifest-path ~/pixi_envs/cellbender/ plink \
    --bfile /mnt/sda/scz_meta_analysis/shin_nowakowski_cellstemcell_2025/genotypes/gtc2vcf/vcfs/donors \
    --genome \
    --out /mnt/sda/scz_meta_analysis/shin_nowakowski_cellstemcell_2025/genotypes/gtc2vcf/vcfs/donors_relatedness


