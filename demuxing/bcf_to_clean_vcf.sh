bcftools view -O v -o /mnt/sda/scz_meta_analysis/shin_nowakowski_cellstemcell_2025/genotypes/gtc2vcf/vcfs/from_idat.vcf /mnt/sda/scz_meta_analysis/shin_nowakowski_cellstemcell_2025/genotypes/gtc2vcf/vcfs/from_idat.bcf

bcftools annotate --threads=16 -Wtbi --rename-chrs /mnt/sda/scz_meta_analysis/shin_nowakowski_cellstemcell_2025/genotypes/gtc2vcf/vcfs/chr_mapping.txt -Oz -o /mnt/sda/scz_meta_analysis/shin_nowakowski_cellstemcell_2025/genotypes/gtc2vcf/vcfs/from_idat_chr_mapped.vcf.gz /mnt/sda/scz_meta_analysis/shin_nowakowski_cellstemcell_2025/genotypes/gtc2vcf/vcfs/from_idat.vcf
