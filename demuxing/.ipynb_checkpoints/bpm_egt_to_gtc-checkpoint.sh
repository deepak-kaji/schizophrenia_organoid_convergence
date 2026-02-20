bpm_manifest_file=/mnt/sda/scz_meta_analysis/shin_nowakowski_cellstemcell_2025/genotypes/gtc2vcf/GSA-24v3-0_A2.bpm
path_to_idat_folder=/mnt/sda/scz_meta_analysis/shin_nowakowski_cellstemcell_2025/genotypes/gtc2vcf/idats
egt_cluster_file=/mnt/sda/scz_meta_analysis/shin_nowakowski_cellstemcell_2025/genotypes/gtc2vcf/GSA-24v3-0_A1_ClusterFile.egt
ref="/home/deepak/GRCh38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna"
#ref=/home/deepak/datasets/annotations/refdata-cellranger-GRCh38-3.0.0/fasta/genome.fa   $HOME/GRCh37/human_g1k_v37.fasta"
out_prefix=/mnt/sda/scz_meta_analysis/shin_nowakowski_cellstemcell_2025/genotypes/gtc2vcf/gtcs

bcftools +idat2gtc \
  --bpm $bpm_manifest_file \
  --egt $egt_cluster_file \
  --idats $path_to_idat_folder \
  --output $out_prefix
