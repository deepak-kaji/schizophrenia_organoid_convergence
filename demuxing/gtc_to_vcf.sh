bpm_manifest_file=/mnt/sda/scz_meta_analysis/shin_nowakowski_cellstemcell_2025/genotypes/gtc2vcf/GSA-24v3-0_A2.bpm
egt_cluster_file=/mnt/sda/scz_meta_analysis/shin_nowakowski_cellstemcell_2025/genotypes/gtc2vcf/GSA-24v3-0_A1_ClusterFile.egt
path_to_gtc_folder=/mnt/sda/scz_meta_analysis/shin_nowakowski_cellstemcell_2025/genotypes/gtc2vcf/gtcs/
ref="/home/deepak/GRCh38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna" # or ref="$HOME/GRCh37/human_g1k_v37.fasta"
out_prefix=/mnt/sda/scz_meta_analysis/shin_nowakowski_cellstemcell_2025/genotypes/gtc2vcf/vcfs/from_idat
bcftools +gtc2vcf \
  --no-version -Ou \
  --bpm $bpm_manifest_file \
  --egt $egt_cluster_file \
  --gtcs $path_to_gtc_folder \
  --fasta-ref $ref \
  --extra $out_prefix.tsv | \
  bcftools sort -Ou -T ./bcftools. | \
  bcftools norm --no-version -o $out_prefix.bcf -Ob -c x -f $ref --write-index
