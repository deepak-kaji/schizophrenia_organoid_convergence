cd /mnt/sda/scz_meta_analysis/rao_gogos_nature_comm_2025/

cat remaining_srr_list.txt | xargs -n 1 -P 1 -I {} bash -c '
  echo "Processing {}"

  fasterq-dump {} \
    --split-files \
    --force \
    -O ./fastq/ \
    -t ./tmp/ \
    -v \
    --threads 4

  pigz ./fastq/{}_*.fastq
'

