cd /mnt/sda/scz_meta_analysis/notaras_colak_molecular_psychiatry_2021/raw_data/

ls -d SRR*/ | sed 's|/||' | \
xargs -n 1 -P 4 -I {} \
bash -c '
  echo "Processing {}"
  fasterq-dump {}/{}.sra \
    --split-files \
    --force \
    -O ../fastq/ \
    -t ../tmp/ \
    --threads 4
'

