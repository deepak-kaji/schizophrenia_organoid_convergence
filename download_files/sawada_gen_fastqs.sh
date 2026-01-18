cd /mnt/sda/scz_meta_analysis/sawada_kato_molecular_psych_2020/raw_data/

ls -d SRR*/ | sed 's|/||' | \
xargs -n 1 -P 4 -I {} \
bash -c '
  echo "Processing {}"
  fasterq-dump {}/{}.sra \
    --split-files \
    --force \
    --include-technical \
    -O ../fastq/ \
    -t ../tmp/ \
    --threads 4
'

