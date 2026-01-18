cd /mnt/sda/scz_meta_analysis/shin_nowakowski_cellstemcell_2025/raw_data/ || exit 1

ls -d SRR*/ | sed 's|/||' | \
xargs -n 1 -P 1 -I {} \
bash -c '
  srr={}
  echo "Processing $srr"

  # Dump FASTQs
  fasterq-dump "$srr/$srr.sra" \
    --split-files \
    --force \
    -O ../fastq/ \
    -t ../tmp/ \
    --threads 4

  # Compress immediately
  pigz ../fastq/${srr}_*.fastq

  # Verify compression
  gzip -t ../fastq/${srr}_*.fastq.gz

  # Remove any leftover uncompressed FASTQs
  rm -f ../fastq/${srr}_*.fastq
'

