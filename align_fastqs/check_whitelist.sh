NREADS=1000000

CR_BARCODES=~/programs/cellranger/cellranger-10.0.0/lib/python/cellranger/barcodes
WH_2014=$CR_BARCODES/737K-april-2014_rc.txt
WH_2016=$CR_BARCODES/737K-august-2016.txt
WH_2018=~/datasets/annotations/3M-february-2018_TRU.txt
WH_2023=~/datasets/annotations/3M-3pgex-may-2023_TRU.txt

# cd to the directory worthy of checking

echo "Sample,CB_len,Total_unique_CBs,Matches_2014,Matches_2016,Matches_2018,Matches_2023"

for fq in *_1*.fastq.gz
do
  for CB_LEN in 16
  do
    zcat "$fq" \
      | awk 'NR%4==2 {print substr($0,1,'"$CB_LEN"')}' \
      | head -n "$NREADS" \
      | sort -u \
      > /tmp/cb_test.$$

    TOTAL=$(wc -l < /tmp/cb_test.$$)
    M2014=$(grep -F -x -f /tmp/cb_test.$$ "$WH_2014" | wc -l)
    M2016=$(grep -F -x -f /tmp/cb_test.$$ "$WH_2016" | wc -l)
    M2018=$(grep -F -x -f /tmp/cb_test.$$ "$WH_2018" | wc -l)
    M2023=$(grep -F -x -f /tmp/cb_test.$$ "$WH_2023" | wc -l)

    echo "$(basename "$fq"),$CB_LEN,$TOTAL,$M2014,$M2016,$M2018,$M2023"
  done
done

rm -f /tmp/cb_test.$$

