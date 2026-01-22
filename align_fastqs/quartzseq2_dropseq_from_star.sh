#!/usr/bin/env bash
set -euo pipefail

# -----------------------------
# CONFIG
# -----------------------------
STAR_DIR=/mnt/sdb/scz_meta_analysis_processed/sawada_kato/star_outputs
OUT_DIR=/mnt/sdb/scz_meta_analysis_processed/sawada_kato/dropseq
GTF=/home/deepak/datasets/annotations/gencode.v43.nochr.gtf

THREADS=25

mkdir -p "${OUT_DIR}"

# -----------------------------
# SRRs to process
# -----------------------------
SRRS=(
  SRR7878531
  SRR7878532
  SRR7878533
  SRR7878534
  SRR7878535
  SRR7878536
  SRR7878537
  SRR7878538
)

# -----------------------------
# PROCESS
# -----------------------------
for SRR in "${SRRS[@]}"; do
  echo "=== Processing ${SRR} ==="

  STAR_BAM="${STAR_DIR}/${SRR}_Aligned.sortedByCoord.out.bam"

  [[ -f "${STAR_BAM}" ]] || {
    echo "Missing STAR BAM for ${SRR}, skipping"
    continue
  }

  cd "${OUT_DIR}"

  # 1. Queryname sort
  samtools sort -n -@ ${THREADS} \
    -o ${SRR}.qname.bam \
    "${STAR_BAM}"

  # 2. Tag UMI (R1 bases 1–8)
  drop-seq TagBamWithReadSequenceExtended \
    INPUT=${SRR}.qname.bam \
    OUTPUT=${SRR}.umi.bam \
    SUMMARY=${SRR}.umi.summary.txt \
    BASE_RANGE=1-8 \
    BASE_QUALITY=10 \
    BARCODED_READ=1 \
    DISCARD_READ=false \
    TAG_NAME=UB \
    NUM_BASES_BELOW_QUALITY=1

  # 3. Tag Cell Barcode (R1 bases 9–22)
  drop-seq TagBamWithReadSequenceExtended \
    INPUT=${SRR}.umi.bam \
    OUTPUT=${SRR}.cb.bam \
    SUMMARY=${SRR}.cb.summary.txt \
    BASE_RANGE=9-22 \
    BASE_QUALITY=10 \
    BARCODED_READ=1 \
    DISCARD_READ=false \
    TAG_NAME=CB \
    NUM_BASES_BELOW_QUALITY=1

  # 4. Assign genes
  drop-seq TagReadWithGeneExonFunction \
    INPUT=${SRR}.cb.bam \
    OUTPUT=${SRR}.gene.bam \
    ANNOTATIONS_FILE=${GTF} 
#    TAG=GE

  # 5. Digital expression
  drop-seq DigitalExpression \
    INPUT=${SRR}.gene.bam \
    OUTPUT=${SRR}.dge.txt.gz \
    SUMMARY=${SRR}.dge.summary.txt \
    NUM_CORE_BARCODES=3000

  echo "=== Done ${SRR} ==="
done

