#!/bin/bash

###############################################################################
# MAGMA SCZ PGC3 pipeline
#
# Goal:
# Map schizophrenia GWAS summary statistics (PGC3 SCZ 2022)
# to genes using MAGMA v1.10 and generate gene-level SCZ association results.
#
# Genome build:
# GRCh37 / hg19
#
# MAGMA:
# v1.10
#
# GWAS:
# PGC3_SCZ_wave3.primary.autosome.public.v3.vcf.tsv.gz
#
###############################################################################

set -e

###############################################################################
# 0. Set working directory
###############################################################################

cd ~/programs/MAGMA


###############################################################################
# 1. Reference files
###############################################################################

# Gene coordinates:
# Downloaded from MAGMA website
# Build 37 gene locations
#
# reference/NCBI37/NCBI37.3.gene.loc

# LD reference:
# 1000 Genomes European panel
#
# reference/LD/g1000_eur.bed
# reference/LD/g1000_eur.bim
# reference/LD/g1000_eur.fam
# reference/LD/g1000_eur.synonyms


###############################################################################
# 2. Input GWAS
###############################################################################

GWAS=/mnt/sdb/scz_meta_analysis_processed/psychad_mapping/gwas_files/PGC3_SCZ_wave3.primary.autosome.public.v3.vcf.tsv.gz


###############################################################################
# 3. Create SNP location file
#
# MAGMA annotation requires:
#
# SNP_ID    chromosome    position
#
# This includes ALL SNPs, not only significant SNPs.
###############################################################################

zcat ${GWAS} \
| awk 'BEGIN{OFS="\t"} $1 !~ /^#/ && $1!="CHROM" {print $2,$1,$3}' \
> gwas/scz2022/SCZ_PGC3.snp.loc


###############################################################################
# 4. Check SNP location file
###############################################################################

head gwas/scz2022/SCZ_PGC3.snp.loc

wc -l gwas/scz2022/SCZ_PGC3.snp.loc


###############################################################################
# 5. Create SNP p-value file
#
# MAGMA requires:
#
# SNP_ID    P_VALUE
#
###############################################################################

zcat ${GWAS} \
| awk 'BEGIN{OFS="\t"} $1 !~ /^#/ && $1!="CHROM" {print $2,$11}' \
> gwas/scz2022/SCZ_PGC3.pval


###############################################################################
# 6. Check p-value file
###############################################################################

head gwas/scz2022/SCZ_PGC3.pval

wc -l gwas/scz2022/SCZ_PGC3.pval


###############################################################################
# 7. Initial annotation test
#
# Purpose:
# Confirm that SNP coordinates and gene coordinates match.
#
# This uses gene-body only annotation.
###############################################################################

./magma \
--annotate \
--snp-loc gwas/scz2022/SCZ_PGC3.snp.loc \
--gene-loc reference/NCBI37/NCBI37.3.gene.loc \
--out results/SCZ_PGC3


###############################################################################
# 8. Final annotation
#
# PGC3 SCZ paper used:
#
# 35 kb upstream
# 10 kb downstream
#
# This captures nearby regulatory regions.
###############################################################################

./magma \
--annotate window=35,10 \
--snp-loc gwas/scz2022/SCZ_PGC3.snp.loc \
--gene-loc reference/NCBI37/NCBI37.3.gene.loc \
--out results/SCZ_PGC3_window35_10


###############################################################################
# 9. Output from annotation
#
# Main file:
#
# results/SCZ_PGC3_window35_10.genes.annot
#
# Contains:
#
# Entrez_gene_ID
# genomic coordinates
# SNPs assigned to each gene
#
###############################################################################

echo "Annotation complete."

# Find N

zcat /mnt/sdb/scz_meta_analysis_processed/psychad_mapping/gwas_files/PGC3_SCZ_wave3.primary.autosome.public.v3.vcf.tsv.gz \ | awk '$1 !~ /^#/ && $1!="CHROM" {sum += 2*$16; n++} END {print "Mean effective N =",sum/n}'

# Use N to calculate MAGMA gene significance scores #

./magma \
--bfile reference/LD/g1000_eur \
--pval gwas/scz2022/SCZ_PGC3.pval N=168175 \
--gene-annot results/SCZ_PGC3_window35_10.genes.annot \
--out results/SCZ_PGC3_gene_results
