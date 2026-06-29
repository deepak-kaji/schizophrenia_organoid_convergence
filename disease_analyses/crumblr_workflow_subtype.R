library(dplyr)
library(zellkonverter)
library(BiocParallel)
library(SingleCellExperiment)
library(dreamlet)
library(crumblr)

CSV_PATH <- '/mnt/sdb/scz_meta_analysis_processed/dge_signatures/dreamlet_dges/disease_analyses/'

BPPARAM <- MulticoreParam(workers = 25, progressbar = TRUE)

sce = readH5AD('/mnt/sdb/scz_meta_analysis_processed/anndata_objs/integrated_adata_full_annotation_dreamlet.h5ad',
	       use_hdf5=TRUE, raw=FALSE, verbose=FALSE, uns=FALSE)

sce$Donor_Sample <- paste(sce$Donor, sce$Sample.Name, sep="_")

pb <- aggregateToPseudoBulk(
  sce,
  assay = "X",
  sample_id = "Donor_Sample",
  cluster_id = "subtype_cluster_annotations_markers",
  BPPARAM = BPPARAM)

colData(pb)$Broad_Genotype <- make.names(colData(pb)$Broad_Genotype)

## Model Broad Genotype As Random Effect ##

# crumblr wont take cell level variables like n_counts or percent_mito #

formula_sep <- ~ (1|Broad_Genotype) + (1|Sample.Name) + (1|Donor) + (1|Chemistry) + (1|Manuscript) +(1|Sex) + Day + (1|Protocol) 

formula_sep_fixed <- ~ 0 + Broad_Genotype + (1|Sample.Name) + (1|Donor) + (1|Chemistry) + (1|Manuscript) + (1|Sex) + Day + (1|Protocol) 

cobj <- crumblr(cellCounts(pb))

# Variance Partitioning 

vp.c.sep <- fitExtractVarPartModel(cobj, formula_sep, colData(pb))
write.csv(as.data.frame(vp.c.sep), file = "/mnt/sdb/scz_meta_analysis_processed/dge_signatures/crumblr_outs/crumblr_broad_genotype_subtype.csv", row.names = TRUE)
write.csv(as.data.frame(cobj$E), file = "/mnt/sdb/scz_meta_analysis_processed/dge_signatures/crumblr_outs/crumblr_broad_genotype_matrix_subtype.csv", row.names = TRUE)

# Differential Testing #

L = makeContrastsDream(formula_sep_fixed, colData(pb),
		       contrasts = c(Broad_Genotype_22q11 = "Broad_GenotypeX22q112del - Broad_GenotypeControl", 
                                     Broad_Genotype_NRXN1 = "Broad_GenotypeNRXN1del - Broad_GenotypeControl", 
                                     Broad_Genotype_3q29 = "Broad_GenotypeX3q29del - Broad_GenotypeControl", 
                                     Broad_Genotype_15q13 = "Broad_GenotypeX15q133del - Broad_GenotypeControl", 
                                     Broad_Genotype_Idiopathic = "Broad_GenotypeIdiopathic_Schizophrenia - Broad_GenotypeControl"))

fit <- dream(cobj, formula_sep_fixed, colData(pb), L=L)

fit <- eBayes(fit)

# Top Tables #
dge_22q11 = topTable(fit, coef='Broad_Genotype_22q11', number=Inf)
dge_NRXN1 = topTable(fit, coef='Broad_Genotype_NRXN1', number=Inf)
dge_3q29 = topTable(fit, coef='Broad_Genotype_3q29', number=Inf)
dge_15q13 = topTable(fit, coef='Broad_Genotype_15q13', number=Inf)
dge_Idiopathic = topTable(fit, coef='Broad_Genotype_Idiopathic', number=Inf)

write.csv(dge_22q11, paste0(CSV_PATH, 'dream_crumblr_22q11_dge_scz_results_subtype.csv'), row.names = TRUE)
write.csv(dge_NRXN1, paste0(CSV_PATH, 'dream_crumblr_NRNXN1_dge_scz_results_subtype.csv'), row.names = TRUE)
write.csv(dge_3q29, paste0(CSV_PATH, 'dream_crumblr_3q29_dge_scz_results_subtype.csv'), row.names = TRUE)
write.csv(dge_15q13, paste0(CSV_PATH, 'dream_crumblr_15q13_dge_scz_results_subtype.csv'), row.names = TRUE)
write.csv(dge_Idiopathic, paste0(CSV_PATH, 'dream_crumblr_Idiopathic_dge_scz_results_subtype.csv'), row.names = TRUE)

