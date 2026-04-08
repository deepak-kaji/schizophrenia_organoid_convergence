library(dplyr)
library(zellkonverter)
library(BiocParallel)
library(SingleCellExperiment)
library(dreamlet)
library(crumblr)

BPPARAM <- MulticoreParam(workers = 25, progressbar = TRUE)

sce = readH5AD('/mnt/sdb/scz_meta_analysis/anndata_objs/integrated_adata_full_annotation_dreamlet.h5add',
	       use_hdf5=TRUE, raw=FALSE, verbose=FALSE, uns=FALSE)

pb <- aggregateToPseudoBulk(
  sce,
  assay = "X",
  sample_id = "Run_Donor_Sample",
  cluster_id = "subclass_annotations",
  BPPARAM = BPPARAM)

formula_sep <- ~ (1|Run) + (1|Donor) + (1|Sample.Name) + (1|Chemistry) + (1|Manuscript) + (1|Sex) + 
	          (1|Broad_Genotype) + (1|Whitelist) + (1|Protocol) + scale(n_counts) + scale(percent_mito)

cobj <- crumblr(cellCounts(pb))
vp.c.sep <- fitExtractVarPartModel(cobj, formula_sep, colData(pb))

write.csv(as.data.frame(vp.c.sep), file = "/mnt/sdb/scz_meta_analysis/dge_signatures/crumblr_outs/crumblr_vp_neo_organoids.csv", row.names = TRUE)



formula_full <- ~ (1|Run) + (1|Donor) + (1|Sample.Name) + (1|Chemistry) + (1|Manuscript) + (1|Sex) + 
	          (1|Psychosis) + (1|Whitelist) + (1|Protocol) + scale(n_counts) + scale(percent_mito)

cobj <- crumblr(cellCounts(pb))
vp.c <- fitExtractVarPartModel(cobj, formula_full, colData(pb))

write.csv(as.data.frame(vp.c), file = "/mnt/sdb/scz_meta_analysis/dge_signatures/crumblr_outs/crumblr_vp_neo_organoids.csv", row.names = TRUE)
