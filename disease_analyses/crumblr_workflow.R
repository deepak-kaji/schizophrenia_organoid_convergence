library(dplyr)
library(zellkonverter)
library(BiocParallel)
library(SingleCellExperiment)
library(dreamlet)
library(crumblr)

BPPARAM <- MulticoreParam(workers = 25, progressbar = TRUE)

sce = readH5AD('/mnt/sdb/scz_meta_analysis_processed/anndata_objs/integrated_adata_full_annotation_dreamlet.h5ad',
	       use_hdf5=TRUE, raw=FALSE, verbose=FALSE, uns=FALSE)

sce$Donor_Sample <- paste(sce$Donor, sce$Sample.Name, sep="_")

pb <- aggregateToPseudoBulk(
  sce,
  assay = "X",
  sample_id = "Donor_Sample",
  cluster_id = "CellType",
#  cluster_id = "subclass_annotations_markers",
  BPPARAM = BPPARAM)

colData(pb)$n_counts <- metadata(pb)$aggr_means$n_counts[
  match(as.character(rownames(colData(pb))),
        as.character(metadata(pb)$aggr_means$Donor_Sample))
]

colData(pb)$percent_mito <- metadata(pb)$aggr_means$percent_mito[
  match(as.character(rownames(colData(pb))),
        as.character(metadata(pb)$aggr_means$Donor_Sample))
]

## Model Broad Genotype As Random Effect ##

#formula_sep <- ~ (1|Broad_Genotype) + (1|Donor) + (1|Sample.Name) + (1|Chemistry) + (1|Manuscript) +(1|Sex) + Day + (1|Protocol) + scale(n_counts) + scale(percent_mito)
formula_sep <- ~ (1|Broad_Genotype) + (1|Donor) + (1|Chemistry) + (1|Manuscript) +(1|Sex) + Day + (1|Protocol) + scale(n_counts) + scale(percent_mito)

cobj <- crumblr(cellCounts(pb))
vp.c.sep <- fitExtractVarPartModel(cobj, formula_sep, colData(pb))

write.csv(as.data.frame(vp.c.sep), file = "/mnt/sdb/scz_meta_analysis_processed/dge_signatures/crumblr_outs/crumblr_broad_genotype.csv", row.names = TRUE)

## Model Broad Genotype But As Fixed Effect #

#formula_sep_fixed <- ~ Broad_Genotype + (1|Donor) + (1|Sample.Name) + (1|Chemistry) + (1|Manuscript) + (1|Sex) + Day + (1|Protocol) + scale(n_counts) + scale(percent_mito)
formula_sep_fixed <- ~ Broad_Genotype + (1|Donor) + (1|Chemistry) + (1|Manuscript) + (1|Sex) + Day + (1|Protocol) + scale(n_counts) + scale(percent_mito)

cobj_fixed <- crumblr(cellCounts(pb))
vp.c.sep.fixed <- fitExtractVarPartModel(cobj_fixed, formula_sep, colData(pb))

write.csv(as.data.frame(vp.c.sep.fixed), file = "/mnt/sdb/scz_meta_analysis_processed/dge_signatures/crumblr_outs/crumblr_broad_genotype_fixed.csv", row.names = TRUE)

## Model effect of Psychosis Overall##

#formula_full <- ~ Psychosis + (1|Donor) + (1|Sample.Name) + (1|Chemistry) + (1|Manuscript) + (1|Sex) + Day + (1|Protocol) + scale(n_counts) + scale(percent_mito)
formula_full <- ~ Psychosis + (1|Donor) + (1|Chemistry) + (1|Manuscript) + (1|Sex) + Day + (1|Protocol) + scale(n_counts) + scale(percent_mito)

cobj <- crumblr(cellCounts(pb))
vp.c <- fitExtractVarPartModel(cobj, formula_full, colData(pb))

write.csv(as.data.frame(vp.c), file = "/mnt/sdb/scz_meta_analysis_processed/dge_signatures/crumblr_outs/crumblr_psychosis.csv", row.names = TRUE)
