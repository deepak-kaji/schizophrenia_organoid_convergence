library(zellkonverter)
library(BiocParallel)
library(SingleCellExperiment)
library(dreamlet)
library(variancePartition)

BPPARAM <- MulticoreParam(workers = 25, progressbar = TRUE)

sce = readH5AD('/mnt/sdb/scz_meta_analysis_processed/anndata_objs/integrated_adata_full_annotation_scanvi_loaded_dreamlet.h5ad', 
	       use_hdf5=TRUE, layers=FALSE, raw=FALSE, verbose=FALSE, uns=FALSE)

Z <- assay(sce, "X")
expr <- as.matrix(Z)

meta <- as.data.frame(colData(sce))

meta$Donor <- factor(meta$Donor)
meta$Sample.Name <- factor(meta$Sample.Name)
meta$Broad_Genotype <- factor(meta$Broad_Genotype)
meta$Chemistry <- factor(meta$Chemistry)
meta$Manuscript <- factor(meta$Manuscript)
meta$Sex <- factor(meta$Sex)
meta$Protocol <- factor(meta$Protocol)
meta$subclass_annotations_markers <- factor(meta$subclass_annotations_markers)

## consider Protocol? is Whitelist colinear with any of the other columns
formula_full <- ~ (1|Broad_Genotype) + (1|Donor) + (1|Sample.Name) + (1|Chemistry) + (1|Manuscript) + Day +
                  (1|Sex) + (1|Protocol) + (1|subclass_annotations_markers) + scale(n_counts) + scale(percent_mito)  


vp <- fitExtractVarPartModel(expr, formula_full, meta, BPPARAM=BPPARAM)

write.csv(vp, file = "/mnt/sdb/scz_meta_analysis_processed/dge_signatures/dreamlet_dges/disease_analyses/scvi_space_variance_partition_long.csv", row.names = FALSE)

# Optional: Visualize variance explained
plot_varpart <- plotVarPart(vp)
ggsave(plot_varpart, file='/mnt/sdb/scz_meta_analysis_processed/dge_signatures/dreamlet_dges/disease_analyses/scvi_space_visualize_variance_full.png', dpi=500)

