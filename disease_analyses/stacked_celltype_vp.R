library(dplyr)
library(zellkonverter)
library(BiocParallel)
library(SingleCellExperiment)
library(dreamlet)

BPPARAM <- MulticoreParam(workers = 25, progressbar = TRUE)

sce = readH5AD('/mnt/sdb/scz_meta_analysis_processed/anndata_objs/integrated_adata_full_annotation_dreamlet.h5ad',
	       use_hdf5=TRUE, layers=TRUE, raw=FALSE, verbose=FALSE, uns=FALSE)

sce$Donor_Sample <- paste(sce$Donor, sce$Sample.Name, sep="_")

pb <- aggregateToPseudoBulk(
  sce,
  assay = "X",
  sample_id = "Donor_Sample",
  cluster_id = "subclass_annotations_markers",
  BPPARAM = BPPARAM
)

pb_stacked <- stackAssays(pb)

formula_full <- ~ (1|Broad_Genotype) + (1|Donor) + (1|Sample.Name) + (1|Chemistry) + (1|Manuscript) + (1|Sex) + Day + (1|Protocol) + scale(n_counts) + scale(percent_mito) + (1|stackedAssay)

res.proc.vp <- processAssays(
  pb_stacked,
  formula = formula_full,
  min.cells = 5,
  min.count = 5,
  min.samples = 5,
  min.prop = 0.2,
  BPPARAM = BPPARAM
)

plot_voom_fig = plotVoom(res.proc.vp, ncol=4)
ggsave(plot_voom_fig, file='/mnt/sdb/scz_meta_analysis_processed/dge_signatures/dreamlet_dges/disease_analyses/plot_voom_stacked.png', dpi=500)

# --- Step 3: Fit variance partition model ---
vp.lst <- fitVarPart(res.proc.vp, formula_full)
write.csv(vp.lst, file = "/mnt/sdb/scz_meta_analysis_processed/dge_signatures/dreamlet_dges/disease_analyses/variance_partition_long_stacked.csv", row.names = FALSE)

# Optional: Visualize variance explained
plot_varpart <- plotVarPart(vp.lst)
ggsave(plot_varpart, file='/mnt/sdb/scz_meta_analysis_processed/dge_signatures/dreamlet_dges/disease_analyses/visualize_variance_stacked.png', dpi=500)

