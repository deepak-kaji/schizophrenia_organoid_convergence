library(zellkonverter)
library(BiocParallel)
library(SingleCellExperiment)
library(dreamlet)

BPPARAM <- MulticoreParam(workers = 25, progressbar = TRUE)

sce = readH5AD('/mnt/sdb/scz_meta_analysis_processed/anndata_objs/integrated_adata_full_annotation_scanvi_loaded_dreamlet.h5ad', 
	       use_hdf5=TRUE, layers=FALSE, raw=FALSE, verbose=FALSE, uns=FALSE)

pb <- aggregateToPseudoBulk(
  sce,
  assay = "X",
  sample_id = "Run_Donor_Sample",
  cluster_id = "subclass_annotations",
  BPPARAM = BPPARAM
)

colData(pb)$Broad_Genotype <- make.names(colData(pb)$Broad_Genotype)

## consider Protocol? is Whitelist colinear with any of the other columns
formula_full <- ~ (1|Run) + (1|Donor) + (1|Sample.Name) + (1|Chemistry) + (1|Manuscript) + (1|Sex) + 
	          (1|Broad_Genotype) + (1|Protocol) + scale(n_counts) + scale(percent_mito)  

res.proc.vp <- processAssays(
  pb,
  formula = formula_full,
  min.cells = 5,
  min.count = 5,
  min.samples = 4,
  min.prop = 0.2,
  BPPARAM = BPPARAM
)

plot_voom_fig = plotVoom(res.proc.vp, ncol=4)
ggsave(plot_voom_fig, file='/mnt/sdb/scz_meta_analysis_processed/dge_signatures/dreamlet_dges/disease_analyses/scvi_space_plot_voom.pdf')

# --- Step 3: Fit variance partition model ---

vp.lst <- fitVarPart(res.proc.vp, formula_full)

write.csv(vp.lst, file = "/mnt/sdb/scz_meta_analysis_processed/dge_signatures/dreamlet_dges/disease_analyses/scvi_space_variance_partition_long.csv", row.names = FALSE)

# Optional: Visualize variance explained
plot_varpart <- plotVarPart(vp.lst)
ggsave(plot_varpart, file='/mnt/sdb/scz_meta_analysis_processed/dge_signatures/dreamlet_dges/disease_analyses/scvi_space_visualize_variance_full.pdf')

