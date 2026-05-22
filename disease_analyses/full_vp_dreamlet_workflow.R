library(zellkonverter)
library(BiocParallel)
library(SingleCellExperiment)
library(dreamlet)

BPPARAM <- MulticoreParam(workers = 25, progressbar = TRUE)

sce = readH5AD('/mnt/sdb/scz_meta_analysis_processed/anndata_objs/integrated_adata_full_annotation_dreamlet.h5ad', 
	       use_hdf5=TRUE, layers=FALSE, raw=FALSE, verbose=FALSE, uns=FALSE)

sce$Donor_Sample <- paste(sce$Donor, sce$Sample.Name, sep = "_") 

pb <- aggregateToPseudoBulk(
  sce,
  assay = "X",
  sample_id = "Donor_Sample",
  cluster_id = "CellType",
  BPPARAM = BPPARAM
)

colData(pb)$Broad_Genotype <- make.names(colData(pb)$Broad_Genotype)

formula_full <- ~ (1|Broad_Genotype) + (1|Donor) + (1|Sample.Name) + (1|Chemistry) + (1|Manuscript) + Day + (1|Sex) + (1|Protocol) + scale(n_counts) + scale(percent_mito)  


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
ggsave(plot_voom_fig, file='/mnt/sdb/scz_meta_analysis_processed/dge_signatures/dreamlet_dges/disease_analyses/plot_voom.png', dpi=500)

# --- Step 3: Fit variance partition model ---

vp.lst <- fitVarPart(res.proc.vp, formula_full)

write.csv(vp.lst, file = "/mnt/sdb/scz_meta_analysis_processed/dge_signatures/dreamlet_dges/disease_analyses/variance_partition_long.csv", row.names = FALSE)

# Optional: Visualize variance explained
plot_varpart <- plotVarPart(vp.lst)
ggsave(plot_varpart, file='/mnt/sdb/scz_meta_analysis_processed/dge_signatures/dreamlet_dges/disease_analyses/visualize_variance_full.png', dpi=500)

# --- Step 4: Process assays again for DE model ---

# Trimmed model for differential expression

# removed Protocol 

formula_trim <- ~ 0 + Broad_Genotype + (1|Donor) + (1|Sample.Name) + (1|Chemistry) + (1|Manuscript) + Day +
            	 (1|Sex) + (1|Protocol) + scale(n_counts) + scale(percent_mito)  

res.proc.de <- processAssays(
  pb,
  formula = formula_trim,
  min.cells = 5,
  min.count = 5,
  min.samples = 4,
  min.prop = 0.2,
  BPPARAM = BPPARAM
)

# --- Step 5: Run dreamlet DE analysis ---

res.dl <- dreamlet(
  res.proc.de,
  formula = formula_trim,
  contrasts = c(Broad_Genotype_22q11 = "Broad_GenotypeX22q112del - Broad_GenotypeControl", 
		Broad_Genotype_NRXN1 = "Broad_GenotypeNRXN1del - Broad_GenotypeControl", 
		Broad_Genotype_3q29 = "Broad_GenotypeX3q29del - Broad_GenotypeControl", 
		Broad_Genotype_15q13 = "Broad_GenotypeX15q133del - Broad_GenotypeControl", 
		Broad_Genotype_Idiopathic = "Broad_GenotypeIdiopathic_Schizophrenia - Broad_GenotypeControl"), BPPARAM = BPPARAM)

# Combine results across all assays

dge_22q11 = topTable(res.dl, coef='Broad_Genotype_22q11', number=Inf)
dge_NRXN1 = topTable(res.dl, coef='Broad_Genotype_NRXN1', number=Inf)
dge_3q29 = topTable(res.dl, coef='Broad_Genotype_3q29', number=Inf)
dge_15q13 = topTable(res.dl, coef='Broad_Genotype_15q13', number=Inf)
dge_Idiopathic = topTable(res.dl, coef='Broad_Genotype_Idiopathic', number=Inf)

# Save as CSV

CSV_PATH <- '/mnt/sdb/scz_meta_analysis_processed/dge_signatures/dreamlet_dges/disease_analyses/'

write.csv(dge_22q11, paste0(CSV_PATH, 'dreamlet_22q11_dge_scz_results.csv'), row.names = FALSE)
write.csv(dge_NRXN1, paste0(CSV_PATH, 'dreamlet_NRNXN1_dge_scz_results.csv'), row.names = FALSE)
write.csv(dge_3q29, paste0(CSV_PATH, 'dreamlet_3q29_dge_scz_results.csv'), row.names = FALSE)
write.csv(dge_15q13, paste0(CSV_PATH, 'dreamlet_15q13_dge_scz_results.csv'), row.names = FALSE)
write.csv(dge_Idiopathic, paste0(CSV_PATH, 'dreamlet_Idiopathic_dge_scz_results.csv'), row.names = FALSE)
