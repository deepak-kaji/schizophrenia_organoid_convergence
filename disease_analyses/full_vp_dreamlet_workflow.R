library(zellkonverter)
library(BiocParallel)
library(SingleCellExperiment)
library(dreamlet)

BPPARAM <- MulticoreParam(workers = 25, progressbar = TRUE)

sce = readH5AD('/mnt/sdb/scz_meta_analysis_processed/anndata_objs/FILLER.h5ad', 
	       use_hdf5=TRUE, layers=FALSE, raw=FALSE, verbose=FALSE, uns=FALSE)

#sce <- sce[rowData(sce)$robust_protein_coding,]

pb <- aggregateToPseudoBulk(
  sce,
  assay = "X",
  sample_id = "Run_Donor_Sample",
  cluster_id = "subclass_annotations",
  BPPARAM = BPPARAM
)

## consider Protocol? is Whitelist colinear with any of the other columns
formula_full <- ~ (1|Run) + (1|Donor) + (1|BioSample) + (1|Chemistry) + (1|Manuscript) + (1|Sex) + (1|Sex) + (1|SCZ_SUBTYPE) + (1|Whitelist) + scale(n_counts) + scale(percent_mito)  

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
ggsave(plot_voom_fig, file='/mnt/sdb/scz_meta_analysis_processed/dge_signatures/dreamlet_dges/disease_analyses/plot_voom.pdf')

# --- Step 3: Fit variance partition model ---
vp.lst <- fitVarPart(res.proc.vp, formula_full)

write.csv(vp.lst, file = "/mnt/sdb/scz_meta_analysis_processed/dge_signatures/dreamlet_dges/disease_analyses/variance_partition_long.csv", row.names = FALSE)

# Optional: Visualize variance explained
plot_varpart <- plotVarPart(vp.lst)
ggsave(plot_varpart, file='/mnt/sdb/scz_meta_analysis_processed/dge_signatures/dreamlet_dges/disease_analyses/visualize_variance_full.pdf')

# --- Step 4: Process assays again for DE model ---

# Trimmed model for differential expression

formula_trim <- ~ (1|Run) + (1|Donor) + (1|BioSample) + (1|Chemistry) + (1|Manuscript) + (1|Sex) + (1|Sex) + (1|Whitelist) + scale(n_counts) + scale(percent_mito) + (1|Diagnosis)   
	  
# as a factor, so you can treat as  fixed variable

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
  contrasts = c(Status_22q11 = "Status22q11 - StatusControl", 
		Status_NRXN1 = "StatusNRXN1 - StatusControl", 
		Status_3q29 = "Status3q29 - StatusControl", 
		Status_15q13 = "Status15q13 - StatusControl", 
		Status_Idiopathic = "StatusIdiopathic - StatusControl"), BPPARAM = BPPARAM)

# Combine results across all assays

22q11_dge = topTable(res.dl, coef='Status_22q11', number=Inf)
NRXN1_dge = topTable(res.dl, coef='Status_NRXN1', number=Inf)
3q29_dge = topTable(res.dl, coef='Status_3q29', number=Inf)
15q13_dge = topTable(res.dl, coef='Status_15q13', number=Inf)
Idiopathic_dge = topTable(res.dl, coef='Status_Idiopathic', number=Inf)

# Save as CSV

CSV_PATH <- '/mnt/sdb/scz_meta_analysis_processed/dge_signatures/dreamlet_dges/disease_analyses/'

write.csv(22q11_dge, paste0(CSV_PATH, 'dreamlet_22q11_dge_scz_results.csv', row.names = FALSE))
write.csv(NRXN1_dge, paste0(CSV_PATH, 'dreamlet_NRNXN1_dge_scz_results.csv', row.names = FALSE))
write.csv(3q29_dge, paste0(CSV_PATH, 'dreamlet_3q29_dge_scz_results.csv', row.names = FALSE))
write.csv(15q13_dge, paste0(CSV_PATH, 'dreamlet_15q13_dge_scz_results.csv', row.names = FALSE))
write.csv(Idiopathic_dge, paste0(CSV_PATH, 'dreamlet_Idiopathic_dge_scz_results.csv', row.names = FALSE))
