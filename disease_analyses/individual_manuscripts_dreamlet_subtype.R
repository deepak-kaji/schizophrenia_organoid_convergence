library(zellkonverter)
library(BiocParallel)
library(SingleCellExperiment)
library(dreamlet)

BPPARAM <- MulticoreParam(workers = 25, progressbar = TRUE)

sce = readH5AD('/mnt/sdb/scz_meta_analysis_processed/anndata_objs/integrated_adata_full_annotation_dreamlet.h5ad', 
	       use_hdf5=TRUE, layers=FALSE, raw=FALSE, verbose=FALSE, uns=FALSE)

sce$Donor_Sample <- paste(sce$Donor, sce$Sample.Name, sep = "_") 

# Stressed Glia Too Small To Test 
sce <- sce[,sce$CellType != 'Stressed Glia FTL+B2M+GLUL+']

colData(sce)$Broad_Genotype <- make.names(colData(sce)$Broad_Genotype)

## Step 2: Subset Out Manuscripts ##

# khan not estimatable due to small donor study (2x2 case control) resulted in unstable variance downstream
fernando <- sce[,colData(sce)$Manuscript == 'Fernando'] 
shin <- sce[,colData(sce)$Manuscript == 'Shin'] 
purcell <- sce[,colData(sce)$Manuscript == 'Purcell'] 
sebastian <- sce[,colData(sce)$Manuscript == 'Sebastian'] 
sawada <- sce[,colData(sce)$Manuscript == 'Sawada'] 
walsh <- sce[,colData(sce)$Manuscript == 'Walsh']
notaras<- sce[,colData(sce)$Manuscript == 'Notaras'] 
rao <- sce[,colData(sce)$Manuscript == 'Rao'] 

## purcell doesnt contribute meaningfully to OPC/Oligo --> drop
purcell_keep <- !grepl("OPC|Oligodendrocyte|Astro", purcell$CellType)
purcell <- purcell[, purcell_keep]

walsh_keep <- !grepl("OPC|Oligodendrocyte", walsh$CellType)
walsh <- walsh[, walsh_keep]

pb_fernando <- aggregateToPseudoBulk(fernando, assay = "X", sample_id = "Donor_Sample", cluster_id = "subtype_cluster_annotations_markers", BPPARAM = BPPARAM)
pb_shin <- aggregateToPseudoBulk(shin, assay = "X", sample_id = "Donor_Sample", cluster_id = "subtype_cluster_annotations_markers", BPPARAM = BPPARAM)
pb_purcell <- aggregateToPseudoBulk(purcell, assay = "X", sample_id = "Donor_Sample", cluster_id = "subtype_cluster_annotations_markers", BPPARAM = BPPARAM)
pb_sebastian <- aggregateToPseudoBulk(sebastian, assay = "X", sample_id = "Donor_Sample", cluster_id = "subtype_cluster_annotations_markers", BPPARAM = BPPARAM)
pb_sawada <- aggregateToPseudoBulk(sawada, assay = "X", sample_id = "Donor_Sample", cluster_id = "subtype_cluster_annotations_markers", BPPARAM = BPPARAM)
pb_walsh <- aggregateToPseudoBulk(walsh, assay = "X", sample_id = "Donor_Sample", cluster_id = "subtype_cluster_annotations_markers", BPPARAM = BPPARAM)
pb_notaras <- aggregateToPseudoBulk(notaras, assay = "X", sample_id = "Donor_Sample", cluster_id = "subtype_cluster_annotations_markers", BPPARAM = BPPARAM)
pb_rao <- aggregateToPseudoBulk(rao, assay = "X", sample_id = "Donor_Sample", cluster_id = "subtype_cluster_annotations_markers", BPPARAM = BPPARAM)

# walsh crashing on this cell type , not estimatable

assays(pb_walsh) <- assays(pb_walsh)[names(assays(pb_walsh)) != "Mesenchymal-like cells VIM+VCAN+SPARC+"]

# Step 3: Process assays again for DE model ---

# Trimmed model for differential expression

# Note in meta-analyses, Protocol is kept as a covariate to ensure proper modeling for Fernando which uses 2 protocols within the same study
# In all other studies, protocol will be dropped since it exhibits no variance 

formula_trim <- ~ 0 + Broad_Genotype + (1|Donor) + (1|Sample.Name) + (1|Chemistry) + (1|Manuscript) +
            	(1|Protocol) + (1|Sex) + Day + scale(n_counts) + scale(percent_mito)  

res.fernando <-processAssays(pb_fernando, formula = formula_trim, min.cells = 5, min.count = 5, min.samples = 4, min.prop = 0.2, BPPARAM = BPPARAM)
res.shin <-processAssays(pb_shin, formula = formula_trim, min.cells = 5, min.count = 5, min.samples = 4, min.prop = 0.2, BPPARAM = BPPARAM)
res.purcell <-processAssays(pb_purcell, formula = formula_trim, min.cells = 5, min.count = 5, min.samples = 4, min.prop = 0.2, BPPARAM = BPPARAM)
res.sebastian <-processAssays(pb_sebastian, formula = formula_trim, min.cells = 5, min.count = 5, min.samples = 4, min.prop = 0.2, BPPARAM = BPPARAM)
res.sawada <-processAssays(pb_sawada, formula = formula_trim, min.cells = 5, min.count = 5, min.samples = 4, min.prop = 0.2, BPPARAM = BPPARAM)
res.walsh <-processAssays(pb_walsh, formula = formula_trim, min.cells = 5, min.count = 5, min.samples = 4, min.prop = 0.2, BPPARAM = BPPARAM)
res.notaras <-processAssays(pb_notaras, formula = formula_trim, min.cells = 5, min.count = 5, min.samples = 4, min.prop = 0.2, BPPARAM = BPPARAM)
res.rao <-processAssays(pb_rao, formula = formula_trim, min.cells = 5, min.count = 5, min.samples = 4, min.prop = 0.2, BPPARAM = BPPARAM)

# Step 4: Run dreamlet DE analysis #Raodl.fernando <- dreamlet(res.fernando, formula = formula_trim, contrasts = c(Broad_Genotype_NRXN1 = "Broad_GenotypeNRXN1del - Broad_GenotypeControl"), BPPARAM = BPPARAM)
res.dl.fernando <- dreamlet(res.fernando, formula = formula_trim, contrasts = c(Broad_Genotype_NRXN1 = "Broad_GenotypeNRXN1del - Broad_GenotypeControl"), BPPARAM = BPPARAM)
res.dl.shin <- dreamlet(res.shin, formula = formula_trim, contrasts = c(Broad_Genotype_22q11 = "Broad_GenotypeX22q112del - Broad_GenotypeControl"), BPPARAM = BPPARAM)
res.dl.purcell <- dreamlet(res.purcell, formula = formula_trim, contrasts = c(Broad_Genotype_3q29 = "Broad_GenotypeX3q29del - Broad_GenotypeControl"), BPPARAM = BPPARAM)
res.dl.sebastian <- dreamlet(res.sebastian, formula = formula_trim, contrasts = c(Broad_Genotype_NRXN1 = "Broad_GenotypeNRXN1del - Broad_GenotypeControl"), BPPARAM = BPPARAM)
res.dl.sawada <- dreamlet(res.sawada, formula = formula_trim,
			  contrasts = c(Broad_Genotype_Idiopathic = "Broad_GenotypeIdiopathic_Schizophrenia - Broad_GenotypeControl"), BPPARAM = BPPARAM)
res.dl.notaras <- dreamlet(res.notaras, formula = formula_trim,
			   contrasts = c(Broad_Genotype_Idiopathic = "Broad_GenotypeIdiopathic_Schizophrenia - Broad_GenotypeControl"), BPPARAM = BPPARAM)
res.dl.rao <- dreamlet(res.rao, formula = formula_trim, contrasts = c(Broad_Genotype_22q11 = "Broad_GenotypeX22q112del - Broad_GenotypeControl"), BPPARAM = BPPARAM)
res.dl.walsh <- dreamlet(res.walsh, formula = formula_trim, 
			 contrasts = c(Broad_Genotype_15q13 = "Broad_GenotypeX15q133del - Broad_GenotypeControl",
				       Broad_Genotype_22q11 = "Broad_GenotypeX22q112del - Broad_GenotypeControl"), BPPARAM = BPPARAM)

# Combine results across all assays

dge_fernando = topTable(res.dl.fernando, coef='Broad_Genotype_NRXN1', number=Inf)
dge_shin = topTable(res.dl.shin, coef='Broad_Genotype_22q11', number=Inf)
dge_purcell = topTable(res.dl.purcell, coef='Broad_Genotype_3q29', number=Inf)
dge_sebastian = topTable(res.dl.sebastian, coef='Broad_Genotype_NRXN1', number=Inf)
dge_sawada = topTable(res.dl.sawada, coef='Broad_Genotype_Idiopathic', number=Inf)
dge_notaras = topTable(res.dl.notaras, coef='Broad_Genotype_Idiopathic', number=Inf)
dge_rao = topTable(res.dl.rao, coef='Broad_Genotype_22q11', number=Inf)
dge_walsh_15q13 = topTable(res.dl.walsh, coef='Broad_Genotype_15q13', number=Inf)
dge_walsh_22q11 = topTable(res.dl.walsh, coef='Broad_Genotype_22q11', number=Inf)

# rbind them for meta-analysis #

dge_fernando$dataset <- 'fernando'
dge_shin$dataset <- 'shin'
dge_purcell$dataset <- 'purcell'
dge_sebastian$dataset <- 'sebastian'
dge_sawada$dataset <- 'sawada'
dge_notaras$dataset <- 'notaras'
dge_rao$dataset <- 'rao'
dge_walsh_15q13$dataset <- 'walsh_15q13'
dge_walsh_22q11$dataset <- 'walsh_22q11'

# Combine all datasets

dataset_nrnx1 <- rbind(dge_fernando, dge_sebastian)
dataset_22q11 <- rbind(dge_shin, dge_rao, dge_walsh_22q11)
dataset_idiopathic <- rbind(dge_sawada, dge_notaras)

# run meta-analysis #

res.nrxn1 <- meta_analysis(dataset_nrnx1, method = 'FE')
res.22q11 <- meta_analysis(dataset_22q11, method = 'FE')
res.idiopathic <- meta_analysis(dataset_idiopathic, method = 'FE')

# Save as CSV

CSV_PATH <- '/mnt/sdb/scz_meta_analysis_processed/dge_signatures/dreamlet_dges/disease_analyses/'

write.csv(dge_fernando, paste0(CSV_PATH, 'dreamlet_fernando_dge_scz_results_subtype.csv'), row.names = FALSE)
#write.csv(dge_khan, paste0(CSV_PATH, 'dreamlet_khan_dge_scz_results_subtype.csv'), row.names = FALSE)
write.csv(dge_shin, paste0(CSV_PATH, 'dreamlet_shin_dge_scz_results_subtype.csv'), row.names = FALSE)
write.csv(dge_purcell, paste0(CSV_PATH, 'dreamlet_purcell_dge_scz_results_subtype.csv'), row.names = FALSE)
write.csv(dge_sebastian, paste0(CSV_PATH, 'dreamlet_sebastian_dge_scz_results_subtype.csv'), row.names = FALSE)
write.csv(dge_sawada, paste0(CSV_PATH, 'dreamlet_sawada_dge_scz_results_subtype.csv'), row.names = FALSE)
write.csv(dge_notaras, paste0(CSV_PATH, 'dreamlet_notaras_dge_scz_results_subtype.csv'), row.names = FALSE)
write.csv(dge_rao, paste0(CSV_PATH, 'dreamlet_rao_dge_scz_results_subtype.csv'), row.names = FALSE)
write.csv(dge_walsh_15q13, paste0(CSV_PATH, 'dreamlet_walsh_15q13_dge_scz_results_subtype.csv'), row.names = FALSE)
write.csv(dge_walsh_22q11, paste0(CSV_PATH, 'dreamlet_walsh_22q11_dge_scz_results_subtype.csv'), row.names = FALSE)

write.csv(res.nrxn1, paste0(CSV_PATH, 'dreamlet_nrxn1_combined_subtype.csv'), row.names = FALSE)
write.csv(res.22q11, paste0(CSV_PATH, 'dreamlet_22q11_combined_subtype.csv'), row.names = FALSE)
write.csv(res.idiopathic, paste0(CSV_PATH, 'dreamlet_idiopathic_combined_subtype.csv'), row.names = FALSE)


