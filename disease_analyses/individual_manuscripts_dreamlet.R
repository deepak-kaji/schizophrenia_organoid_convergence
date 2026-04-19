library(zellkonverter)
library(BiocParallel)
library(SingleCellExperiment)
library(dreamlet)

BPPARAM <- MulticoreParam(workers = 25, progressbar = TRUE)

sce = readH5AD('/mnt/sdb/scz_meta_analysis_processed/anndata_objs/integrated_adata_full_annotation_dreamlet.h5ad', 
	       use_hdf5=TRUE, layers=FALSE, raw=FALSE, verbose=FALSE, uns=FALSE)

sce$Donor_Sample <- paste(sce$Donor, sce$Sample.Name, sep = "_") 

# Stressed Glia Too Small To Test 
sce <- sce[,sce$subclass_annotations_markers != 'Stressed Glia FTL+B2M+GLUL+']

pb <- aggregateToPseudoBulk(
  sce,
  assay = "X",
  sample_id = "Donor_Sample",
  cluster_id = "subclass_annotations_markers",
  BPPARAM = BPPARAM
)

colData(pb)$Broad_Genotype <- make.names(colData(pb)$Broad_Genotype)

## Step 2: Subset Out Manuscripts ##

fernando <- pb[,colData(pb)$Manuscript == 'Fernando'] 
#khan <- pb[,colData(pb)$Manuscript == 'Khan'] 
shin <- pb[,colData(pb)$Manuscript == 'Shin'] 
purcell <- pb[,colData(pb)$Manuscript == 'Purcell'] 
sebastian <- pb[,colData(pb)$Manuscript == 'Sebastian'] 
sawada <- pb[,colData(pb)$Manuscript == 'Sawada'] 
walsh <- pb[,colData(pb)$Manuscript == 'Walsh']
notaras<- pb[,colData(pb)$Manuscript == 'Notaras'] 
rao <- pb[,colData(pb)$Manuscript == 'Rao'] 

# walsh crashing on this cell type , not estimatable
walsh <- walsh[, names(assays(walsh)) != "Mesenchymal-like cells VIM+VCAN+SPARC+"]

# Step 3: Process assays again for DE model ---

# Trimmed model for differential expression

formula_trim <- ~ 0 + Broad_Genotype + (1|Donor) + (1|Sample.Name) + (1|Chemistry) + (1|Manuscript) +
            	(1|Protocol) + (1|Sex) + Day + scale(n_counts) + scale(percent_mito)  

#formula_trim <- ~ 0 + Broad_Genotype  

res.fernando <-processAssays(fernando, formula = formula_trim, min.cells = 5, min.count = 5, min.samples = 4, min.prop = 0.2, BPPARAM = BPPARAM)
#res.khan <-processAssays(khan, formula = formula_trim, min.cells = 5, min.count = 5, min.samples = 4, min.prop = 0.2, BPPARAM = BPPARAM)
res.shin <-processAssays(shin, formula = formula_trim, min.cells = 5, min.count = 5, min.samples = 4, min.prop = 0.2, BPPARAM = BPPARAM)
res.purcell <-processAssays(purcell, formula = formula_trim, min.cells = 5, min.count = 5, min.samples = 4, min.prop = 0.2, BPPARAM = BPPARAM)
res.sebastian <-processAssays(sebastian, formula = formula_trim, min.cells = 5, min.count = 5, min.samples = 4, min.prop = 0.2, BPPARAM = BPPARAM)
res.sawada <-processAssays(sawada, formula = formula_trim, min.cells = 5, min.count = 5, min.samples = 4, min.prop = 0.2, BPPARAM = BPPARAM)
res.walsh <-processAssays(walsh, formula = formula_trim, min.cells = 5, min.count = 5, min.samples = 4, min.prop = 0.2, BPPARAM = BPPARAM)
res.notaras <-processAssays(notaras, formula = formula_trim, min.cells = 5, min.count = 5, min.samples = 4, min.prop = 0.2, BPPARAM = BPPARAM)
res.rao <-processAssays(rao, formula = formula_trim, min.cells = 5, min.count = 5, min.samples = 4, min.prop = 0.2, BPPARAM = BPPARAM)

# Step 4: Run dreamlet DE analysis #Raodl.fernando <- dreamlet(res.fernando, formula = formula_trim, contrasts = c(Broad_Genotype_NRXN1 = "Broad_GenotypeNRXN1del - Broad_GenotypeControl"), BPPARAM = BPPARAM)
res.dl.fernando <- dreamlet(res.fernando, formula = formula_trim, contrasts = c(Broad_Genotype_NRXN1 = "Broad_GenotypeNRXN1del - Broad_GenotypeControl"), BPPARAM = BPPARAM)
#res.dl.khan <- dreamlet(res.khan, formula = formula_trim, contrasts = c(Broad_Genotype_22q11 = "Broad_GenotypeX22q112del - Broad_GenotypeControl"), BPPARAM = BPPARAM)
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
#dge_khan = topTable(res.dl.khan, coef='Broad_Genotype_22q11', number=Inf)
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
#dge_khan$dataset <- 'khan'
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

write.csv(dge_fernando, paste0(CSV_PATH, 'dreamlet_fernando_dge_scz_results.csv'), row.names = FALSE)
#write.csv(dge_khan, paste0(CSV_PATH, 'dreamlet_khan_dge_scz_results.csv'), row.names = FALSE)
write.csv(dge_shin, paste0(CSV_PATH, 'dreamlet_shin_dge_scz_results.csv'), row.names = FALSE)
write.csv(dge_purcell, paste0(CSV_PATH, 'dreamlet_purcell_dge_scz_results.csv'), row.names = FALSE)
write.csv(dge_sebastian, paste0(CSV_PATH, 'dreamlet_sebastian_dge_scz_results.csv'), row.names = FALSE)
write.csv(dge_sawada, paste0(CSV_PATH, 'dreamlet_sawada_dge_scz_results.csv'), row.names = FALSE)
write.csv(dge_notaras, paste0(CSV_PATH, 'dreamlet_notaras_dge_scz_results.csv'), row.names = FALSE)
write.csv(dge_rao, paste0(CSV_PATH, 'dreamlet_rao_dge_scz_results.csv'), row.names = FALSE)
write.csv(dge_walsh_15q13, paste0(CSV_PATH, 'dreamlet_walsh_15q13_dge_scz_results.csv'), row.names = FALSE)
write.csv(dge_walsh_22q11, paste0(CSV_PATH, 'dreamlet_walsh_22q11_dge_scz_results.csv'), row.names = FALSE)

write.csv(res.nrxn1, paste0(CSV_PATH, 'dreamlet_nrxn1_combined.csv'), row.names = FALSE)
write.csv(res.22q11, paste0(CSV_PATH, 'dreamlet_22q11_combined.csv'), row.names = FALSE)
write.csv(res.idiopathic, paste0(CSV_PATH, 'dreamlet_idiopathic_combined.csv'), row.names = FALSE)


