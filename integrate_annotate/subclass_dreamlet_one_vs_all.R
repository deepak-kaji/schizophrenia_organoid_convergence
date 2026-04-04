library(zellkonverter)
library(BiocParallel)
library(SingleCellExperiment)
library(dreamlet)
library(dplyr)

BPPARAM <- MulticoreParam(workers = 25, progressbar = TRUE)

output_dir <- '/mnt/sdb/scz_meta_analysis_processed/dge_signatures/dreamlet_dges/one_vs_all/subclass'

sce = readH5AD('/mnt/sdb/scz_meta_analysis_processed/anndata_objs/one_versus_all_for_dreamlet.h5ad',
	       use_hdf5=TRUE, layers=FALSE, raw=FALSE, verbose=FALSE, uns=FALSE)

pbObj <- aggregateToPseudoBulk(
  sce,
  assay = "X",
  sample_id = "Run_Donor_Sample",
  cluster_id = "subclass",
  BPPARAM = BPPARAM)

# Get Assay Names  ---

mg = assayNames(pbObj)

# Evaluate the specificity of each gene for each cluster

df_cts = cellTypeSpecificity(pbObj)
write.csv(df_cts, file.path(output_dir, "cellTypeSpecificity.csv"), row.names = TRUE)

# Make a list of fits

fitList = list()

for(i in c(1:length(mg))){
    
    print(mg[i])
    sce_limited <- sce[, sce$class == sub("_.*", "", mg[i])]

    # run comparison of 1 cluster vs rest
    mg_one_minus <- unique(sce_limited$subclass)
    mg_one_minus <- mg_one_minus[mg_one_minus != mg[i]]
    fit = dreamletCompareClusters(pbObj, list(test = mg[i], baseline = mg_one_minus), method='none', min.cells=5, min.count=1, errorsAsWarnings = TRUE) # method='fixed'

    # store the fit obj
    fitList[[mg[i]]] = fit
    
    # up- or down-reg genes
    res = topTable(fit, coef='compare', number=Inf, sort.by='logFC', lfc=0)
    res$subclass = mg[i]

    # save table as CSF for python import
    out_file <- file.path(output_dir, paste0("DE_", mg[i], ".csv"))
    write.csv(res, out_file, row.names=TRUE)

    # specificity
    options(repr.plot.width=5, repr.plot.height=5)
}

# create a dreamletResult form this list
res.compare = new("dreamletResult", fitList)
saveRDS(res.compare, file=file.path(output_dir, 'one_vs_all_fit_dreamletResult.rds'))
