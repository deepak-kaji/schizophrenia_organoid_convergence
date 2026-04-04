library(zellkonverter)
library(BiocParallel)
library(SingleCellExperiment)
library(dreamlet)
library(dplyr)

BPPARAM <- MulticoreParam(workers = 25, progressbar = TRUE)

output_dir <- '/mnt/sdb/scz_meta_analysis_processed/dge_signatures/dreamlet_dges/one_vs_all/subtypes_remerged'

sce = readH5AD('/mnt/sdb/scz_meta_analysis_processed/anndata_objs/one_versus_all_for_dreamlet_remerged.h5ad',
	       use_hdf5=TRUE, layers=FALSE, raw=FALSE, verbose=FALSE, uns=FALSE)

pbObj <- aggregateToPseudoBulk(
  sce,
  assay = "X",
  sample_id = "Run_Donor_Sample",
  cluster_id = "subtypes_remerged",
  BPPARAM = BPPARAM)

# Get Assay Names  ---

mg = assayNames(pbObj)

# Evaluate the specificity of each gene for each cluster

df_cts = cellTypeSpecificity(pbObj)
write.csv(df_cts, file.path(output_dir, "cellTypeSpecificity.csv"), row.names = TRUE)

# Make a list of fits

fitList = list()

drop_log <- data.frame(
  subtype = character(),
  dropped = integer(),
  total = integer(),
  stringsAsFactors = FALSE
)

for(i in c(1:length(mg))){
    
    print(mg[i])
    sce_limited <- sce[, sce$subclass == sub("^(([^_]+_[^_]+)).*$", "\\1", mg[i])]

    # run comparison of 1 cluster vs rest
    mg_one_minus <- unique(sce_limited$subtypes_remerged)
    mg_one_minus <- mg_one_minus[mg_one_minus != mg[i]]
    
    if(length(mg_one_minus) == 0){
    message(paste("Skipping", mg[i], "- no baseline clusters"))

    drop_log <- rbind(drop_log, data.frame(
        subtype = mg[i],
        dropped = NA,
        total = NA
    ))

    next  # <-- skips dreamlet AND everything below
    }

    # Capture output
    captured <- capture.output({
        fit = dreamletCompareClusters(
            pbObj,
            list(test = mg[i], baseline = mg_one_minus),
            method='none',
            min.cells=5,
            min.count=1,
            errorsAsWarnings = TRUE
        )
    }, type='message')

    # Extract "Dropped X/Y samples"
    drop_line <- grep("Dropped", captured, value = TRUE)

    if(length(drop_line) > 0){
        nums <- regmatches(drop_line, gregexpr("[0-9]+", drop_line))[[1]]
        dropped <- as.integer(nums[1])
        total <- as.integer(nums[2])
    } else {
        dropped <- NA
        total <- NA
    }

    print(dropped)
    print(total)
    
    # Store log
    drop_log <- rbind(drop_log, data.frame(
        subtype = mg[i],
        dropped = dropped,
        total = total
    ))
    # store the fit obj
    fitList[[mg[i]]] = fit
    
    # up- or down-reg genes
    res = topTable(fit, coef='compare', number=Inf, sort.by='logFC', lfc=0)
    res$subtype = mg[i]

    # save table as CSF for python import
    out_file <- file.path(output_dir, paste0("DE_", mg[i], ".csv"))
    write.csv(res, out_file, row.names=TRUE)

    # specificity
    options(repr.plot.width=5, repr.plot.height=5)
}

# create a dreamletResult form this list
res.compare = new("dreamletResult", fitList)
write.csv(drop_log, file.path(output_dir, "dropped_samples_log.csv"), row.names = FALSE)
saveRDS(res.compare, file=file.path(output_dir, 'one_vs_all_fit_dreamletResult.rds'))
