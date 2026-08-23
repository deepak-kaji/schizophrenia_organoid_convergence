library(BiocParallel)
library(SingleCellExperiment)
library(dreamlet)
library(mashr)

# Read in Meta-Analytic Estimates and Dreamlet Results #  

CSV_PATH <- '/mnt/sdb/scz_meta_analysis_processed/dge_signatures/dreamlet_dges/disease_analyses/'

dge_purcell <- read.csv(paste0(CSV_PATH, 'dreamlet_purcell_dge_scz_results.csv'))
#dge_walsh_15q13 <- read.csv(paste0(CSV_PATH, 'dreamlet_walsh_15q13_dge_scz_results.csv'))
dge_nrxn1 <- read.csv(paste0(CSV_PATH, 'dreamlet_nrxn1_combined.csv'))
dge_22q11 <- read.csv(paste0(CSV_PATH, 'dreamlet_22q11_combined.csv'))
dge_idiopathic <- read.csv(paste0(CSV_PATH, 'dreamlet_idiopathic_combined.csv'))

# add se to dreamlet outputs
dge_purcell$std.error <- abs(dge_purcell$logFC / dge_purcell$t)

# make estimate column for dreamlet outputs so they resemble meta_analysis outputs
dge_purcell$estimate <- dge_purcell$logFC

# get intersection fo all genes and get all celltypes

common_genes <- Reduce(intersect, list(dge_22q11$ID, dge_nrxn1$ID, dge_idiopathic$ID,
				       dge_purcell$ID))

celltypes <- unique(c(dge_22q11$assay, dge_nrxn1$assay,
                      dge_idiopathic$assay, dge_purcell$assay))

celltypes <- c(
  "Dorsal Forebrain Neuron FOXG1+EMX1+NEUROG1+",
  "GABAergic GAD1+GAD2+CALB2+",
  "Hindbrain Neurons NR2F2+PBX3+LHX1+",
  "IPC EOMES+NEUROG2+PAX6+",
  "Astrocyte GFAP+AQP4+HOPX+",
  "Mesenchymal-like cells VIM+VCAN+SPARC+",
  "Mixed Neurons FGF12+GRIN2B+CAMK2B+",
  "Proliferative Radial Glia SOX2+HES6+TOP2A+",
  "Radial Glia SOX2+PAX6+FABP7+",
  "Radial Glia SOX2+VIM+FABP7+")

make_vector <- function(df, genes, value_col) {

  out <- rep(NA_real_, length(genes))
  names(out) <- genes

  idx <- match(df$ID, genes)

  keep <- !is.na(idx)

  out[idx[keep]] <- df[[value_col]][keep]

  as.numeric(out)
}

conditions_B <- list()
conditions_S <- list()

for (ct in celltypes) {

    d22_ct <- subset(dge_22q11, assay == ct)
    dnr_ct <- subset(dge_nrxn1, assay == ct)
    did_ct <- subset(dge_idiopathic, assay == ct)
    d3q_ct <- subset(dge_purcell, assay == ct)

    Bhat_ct <- cbind(
        idiopathic = make_vector(did_ct, common_genes, "estimate"),
        nrxn1      = make_vector(dnr_ct, common_genes, "estimate"),
        del22q11   = make_vector(d22_ct, common_genes, "estimate"),
        del3q29    = make_vector(d3q_ct, common_genes, "estimate")
    )

    Shat_ct <- cbind(
        idiopathic = make_vector(did_ct, common_genes, "std.error"),
        nrxn1      = make_vector(dnr_ct, common_genes, "std.error"),
        del22q11   = make_vector(d22_ct, common_genes, "std.error"),
        del3q29    = make_vector(d3q_ct, common_genes, "std.error")
    )

    # rename columns to encode celltype
    colnames(Bhat_ct) <- paste0(colnames(Bhat_ct), "_", ct)
    colnames(Shat_ct) <- paste0(colnames(Shat_ct), "_", ct)

    conditions_B[[ct]] <- Bhat_ct
    conditions_S[[ct]] <- Shat_ct
}
## mashR inference ##

# combine matrices
Bhat <- do.call(cbind, conditions_B)
Shat <- do.call(cbind, conditions_S)

dim(Bhat)
dim(Shat)

# -------------------------------------------------------
# 1. initial gene filter (keep genes with any signal)
# -------------------------------------------------------
good <- rowSums(!is.na(Bhat)) > 0
Bhat <- Bhat[good, ]
Shat <- Shat[good, ]

# -------------------------------------------------------
# 2. enforce numeric validity BEFORE mash
# -------------------------------------------------------

# remove non-finite values
Bhat[!is.finite(Bhat)] <- NA
Shat[!is.finite(Shat)] <- NA

# remove invalid SEs
Shat[Shat <= 0] <- NA

# -------------------------------------------------------
# 3. ensure row alignment stays identical
# -------------------------------------------------------
stopifnot(identical(rownames(Bhat), rownames(Shat)))

# -------------------------------------------------------
# Ensure Alignment And Remove Dead Columns Like Gabriel 
# -------------------------------------------------------
stopifnot(all(dim(Bhat) == dim(Shat)))
stopifnot(rownames(Bhat) == rownames(Shat))
stopifnot(colnames(Bhat) == colnames(Shat))

keep_cols <- colSums(!is.na(Bhat)) > 0 &
             colSums(!is.na(Shat)) > 0

Bhat <- Bhat[, keep_cols]
Shat <- Shat[, keep_cols]
# -------------------------------------------------------
# 5. final sanity checks (important for debugging)
# -------------------------------------------------------
keep_cols <- colSums(!is.na(Bhat)) >= 500 &
             colSums(!is.na(Shat)) >= 500

Bhat <- Bhat[, keep_cols]
Shat <- Shat[, keep_cols]
# -------------------------------------------------------
# 6. mash object
# -------------------------------------------------------
m_data <- mash_set_data(Bhat, Shat)

# canonical covariance
U.c <- cov_canonical(m_data)

# estimate residual correlation
V.em <- mash_estimate_corr_em(m_data, U.c, details = TRUE)
pm_em <- get_pm(V.em$mash.model)
m_model <- V.em$mash.model

## try canonical only ##
data <- mash_set_data(Bhat, Shat)
U.c <- cov_canonical(data)
fit_c <- mash(data, Ulist = U.c)
pm_c <- get_pm(fit_c)

print(cor(as.vector(pm_em), as.vector(pm_c), use = "pairwise.complete.obs"))

# EXPORT MASH RESULTS 

outdir <- "/mnt/sdb/scz_meta_analysis_processed/mash_results/"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

model <- V.em$mash.model

# Core results
pm   <- get_pm(model)
lfsr <- get_lfsr(model)

# Original (pre-mash) effects
orig <- m_model$logFC.original

# Uncertainty / diagnostics
post_sd <- model$result$PosteriorSD
neg_prob <- model$result$NegativeProb

## Defining Composite Posterior Tests ##

support <- 1 - model$result$lfsr
rownames(support) <- common_genes

# test whether all diseases exhibit any shared DGE signature by celltype

cols <- colnames(support)
celltype <- sub("^[^_]+_", "", cols)
disease <- sub("_.*", "", cols)

celltype_groups <- split(cols, celltype)

outdir <- "/mnt/sdb/scz_meta_analysis_processed/mash_results/"

run_celltype_composite <- function(support, cols, name, outdir, test = "at least 1") {

  res <- compositePosteriorTest(
    support,
    include = cols,
    test = test
  )

  df <- data.frame(
    gene = rownames(support),
    score = as.numeric(res)
  )

  write.csv(
    df,
    file.path(outdir, paste0("minus_15q13_composite_", name, "_", test, ".csv")),
    row.names = FALSE
  )

  return(res)
}

results <- list()

# test each celltype to see if theres a convergeant disease signature #

for (ct in names(celltype_groups)) {

  message("Running cell type: ", ct)

  cols <- celltype_groups[[ct]]

  results[[ct]] <- run_celltype_composite(
    support = support,
    cols = cols,
    name = gsub(" ", "_", ct),
    outdir = outdir,
    test = "all"
  )
}

# -----------------------
# Write matrices to disk
# -----------------------
write.csv(pm, file.path(outdir, "minus_15q13_mash_posterior_mean.csv"))
write.csv(lfsr, file.path(outdir, "minus_15q13_mash_lfsr.csv"))
write.csv(orig, file.path(outdir, "minus_15q13_mash_logFC_original.csv"))
write.csv(post_sd, file.path(outdir, "minus_15q13_mash_posterior_sd.csv"))
write.csv(neg_prob, file.path(outdir, "minus_15q13_mash_negative_prob.csv"))

# -----------------------
# Simple summaries
# -----------------------

# genes significant in at least one cell type
sig_any_gene <- apply(lfsr, 1, function(x) any(x < 0.05, na.rm = TRUE))

# number of significant genes per cell type
sig_per_celltype <- colSums(lfsr < 0.05, na.rm = TRUE)

write.csv(
  data.frame(gene = common_genes, sig_any = sig_any_gene),
  file.path(outdir, "minus_15q13_mash_sig_any_gene.csv"),
  row.names = FALSE
)

write.csv(
  data.frame(celltype = names(sig_per_celltype), n_sig = sig_per_celltype),
  file.path(outdir, "minus_15q13_mash_sig_per_celltype.csv"),
  row.names = FALSE
)

# -----------------------
# Metadata for QC in Python
# -----------------------

write.csv(
  data.frame(
    gene = common_genes,
    missing_fraction = rowMeans(is.na(pm))
  ),
  file.path(outdir, "minus_15q13_gene_missingness.csv"),
  row.names = FALSE
)

write.csv(
  data.frame(
    celltype = colnames(pm),
    missing_fraction = colMeans(is.na(pm))
  ),
  file.path(outdir, "minus_15q13_celltype_missingness.csv"),
  row.names = FALSE
)

# Save session info for reproducibility
writeLines(capture.output(sessionInfo()),
           file.path(outdir, "minus_15q13_sessionInfo.txt"))
