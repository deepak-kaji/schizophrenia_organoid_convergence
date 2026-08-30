library(BiocParallel)
library(SingleCellExperiment)
library(dreamlet)
library(mashr)
library(dplyr)
library(parallel)
library(doParallel)
library(foreach)
library(progressr)

# ============================================================
# 1. READ IN DATA
# ============================================================

CSV_PATH <- "/mnt/sdb/scz_meta_analysis_processed/dge_signatures/dreamlet_dges/disease_analyses/"

dge_purcell <- read.csv(paste0(CSV_PATH, "dreamlet_purcell_dge_scz_results.csv"))
dge_walsh_15q13 <- read.csv(paste0(CSV_PATH, "dreamlet_walsh_15q13_dge_scz_results.csv"))
dge_nrxn1 <- read.csv(paste0(CSV_PATH, "dreamlet_nrxn1_combined.csv"))
dge_22q11 <- read.csv(paste0(CSV_PATH, "dreamlet_22q11_combined.csv"))
dge_idiopathic <- read.csv(paste0(CSV_PATH, "dreamlet_idiopathic_combined.csv"))

# ============================================================
# 2. MAKE DREAMLET OUTPUTS MATCH META-ANALYSIS STRUCTURE
# ============================================================

dge_purcell$std.error <- abs(dge_purcell$logFC / dge_purcell$t)
dge_walsh_15q13$std.error <- abs(dge_walsh_15q13$logFC / dge_walsh_15q13$t)

dge_purcell$estimate <- dge_purcell$logFC
dge_walsh_15q13$estimate <- dge_walsh_15q13$logFC

dge_purcell$p.value <- dge_purcell$P.Value
dge_walsh_15q13$p.value <- dge_walsh_15q13$P.Value

# ============================================================
# 3. COMMON GENES
# ============================================================

common_genes <- Reduce(intersect, list(
  dge_22q11$ID,
  dge_nrxn1$ID,
  dge_idiopathic$ID,
  dge_purcell$ID,
  dge_walsh_15q13$ID
))

cat("Number of common genes:", length(common_genes), "\n")

# ============================================================
# 4. CELL TYPES
# ============================================================

celltypes <- unique(c(
  dge_22q11$assay,
  dge_nrxn1$assay,
  dge_idiopathic$assay,
  dge_purcell$assay,
  dge_walsh_15q13$assay
))

# ============================================================
# 5. HELPER FUNCTION
# ============================================================

make_vector <- function(df, genes, value_col) {
  out <- rep(NA_real_, length(genes))
  names(out) <- genes
  idx <- match(df$ID, genes)
  keep <- !is.na(idx)
  out[idx[keep]] <- df[[value_col]][keep]
  as.numeric(out)
}

# ============================================================
# 6. CONSTRUCT FULL CONDITION MATRICES
# ============================================================

conditions_B <- list()
conditions_S <- list()

for (ct in celltypes) {

  d22_ct <- subset(dge_22q11, assay == ct)
  dnr_ct <- subset(dge_nrxn1, assay == ct)
  did_ct <- subset(dge_idiopathic, assay == ct)
  d3q_ct <- subset(dge_purcell, assay == ct)
  d15_ct <- subset(dge_walsh_15q13, assay == ct)

  Bhat_ct <- cbind(
    idiopathic = make_vector(did_ct, common_genes, "estimate"),
    nrxn1 = make_vector(dnr_ct, common_genes, "estimate"),
    del22q11 = make_vector(d22_ct, common_genes, "estimate"),
    del3q29 = make_vector(d3q_ct, common_genes, "estimate"),
    del15q13 = make_vector(d15_ct, common_genes, "estimate")
  )

  Shat_ct <- cbind(
    idiopathic = make_vector(did_ct, common_genes, "std.error"),
    nrxn1 = make_vector(dnr_ct, common_genes, "std.error"),
    del22q11 = make_vector(d22_ct, common_genes, "std.error"),
    del3q29 = make_vector(d3q_ct, common_genes, "std.error"),
    del15q13 = make_vector(d15_ct, common_genes, "std.error")
  )

  colnames(Bhat_ct) <- paste0(colnames(Bhat_ct), "_", ct)
  colnames(Shat_ct) <- paste0(colnames(Shat_ct), "_", ct)

  conditions_B[[ct]] <- Bhat_ct
  conditions_S[[ct]] <- Shat_ct
}

Bhat <- do.call(cbind, conditions_B)
Shat <- do.call(cbind, conditions_S)

rownames(Bhat) <- common_genes
rownames(Shat) <- common_genes

cat("Bhat dimensions:", dim(Bhat), "\n")
cat("Shat dimensions:", dim(Shat), "\n")

# ============================================================
# 7. CLEAN DATA EXACTLY AS IN REAL ANALYSIS
# ============================================================

good <- rowSums(!is.na(Bhat)) > 0

Bhat <- Bhat[good, , drop = FALSE]
Shat <- Shat[good, , drop = FALSE]

Bhat[!is.finite(Bhat)] <- NA
Shat[!is.finite(Shat)] <- NA
Shat[Shat <= 0] <- NA

stopifnot(identical(rownames(Bhat), rownames(Shat)))
stopifnot(all(dim(Bhat) == dim(Shat)))
stopifnot(identical(colnames(Bhat), colnames(Shat)))

keep_cols <- colSums(!is.na(Bhat)) > 0 & colSums(!is.na(Shat)) > 0

Bhat <- Bhat[, keep_cols, drop = FALSE]
Shat <- Shat[, keep_cols, drop = FALSE]

stopifnot(identical(rownames(Bhat), rownames(Shat)))
stopifnot(identical(colnames(Bhat), colnames(Shat)))

cat("Final Bhat dimensions:", dim(Bhat), "\n")
cat("Final Shat dimensions:", dim(Shat), "\n")

# ============================================================
# 8. DISEASES AND CELL TYPES
# ============================================================

disease_names <- c(
  "idiopathic",
  "nrxn1",
  "del22q11",
  "del3q29",
  "del15q13"
)

column_metadata <- data.frame(
  column = colnames(Bhat),
  stringsAsFactors = FALSE
)

column_metadata$disease <- sub("_.*$", "", column_metadata$column)
column_metadata$celltype <- sub("^[^_]+_", "", column_metadata$column)

valid_celltypes <- unique(column_metadata$celltype)

valid_celltypes <- valid_celltypes[
  sapply(valid_celltypes, function(ct) {
    all(disease_names %in% column_metadata$disease[column_metadata$celltype == ct])
  })
]

cat("Number of cell types with all 5 diseases:", length(valid_celltypes), "\n")
print(valid_celltypes)

cat("Number of usable disease-celltype conditions:", ncol(Bhat), "\n")
cat("Number of complete 5-disease cell types:", length(valid_celltypes), "\n")
cat("Number of complete 5-disease conditions:", length(valid_celltypes) * length(disease_names), "\n")

# ============================================================
# 9. SIMULATION PARAMETERS
# ============================================================

shared_fraction_range <- c(0.25, 0.50, 0.75)
#logFC_range <- c(0.25, 0.50, 0.75, 1.00, 1.25, 1.50, 2.00)
logFC_range <- c(0.1, 0.3, 0.50, 0.75, 1.00, 1.25, 1.50, 2.00, 4.00)

n_true_genes <- 500
n_replicates <- 10

# ============================================================
# 10. PRECOMPUTE OBSERVED MASH COVARIANCE STRUCTURE
# ============================================================

cat("\nEstimating observed MASH covariance structure...\n")

real_mash_data <- mash_set_data(Bhat, Shat)

U_c <- cov_canonical(real_mash_data)

cat("Canonical covariance matrices:", length(U_c), "\n")

# ============================================================
# 11. POWER SIMULATION FUNCTION
# ============================================================

run_mash_power_simulation <- function(
  logFC,
  shared_fraction,
  n_true_genes,
  Bhat,
  Shat,
  U_c,
  seed = NULL
) {

  if (!is.null(seed)) set.seed(seed)

  genes <- rownames(Bhat)
  n_genes <- length(genes)
  n_conditions <- ncol(Bhat)

  if (is.null(genes) || length(genes) == 0) {
    stop("Bhat has no gene rownames.")
  }

  if (n_true_genes > n_genes) {
    stop(
      paste0(
        "n_true_genes = ",
        n_true_genes,
        " exceeds the available gene universe of ",
        n_genes,
        "."
      )
    )
  }

  if (shared_fraction < 0 || shared_fraction > 1) {
    stop("shared_fraction must be between 0 and 1.")
  }

  if (logFC < 0) {
    stop("logFC must be non-negative.")
  }

  # ----------------------------------------------------------
  # Determine number of shared and disease-specific genes
  # ----------------------------------------------------------

  n_shared <- round(n_true_genes * shared_fraction)

  n_specific_total <- n_true_genes - n_shared

  specific_counts <- rep(
    floor(n_specific_total / length(disease_names)),
    length(disease_names)
  )

  remainder <- n_specific_total - sum(specific_counts)

  if (remainder > 0) {
    specific_counts[seq_len(remainder)] <-
      specific_counts[seq_len(remainder)] + 1
  }

  names(specific_counts) <- disease_names

  # ----------------------------------------------------------
  # Select exactly n_true_genes unique genes
  # ----------------------------------------------------------

  selected_genes <- sample(
    genes,
    size = n_true_genes,
    replace = FALSE
  )

  shared_genes <- if (n_shared > 0) {
    selected_genes[seq_len(n_shared)]
  } else {
    character(0)
  }

  specific_genes <- list()

  cursor <- n_shared + 1

  for (d in disease_names) {

    n_d <- specific_counts[d]

    if (n_d > 0) {

      specific_genes[[d]] <-
        selected_genes[cursor:(cursor + n_d - 1)]

      cursor <- cursor + n_d

    } else {

      specific_genes[[d]] <- character(0)

    }
  }

  # ----------------------------------------------------------
  # TRUE EFFECT MATRIX
  # ----------------------------------------------------------

  true_Bhat <- matrix(
    0,
    nrow = n_genes,
    ncol = n_conditions,
    dimnames = list(
      genes,
      colnames(Bhat)
    )
  )

  # ----------------------------------------------------------
  # Shared effects
  #
  # Shared genes have the same effect in all five diseases
  # within every cell type where all five diseases exist.
  # ----------------------------------------------------------

  for (ct in valid_celltypes) {

    cols_ct <- paste0(
      disease_names,
      "_",
      ct
    )

    cols_ct <- cols_ct[
      cols_ct %in% colnames(true_Bhat)
    ]

    if (
      length(shared_genes) > 0 &&
      length(cols_ct) > 0
    ) {

      true_Bhat[
        shared_genes,
        cols_ct
      ] <- logFC

    }
  }

  # ----------------------------------------------------------
  # Disease-specific effects
  # ----------------------------------------------------------

  for (d in disease_names) {

    genes_d <- specific_genes[[d]]

    if (length(genes_d) == 0) next

    for (ct in valid_celltypes) {

      col_d <- paste0(
        d,
        "_",
        ct
      )

      if (col_d %in% colnames(true_Bhat)) {

        true_Bhat[
          genes_d,
          col_d
        ] <- logFC

      }
    }
  }

  # ----------------------------------------------------------
  # SIMULATE OBSERVED BHAT
  #
  # Preserve the EXACT observed Shat structure.
  #
  # Observed:
  #   Bhat_sim ~ Normal(true effect, observed SE)
  #
  # Missing:
  #   Bhat_sim = NA
  # ----------------------------------------------------------

  sim_Bhat <- matrix(
    NA_real_,
    nrow = n_genes,
    ncol = n_conditions,
    dimnames = list(
      genes,
      colnames(Bhat)
    )
  )

  for (j in seq_len(n_conditions)) {

    observed <- is.finite(Shat[, j]) & Shat[, j] > 0

    if (any(observed)) {

      sim_Bhat[observed, j] <-
        rnorm(
          sum(observed),
          mean = true_Bhat[observed, j],
          sd = Shat[observed, j]
        )

    }
  }

  # ----------------------------------------------------------
  # MASH FIT
  # ----------------------------------------------------------

  sim_data <- mash_set_data(
    sim_Bhat,
    Shat
  )

  sim_fit <- mash(
    sim_data,
    Ulist = U_c
  )

  # ----------------------------------------------------------
  # POSTERIOR QUANTITIES
  # ----------------------------------------------------------

  sim_lfsr <- get_lfsr(sim_fit)
  sim_pm <- get_pm(sim_fit)

  # ----------------------------------------------------------
  # COMPOSITE POSTERIOR TEST
  #
  # For each complete cell type:
  #
  #   test whether ALL five diseases
  #   show posterior support.
  # ----------------------------------------------------------

  celltype_results <- list()

  for (ct in valid_celltypes) {

    cols_ct <- paste0(
      disease_names,
      "_",
      ct
    )

    cols_ct <- cols_ct[
      cols_ct %in% colnames(sim_lfsr)
    ]

    support_ct <- 1 - sim_lfsr[
      ,
      cols_ct,
      drop = FALSE
    ]

    composite_score <- compositePosteriorTest(
      support_ct,
      include = colnames(support_ct),
      test = "all"
    )

    names(composite_score) <-
      rownames(sim_lfsr)

    celltype_results[[ct]] <-
      composite_score
  }

  # ----------------------------------------------------------
  # RECOVERY OF SHARED GENES
  # ----------------------------------------------------------

  shared_recovery_by_celltype <-
    numeric(length(valid_celltypes))

  names(shared_recovery_by_celltype) <-
    valid_celltypes

  n_shared_recovered_by_celltype <-
    numeric(length(valid_celltypes))

  names(n_shared_recovered_by_celltype) <-
    valid_celltypes

  for (ct in valid_celltypes) {

    scores <-
      celltype_results[[ct]][shared_genes]

    recovered <-
      scores > 0.95

    n_shared_recovered_by_celltype[ct] <-
      sum(
        recovered,
        na.rm = TRUE
      )

    shared_recovery_by_celltype[ct] <-
      mean(
        recovered,
        na.rm = TRUE
      )
  }

  # ----------------------------------------------------------
  # GLOBAL SHARED RECOVERY
  # ----------------------------------------------------------

  all_shared_scores <- unlist(
    lapply(
      celltype_results,
      function(x) x[shared_genes]
    )
  )

  n_shared_recovered <-
    sum(
      all_shared_scores > 0.95,
      na.rm = TRUE
    )

  total_shared_tests <-
    length(shared_genes) *
    length(valid_celltypes)

  shared_recovery <-
    if (total_shared_tests > 0) {
      n_shared_recovered /
        total_shared_tests
    } else {
      NA_real_
    }

  # ----------------------------------------------------------
  # PER-CONDITION SIGNIFICANCE
  # ----------------------------------------------------------

  n_sig <- colSums(
    sim_lfsr < 0.05,
    na.rm = TRUE
  )

  # ----------------------------------------------------------
  # RETURN RESULTS
  # ----------------------------------------------------------

  result <- data.frame(
    logFC = logFC,
    shared_fraction = shared_fraction,
    n_genes_universe = n_genes,
    n_conditions = n_conditions,
    n_complete_celltypes = length(valid_celltypes),
    n_true_genes = n_true_genes,
    n_true_shared = n_shared,
    n_true_specific = n_specific_total,
    n_shared_tests = total_shared_tests,
    n_shared_recovered = n_shared_recovered,
    shared_recovery = shared_recovery,
    mean_shared_recovery_across_celltypes = mean(
      shared_recovery_by_celltype,
      na.rm = TRUE
    ),
    seed = seed,
    simulation_status = "success",
    stringsAsFactors = FALSE
  )

  # ----------------------------------------------------------
  # Add condition-specific significant gene counts
  # ----------------------------------------------------------

  for (nm in names(n_sig)) {

    result[[paste0("n_sig_", nm)]] <- n_sig[nm]

  }

  # ----------------------------------------------------------
  # Add cell-type-specific shared recovery
  # ----------------------------------------------------------

  for (ct in valid_celltypes) {

    safe_name <- gsub("[^A-Za-z0-9]+", "_", ct)

    result[[paste0("shared_recovery_", safe_name)]] <- shared_recovery_by_celltype[ct]

  }
  return(result)
    
}

# # ============================================================
# # 12. TEST ONE SIMULATION
# # ============================================================

# cat("\nTesting one full MASH simulation...\n")

# start_time <- Sys.time()

# test_result <- run_mash_power_simulation(
#   logFC = 1.0,
#   shared_fraction = 0.50,
#   n_true_genes = n_true_genes,
#   Bhat = Bhat,
#   Shat = Shat,
#   U_c = U_c,
#   seed = 1
# )

# print(test_result)

# cat(
#   "\nSingle simulation runtime:",
#   round(
#     as.numeric(
#       difftime(
#         Sys.time(),
#         start_time,
#         units = "mins"
#       )
#     ),
#     2
#   ),
#   "minutes\n"
# )

# ============================================================
# 13. BUILD SIMULATION GRID
# ============================================================

simulation_grid <- expand.grid(
  shared_fraction = shared_fraction_range,
  logFC = logFC_range,
  replicate = seq_len(n_replicates),
  stringsAsFactors = FALSE
)

rownames(simulation_grid) <- NULL

cat(
  "\nNumber of simulations:",
  nrow(simulation_grid),
  "\n"
)

outdir <- "/mnt/sdb/scz_meta_analysis_processed/mash_results/power_simulation_45condition/"

dir.create(
  outdir,
  recursive = TRUE,
  showWarnings = FALSE
)

# ============================================================
# 14. PARALLEL SIMULATION WITH LIVE PROGRESS
# ============================================================

start_time <- Sys.time()

n_cores <- 25
total <- nrow(simulation_grid)

cl <- makeCluster(n_cores)
registerDoParallel(cl)

# ------------------------------------------------------------
# Progress-aware combine function
# ------------------------------------------------------------

completed <- 0

progress_combine <- function(x, y) {

  completed <<- completed + 1

  elapsed <- as.numeric(
    difftime(
      Sys.time(),
      start_time,
      units = "secs"
    )
  )

  rate <- completed / elapsed
  remaining <- total - completed

  eta_seconds <- if (rate > 0) remaining / rate else NA_real_

  cat(
    sprintf(
      "\rCompleted: %d/%d (%.1f%%) | %.3f sims/sec | ETA: %.1f min",
      completed,
      total,
      100 * completed / total,
      rate,
      eta_seconds / 60
    ),
    flush = TRUE
  )

  rbind(x, y)
}


# ------------------------------------------------------------
# Run simulations
# ------------------------------------------------------------

results <- foreach(
  i = seq_len(total),
  .packages = c("mashr", "dreamlet"),
  .combine = progress_combine,
  .multicombine = FALSE
) %dopar% {

  row <- simulation_grid[i, ]

  run_mash_power_simulation(
    logFC = row$logFC,
    shared_fraction = row$shared_fraction,
    n_true_genes = n_true_genes,
    Bhat = Bhat,
    Shat = Shat,
    U_c = U_c,
    seed = i
  )
}

stopCluster(cl)

cat("\n\nSimulation complete.\n")

elapsed_minutes <- as.numeric(
  difftime(
    Sys.time(),
    start_time,
    units = "mins"
  )
)

cat(
  "Total runtime:",
  round(elapsed_minutes, 2),
  "minutes\n"
)

# ============================================================
# 15. SAVE RAW RESULTS
# ============================================================

outdir <- "/mnt/sdb/scz_meta_analysis_processed/mash_results/power_simulation_45condition/"

dir.create(
  outdir,
  recursive = TRUE,
  showWarnings = FALSE
)

write.csv(
  results,
  file.path(
    outdir,
    "mash_power_simulation_45condition_results.csv"
  ),
  row.names = FALSE
)

# ============================================================
# 16. SUMMARIZE POWER
# ============================================================

summary_results <- results %>%
  group_by(
    shared_fraction,
    logFC
  ) %>%
  summarise(
    n_success = sum(
      simulation_status == "success"
    ),
    mean_shared_recovery = mean(
      shared_recovery,
      na.rm = TRUE
    ),
    sd_shared_recovery = sd(
      shared_recovery,
      na.rm = TRUE
    ),
    mean_shared_recovered = mean(
      n_shared_recovered,
      na.rm = TRUE
    ),
    sd_shared_recovered = sd(
      n_shared_recovered,
      na.rm = TRUE
    ),
    mean_recovery_across_celltypes = mean(
      mean_shared_recovery_across_celltypes,
      na.rm = TRUE
    ),
    sd_recovery_across_celltypes = sd(
      mean_shared_recovery_across_celltypes,
      na.rm = TRUE
    ),
    n_replicates = n(),
    .groups = "drop"
  )

write.csv(
  summary_results,
  file.path(
    outdir,
    "mash_power_simulation_45condition_summary.csv"
  ),
  row.names = FALSE
)

# ============================================================
# 17. SAVE SIMULATION SETTINGS
# ============================================================

writeLines(
  c(
    paste(
      "Number of genes:",
      nrow(Bhat)
    ),
    paste(
      "Number of MASH conditions:",
      ncol(Bhat)
    ),
    paste(
      "Number of diseases:",
      length(disease_names)
    ),
    paste(
      "Number of complete cell types:",
      length(valid_celltypes)
    ),
    paste(
      "Number of complete disease-celltype conditions:",
      length(valid_celltypes) *
        length(disease_names)
    ),
    paste(
      "True genes:",
      n_true_genes
    ),
    paste(
      "Shared fractions:",
      paste(
        shared_fraction_range,
        collapse = ", "
      )
    ),
    paste(
      "LogFC values:",
      paste(
        logFC_range,
        collapse = ", "
      )
    ),
    paste(
      "Replicates:",
      n_replicates
    ),
    paste(
      "MASH covariance:",
      "Observed-data canonical covariance U.c"
    ),
    paste(
      "Missingness:",
      "Preserved from observed Shat"
    ),
    paste(
      "True shared effects:",
      "Identical magnitude and direction across all 5 diseases and all complete cell types"
    ),
    paste(
      "True null effects:",
      "Exactly zero"
    ),
    paste(
      "True disease-specific effects:",
      "Identical magnitude and direction within disease across all complete cell types"
    ),
    paste(
      "MASH fitting:",
      "mash(..., Ulist = U_c)"
    )
  ),
  file.path(
    outdir,
    "simulation_settings.txt"
  )
)

# ============================================================
# 18. FINAL OUTPUT
# ============================================================

cat(
  "\n============================================\n"
)

cat(
  "MASH POWER SIMULATION COMPLETE\n"
)

cat(
  "============================================\n"
)

cat(
  "Raw results:",
  file.path(
    outdir,
    "mash_power_simulation_45condition_results.csv"
  ),
  "\n"
)

cat(
  "Summary:",
  file.path(
    outdir,
    "mash_power_simulation_45condition_summary.csv"
  ),
  "\n"
)

cat(
  "Successful simulations:",
  sum(
    results$simulation_status == "success"
  ),
  "\n"
)

print(summary_results)
