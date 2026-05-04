library(zellkonverter)
library(SingleCellExperiment)
library(Matrix)
library(dplyr)
library(lme4)
library(broom.mixed)

# ----------------------------
# Load data
# ----------------------------
sce <- readH5AD('/mnt/sdb/scz_meta_analysis_processed/anndata_objs/dorsal_neurogenesis_dreamlet.h5ad',
                use_hdf5 = TRUE, layers = FALSE, raw = FALSE, verbose = FALSE, uns = FALSE)

# ----------------------------
# Metadata cleanup
# ----------------------------

sce$Donor_Sample <- paste(sce$Donor, sce$Sample.Name, sep = "_")

sce$Path <- factor(
  sce$Path,
  levels = c(
    "25",
    "30"
  ),
  labels = c("Indirect", "Direct")
)

sce$Broad_Genotype <- factor(make.names(sce$Broad_Genotype))
sce$Broad_Genotype <- relevel(sce$Broad_Genotype, ref = "Control")

# ----------------------------
# Convert to data.frame
# ----------------------------
df <- as.data.frame(colData(sce))

# ----------------------------
# Pseudobulk at Donor_Sample level
# ----------------------------
pb <- df %>%
  group_by(
    Donor_Sample,
    Donor,
    Broad_Genotype,
    Protocol,
    Chemistry,
    Manuscript,
    Sex
  ) %>%
  summarise(
    n_total = n(),
    n_Direct = sum(Path == "Direct"),
    n_Indirect = sum(Path == "Indirect"),
    n_counts = mean(n_counts),
    percent_mito = mean(percent_mito),
    .groups = "drop"
  )

# ----------------------------
# Scale covariates
# ----------------------------
pb$n_counts_scaled <- scale(pb$n_counts)
pb$percent_mito_scaled <- scale(pb$percent_mito)

# ----------------------------
# Fit binomial mixed model
# ----------------------------
model_pb <- glmer(
  cbind(n_Direct, n_Indirect) ~
    Broad_Genotype +
    n_counts_scaled +
    percent_mito_scaled +
    (1 | Sex) +
    (1 | Chemistry) +
    (1 | Manuscript) +
    (1 | Donor),

  data = pb,
  family = binomial(link = "logit"),
  control = glmerControl(
    optimizer = "bobyqa",
    optCtrl = list(maxfun = 2e5)
  )
)

# ----------------------------
# Output
# ----------------------------
summary(model_pb)

## Save Output ##

coef_df <- broom.mixed::tidy(model_pb, effects = "fixed", conf.int = TRUE)
coef_df_or <- coef_df %>% mutate(odds_ratio = exp(estimate), conf.low.or = exp(conf.low), conf.high.or = exp(conf.high))
coef_table <- coef_df_or %>% filter(term != "(Intercept)") %>% arrange(p.value)
write.csv(coef_table, "/mnt/sdb/scz_meta_analysis_processed/trajectory_analyses/lmer_di_pseudobulk_coefficients.csv", row.names = FALSE)
