library(zellkonverter)
library(SingleCellExperiment)
library(Matrix)
library(limma)
library(dplyr)
library(lme4)
library(lmerTest)
library(purrr)
library(metafor)

# ----------------------------
# Load data
# ----------------------------
sce <- readH5AD('/mnt/sdb/scz_meta_analysis_processed/anndata_objs/dorsal_neurogenesis_dreamlet.h5ad',
                use_hdf5 = TRUE, layers = FALSE, raw = FALSE, verbose = FALSE, uns = FALSE)

# ----------------------------
# Metadata cleanup
# ----------------------------

sce$Donor_Sample <- paste(sce$Donor, sce$Sample.Name, sep = "_")

sce$Broad_Genotype <- factor(make.names(sce$Broad_Genotype))
sce$Broad_Genotype <- relevel(sce$Broad_Genotype, ref = "Control")

## Extract Manuscripts from ColData As List ##
pt <- as.data.frame(colData(sce))

pt$Donor_Sample <- paste(pt$Donor, pt$Sample.Name, sep = "_")

pt_summary <- aggregate(
    palantir_pseudotime ~ Manuscript + Donor_Sample + Donor + Sample.Name + Broad_Genotype + Day + Chemistry + Protocol + Sex,
    data = pt,
    FUN = mean
)

pt_by_study <- split(pt_summary, pt_summary$Manuscript)

## Define Function To Fit Linear Mixed Models to Individual Manuscripts ##

safe_lmer_fit <- function(df, form,
                          re_candidates = c("Donor", "Sample.Name", "Sex", "Protocol", "Chemistry"),
                          min_obs_per_level = 2) {

    df <- df[is.finite(df$pt_logit), ]
    df <- droplevels(df)

    if (nrow(df) < 5) return(NULL)

    # ----------------------------
    # 1. validate random effects
    # ----------------------------
    valid_re <- function(var) {

        if (!var %in% names(df)) return(NULL)

        x <- factor(df[[var]])
        tab <- table(x)

        if (length(tab) < 2) return(NULL)
        if (any(tab < min_obs_per_level)) return(NULL)

        paste0("(1|", var, ")")
    }

    re_terms <- Filter(Negate(is.null),
                       lapply(re_candidates, valid_re))

    # ----------------------------
    # 2. build formula
    # ----------------------------
    if (length(re_terms) > 0) {

        form2 <- update(form,
                        paste(". ~ . +", paste(re_terms, collapse = " + ")))

        message("FORMULA (lmer): ", deparse(form2))

        fit <- tryCatch(
            lmer(form2, data = df, REML = TRUE),
            error = function(e) {
                message("lmer failed: ", conditionMessage(e))
                NULL
            }
        )

    } else {

        message("FORMULA (lm fallback): ", deparse(form))

        fit <- tryCatch(
            lm(form, data = df),
            error = function(e) {
                message("lm failed: ", conditionMessage(e))
                NULL
            }
        )
    }

    return(fit)
}

## Apply Function To Each Manuscript ##
lmm_results <- lapply(names(pt_by_study), function(study) {

    df <- pt_by_study[[study]]

    # Reveal Manuscript #
    print(df[,1][[1]])

    eps <- 1e-4
    df$pt_logit <- qlogis(pmin(pmax(df$palantir_pseudotime, eps), 1 - eps))

    use_day <- length(unique(df$Day)) > 1

    # build formula conditionally
    
    if (use_day) {
      form <- pt_logit ~ Broad_Genotype + Day
    } else {
      form <- pt_logit ~ Broad_Genotype
    }

    # fit model
    fit <- safe_lmer_fit(df, form) #, data = df, REML = FALSE)
    if (is.null(fit)) return(NULL)

    coefs <- summary(fit)$coefficients
    coefs <- coefs[grep("^Broad_Genotype", rownames(coefs)), , drop = FALSE]

    data.frame(
        term = rownames(coefs),
        beta = coefs[, "Estimate"],
        se = coefs[, "Std. Error"],
        t = coefs[, "t value"],
        p = coefs[, "Pr(>|t|)"],
        Manuscript = study
    )
})


meta_df <- bind_rows(lmm_results)
rownames(meta_df) <- NULL
meta_df <- meta_df %>%
  filter(!is.na(beta), !is.na(se)) %>%
  mutate(
    vi = se^2
  )


## Prep for Metafor ##

walsh_df <- meta_df %>%
  filter(Manuscript == "Walsh") %>%
  mutate(genotype = case_when(
    grepl("22q11", term) ~ "22q11",
    grepl("15q13", term) ~ "15q13",
    TRUE ~ NA_character_
  ))

df_22q11 <- meta_df %>%
  filter(
    Manuscript %in% c("Rao", "Shin") |
    (Manuscript == "Walsh" & grepl("22q11", term))
  ) %>%
  mutate(vi = se^2)

df_nrxn1 <- meta_df %>%
  filter(
    Manuscript %in% c("Sebastian", "Fernando")
  ) %>%
  mutate(vi = se^2)

df_idio <- meta_df %>%
  filter(
    Manuscript %in% c("Notaras", "Sawada")
  ) %>%
  mutate(vi = se^2)

## Run Metafor ##

res_22q11 <- rma(
  yi = beta,
  vi = vi,
  data = df_22q11,
  method = "FE"
)

res_nrxn1 <- rma(
  yi = beta,
  vi = vi,
  data = df_nrxn1,
  method = "FE"
)

res_idio <- rma(
  yi = beta,
  vi = vi,
  data = df_idio,
  method = "FE"
)

## Extract CSVs ##

write.csv(meta_df, "/mnt/sdb/scz_meta_analysis_processed/trajectory_analyses/meta_df.csv", row.names = FALSE)

summary_idio <- data.frame(
    estimate = res_idio$b,
    se = res_idio$se,
    z = res_idio$zval,
    p = res_idio$pval,
    ci.lb = res_idio$ci.lb,
    ci.ub = res_idio$ci.ub,
    tau2 = res_idio$tau2,
    I2 = res_idio$I2,
    H2 = res_idio$H2
)

write.csv(summary_idio, "/mnt/sdb/scz_meta_analysis_processed/trajectory_analyses/meta_idio.csv", row.names = FALSE)

summary_22q11 <- data.frame(
    estimate = res_22q11$b,
    se = res_22q11$se,
    z = res_22q11$zval,
    p = res_22q11$pval,
    ci.lb = res_22q11$ci.lb,
    ci.ub = res_22q11$ci.ub,
    tau2 = res_22q11$tau2,
    I2 = res_22q11$I2,
    H2 = res_22q11$H2
)

write.csv(summary_22q11, "/mnt/sdb/scz_meta_analysis_processed/trajectory_analyses/meta_22q11.csv", row.names = FALSE)

summary_nrxn1 <- data.frame(
    estimate = res_nrxn1$b,
    se = res_nrxn1$se,
    z = res_nrxn1$zval,
    p = res_nrxn1$pval,
    ci.lb = res_nrxn1$ci.lb,
    ci.ub = res_nrxn1$ci.ub,
    tau2 = res_nrxn1$tau2,
    I2 = res_nrxn1$I2,
    H2 = res_nrxn1$H2
)

write.csv(summary_nrxn1, "/mnt/sdb/scz_meta_analysis_processed/trajectory_analyses/meta_nrxn1.csv", row.names = FALSE)

