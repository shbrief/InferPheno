## Bundle pre-fit supervised heads into inst/extdata/ so the
## installed package can serve `inferpheno()` without retraining.
##
## This is the Aim 4 step in the targets pipeline. It depends on the
## aim1 benchmark targets having already populated `results/aim1_*.rds`
## with the fitted full-data models (not just OOF predictions).
##
## Run from the InferPhenoCuration project root (sibling of the
## InferPheno package directory):
##   targets::tar_make(inferpheno_extdata)
## or directly:
##   source("../InferPheno/inst/scripts/bundle_heads.R")

suppressPackageStartupMessages({
  library(InferPheno)
})

extdata_dir <- file.path("..", "InferPheno", "inst", "extdata")
dir.create(extdata_dir, showWarnings = FALSE, recursive = TRUE)

bundle_one <- function(rds_name, phenotype, kind, rav_cols, model,
                      training_meta) {
  obj <- methods::new("RAVPhenoModel",
                      phenotype     = phenotype,
                      rav_prior     = getOption("inferpheno.rav_prior", "C2"),
                      rav_cols      = rav_cols,
                      model         = model,
                      kind          = kind,
                      training_meta = training_meta)
  out <- file.path(extdata_dir, rds_name)
  saveRDS(obj, out, compress = "xz")
  message("Bundled ", phenotype, " -> ", out, " (",
          format(file.info(out)$size, big.mark = ","), " bytes)")
  invisible(out)
}

## --- sex -----------------------------------------------------------------
sex_cv <- readRDS("results/aim1_sex_cv.rds")
sex_metrics <- readRDS("results/aim1_sex_metrics.rds")
sex_full <- readRDS("results/aim1_sex_full_model.rds")  # full-data RF
bundle_one(
  rds_name = "sex.rds",
  phenotype = "sex",
  kind = "binary",
  rav_cols = paste0("RAV", c(1076, 1241, 152, 220, 221, 2514)),
  model = sex_full,
  training_meta = list(
    cohorts   = c("BRCA", "UCEC", "OV", "PRAD"),
    n_train   = length(sex_cv$rf$truth),
    auc_oof   = sex_metrics$rf$auc,
    source    = "TCGA via curatedTCGAData v2.1.1"
  )
)

## --- msi_status ----------------------------------------------------------
msi_status_cv <- readRDS("results/aim1_msi_status_cv.rds")
msi_status_metrics <- readRDS("results/aim1_msi_status_metrics.rds")
msi_status_full <- readRDS("results/aim1_msi_status_full_model.rds")
bundle_one(
  rds_name = "msi_status.rds",
  phenotype = "msi_status",
  kind = "binary",
  rav_cols = msi_status_full$rav_cols %||%
    rownames(msi_status_full$importance),
  model = msi_status_full,
  training_meta = list(
    cohorts   = c("COAD", "STAD", "UCEC"),
    n_train   = length(msi_status_cv$rf$truth),
    auc_oof   = msi_status_metrics$rf$auc,
    source    = "TCGA combined MAE"
  )
)

## --- msi_score (regression) ----------------------------------------------
msi_score_cv <- readRDS("results/aim1_msi_score_cv.rds")
msi_score_metrics <- readRDS("results/aim1_msi_score_metrics.rds")
msi_score_full <- readRDS("MSI_score_prediction/results/msi_score_vignette_rfmodel.rds")
bundle_one(
  rds_name = "msi_score.rds",
  phenotype = "msi_score",
  kind = "regression",
  rav_cols = paste0("RAV", c(1166, 192, 178, 188)),
  model = msi_score_full,
  training_meta = list(
    cohort    = "COAD_READ",
    n_train   = length(msi_score_cv$truth),
    rmse_oof  = msi_score_metrics$rmse,
    r2_oof    = msi_score_metrics$r2,
    spearman_oof = msi_score_metrics$spearman,
    source    = "cBioPortal coadread_tcga_pan_can_atlas_2018"
  )
)

## --- cms_subtype (multiclass) -------------------------------------------
if (file.exists("results/aim1_crc_subtype_full_model.rds")) {
  crc_full    <- readRDS("results/aim1_crc_subtype_full_model.rds")
  crc_cv      <- readRDS("results/aim1_crc_subtype_cv.rds")
  crc_metrics <- readRDS("results/aim1_crc_subtype_metrics.rds")
  bundle_one(
    rds_name = "cms_subtype.rds",
    phenotype = "cms_subtype",
    kind = "multiclass",
    rav_cols = paste0("RAV", c(188, 1575, 834, 832, 833, 438, 324, 192,
                                981, 220)),
    model = crc_full,
    training_meta = list(
      cohorts   = c("GSE13067", "GSE13294", "GSE14333", "GSE17536",
                    "GSE17537", "GSE2109", "GSE33113", "GSE39582",
                    "TCGA.RNASeqV2"),
      n_train   = length(crc_cv$truth),
      accuracy_oof = crc_metrics$rf_rav$accuracy,
      macro_f1_oof = crc_metrics$rf_rav$macro_f1,
      source    = paste("curatedCRCData 2.34.0 eSets; CMS labels from",
                        "Synapse syn4978511 (Guinney 2015 CRCSC",
                        "consensus); see",
                        "CRC_subtype_pred/data/eSets/cms_label_source.txt")
    )
  )
} else {
  message("Skipping cms_subtype: aim1_crc_subtype is still blocked ",
          "(training eSets/ missing).")
}

`%||%` <- function(a, b) if (is.null(a)) b else a
