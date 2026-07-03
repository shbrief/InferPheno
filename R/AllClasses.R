#' RAVPhenoModel
#'
#' A supervised head fit on top of a fixed RAV basis. Used to predict one
#' phenotype from new gene expression data without retraining the basis.
#'
#' @slot phenotype character. The phenotype name (e.g., "msi_status").
#' @slot rav_prior character. The GSS prior used for the RAV basis (e.g., "C2").
#' @slot rav_cols character. The RAV columns used as features.
#' @slot model ANY. The fitted classifier/regressor object.
#' @slot kind character. One of "binary", "multiclass", "regression", "cox".
#' @slot training_meta list. Cohort composition, n, source paper.
#'
#' @return An object of class `RAVPhenoModel`. Shipped heads are serialized
#'   instances loaded by [inferpheno()]; construct new ones with
#'   `InferPhenoBuilder::build_pheno_head()`.
#' @examples
#' # Each shipped head is a serialized RAVPhenoModel:
#' m <- readRDS(system.file("extdata", "sex.rds", package = "InferPheno"))
#' m@phenotype
#' m@kind
#' m@rav_cols
#' @export
setClass("RAVPhenoModel",
  representation(
    phenotype     = "character",
    rav_prior     = "character",
    rav_cols      = "character",
    model         = "ANY",
    kind          = "character",
    training_meta = "list"))
