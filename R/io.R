## Data loading utilities shared across all aims.

#' Resolve a project-relative path
#'
#' Always use this instead of hard-coded paths. The root is read from the
#' `INFERPHENO_ROOT` environment variable, falling back to
#' `~/GSS/InferPhenoCuration`.
#'
#' @param ... Path components joined with [file.path()].
#' @return A single absolute path string.
#' @export
proj_path <- function(...) {
    root <- Sys.getenv("INFERPHENO_ROOT", unset = "~/GSS/InferPhenoCuration")
    file.path(path.expand(root), ...)
}

#' Load the GSS RAV model used throughout the paper
#'
#' The prior is pinned at the project level (default `"C2"`) so every aim
#' uses the same basis. Override via `options(inferpheno.rav_prior = ...)`.
#'
#' @param prior Character. GSS prior name.
#' @return A `GenomicSuperSignature::RAVmodel` object.
#' @export
load_rav_model <- function(prior = getOption("inferpheno.rav_prior", "C2")) {
    GenomicSuperSignature::getModel(prior, load = TRUE)
}

#' Project an expression matrix onto RAV space
#'
#' @param expr A genes-x-samples numeric matrix or [SummarizedExperiment::SummarizedExperiment].
#' @param rav_model A `RAVmodel` from [load_rav_model()].
#' @return A samples-x-RAVs matrix with hyphens restored in TCGA barcodes.
#' @export
project_to_rav <- function(expr, rav_model) {
    if (methods::is(expr, "SummarizedExperiment")) {
        expr <- SummarizedExperiment::assay(expr)
    }
    score <- GenomicSuperSignature::calculateScore(expr, rav_model)
    rownames(score) <- gsub("\\.", "-", rownames(score))
    score
}

#' Standardize a sample annotation table
#'
#' Coerces a heterogeneous metadata table to the canonical
#' `(sample_id, cohort, phenotype, ...)` shape used throughout the package.
#'
#' @param df A data.frame.
#' @param sample_id Name of the column to use as `sample_id`.
#' @param cohort Name of the column to use as `cohort`.
#' @param phenotype Name of the column to use as `phenotype`.
#' @param ... Extra named arguments mapping new column names to existing
#'   column names in `df` (e.g. `age = "AGE_AT_DX"`).
#' @return A data.frame with `NA` phenotypes dropped.
#' @export
standardize_pheno <- function(df, sample_id, cohort, phenotype, ...) {
    out <- data.frame(
        sample_id = as.character(df[[sample_id]]),
        cohort    = as.character(df[[cohort]]),
        phenotype = df[[phenotype]],
        stringsAsFactors = FALSE
    )
    extras <- list(...)
    for (nm in names(extras)) out[[nm]] <- df[[extras[[nm]]]]
    out[!is.na(out$phenotype), , drop = FALSE]
}

#' Join a RAV score matrix to a phenotype table
#'
#' Aligns `rav` and `pheno` by `sample_id`, producing one data.frame ready
#' for [fit_rf()] / [fit_logit()] / [cv_oof()].
#'
#' @param rav A samples-x-RAVs matrix (from [project_to_rav()]).
#' @param pheno A standardized phenotype table (from [standardize_pheno()]).
#' @return A data.frame with the phenotype columns followed by RAV columns.
#' @export
join_rav_pheno <- function(rav, pheno) {
    keep <- intersect(rownames(rav), pheno$sample_id)
    if (length(keep) == 0L)
        stop("No overlapping sample IDs between RAV scores and phenotype table.")
    pheno <- pheno[match(keep, pheno$sample_id), , drop = FALSE]
    rav   <- rav[keep, , drop = FALSE]
    cbind(pheno, as.data.frame(rav))
}
