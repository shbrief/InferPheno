#' Infer a clinical phenotype from a gene expression matrix
#'
#' Projects `expr` onto the same RAV basis used to train the chosen
#' `phenotype` head, then returns class probabilities (binary/multiclass)
#' or predicted values (regression), plus a per-sample `confidence`
#' column and a head-level performance summary attached as attributes.
#'
#' @param expr A genes-x-samples numeric matrix or SummarizedExperiment.
#'   Rownames must be HGNC symbols.
#' @param phenotype Character. Phenotype to infer; one of
#'   `availablePheno()$phenotype`.
#' @param rav_model Optional pre-loaded GSS RAVmodel. If `NULL`, the basis
#'   recorded with the phenotype (`rav_prior`) is resolved automatically:
#'   published priors (e.g. `"C2"`) are fetched with
#'   `GenomicSuperSignature::getModel()`, while the custom cohesion basis
#'   resolves to a compact loadings subset bundled with the package (it carries
#'   only the RAVs the supported heads use, so built-in predictions work out of
#'   the box). To use the full basis instead -- e.g. to add or refit heads --
#'   set `getOption("inferpheno.rav_model_path")` or the `INFERPHENO_RAV_MODEL`
#'   environment variable to a full RAVmodel `.rds`; an explicit override always
#'   wins. Pass `rav_model` directly to override either.
#' @param log_input One of `"auto"` (default), `"yes"`, or `"no"`. The RAV basis
#'   is log-native and the heads are trained on log2-scale expression, so linear
#'   input (raw counts / RSEM / TPM) must be log2-transformed before projection.
#'   `"auto"` logs only when the input looks linear (leaving already-log
#'   microarray / log2-CPM input untouched), `"yes"` always applies `log2(x + 1)`,
#'   `"no"` never transforms. Use `"no"` only if you have already log-transformed.
#'
#' @return An `inferpheno_result` (data.frame subclass) with one row per
#'   sample. Always carries a `confidence` column (max class probability
#'   for binary/multiclass; `NA_real_` for regression/cox). The
#'   `phenotype`, `kind`, and `head_perf` attributes report which head
#'   was used and its out-of-fold performance.
#'
#' @seealso [availablePheno()]
#' @examples
#' # Which phenotypes can be inferred today?
#' availablePheno()
#'
#' # Project synthetic expression through the bundled basis and predict sex.
#' # Real use passes a genes-x-samples matrix with HGNC symbol rownames.
#' basis <- readRDS(system.file("extdata", "cohesion_subset.rds",
#'                              package = "InferPheno"))
#' set.seed(1)
#' expr <- matrix(stats::rnorm(nrow(basis) * 4), nrow = nrow(basis),
#'                dimnames = list(rownames(basis), paste0("S", 1:4)))
#' inferpheno(expr, "sex")
#' @export
inferpheno <- function(expr, phenotype, rav_model = NULL,
                       log_input = c("auto", "yes", "no")) {
  log_input <- match.arg(log_input)
  m <- .load_pheno_model(phenotype)
  reg_ph <- .pheno_registry()[[phenotype]]
  if (!is.null(reg_ph) && identical(reg_ph$transfer, "fails")) {
    warning("The '", phenotype, "' head FAILS external transfer: ",
            reg_ph$transfer_note, " In-cohort perf does not carry over to a ",
            "new cohort. See availablePheno().", call. = FALSE)
  }
  if (is.null(rav_model)) {
    rav_model <- .resolve_rav_model(m@rav_prior)
  }
  if (methods::is(expr, "SummarizedExperiment")) {
    expr <- SummarizedExperiment::assay(expr)
  }
  expr <- .to_log_scale(expr, log_input)
  rav <- GenomicSuperSignature::calculateScore(expr, rav_model)
  rownames(rav) <- gsub("\\.", "-", rownames(rav))
  newx <- as.data.frame(rav[, m@rav_cols, drop = FALSE])
  out <- .predict_head(m, newx)

  perf <- .head_perf(m)
  attr(out, "phenotype") <- m@phenotype
  attr(out, "kind")      <- m@kind
  attr(out, "head_perf") <- perf
  class(out) <- c("inferpheno_result", class(out))
  out
}

#' Project an expression matrix onto the RAV basis
#' @param expr A genes-x-samples matrix or SummarizedExperiment.
#' @param prior Character. GSS prior name, resolved with
#'   `GenomicSuperSignature::getModel()`. Default `"C2"`. Ignored when
#'   `rav_model` is supplied.
#' @param rav_model Optional pre-loaded GSS RAV model. When supplied it is used
#'   directly, skipping the `getModel(prior)` fetch -- useful for reusing an
#'   already-loaded basis or for lightweight examples.
#' @param log_input One of `"auto"` (default), `"yes"`, or `"no"` -- see
#'   [inferpheno()]. Linear input is log2-transformed before projection because
#'   the RAV basis is log-native; `"auto"` leaves already-log input untouched.
#' @return A samples-x-RAVs numeric matrix.
#' @examples
#' # Reuse any pre-loaded basis; the mini model ships with
#' # GenomicSuperSignature for lightweight examples.
#' data("miniRAVmodel", package = "GenomicSuperSignature")
#' genes <- rownames(GenomicSuperSignature::RAVindex(miniRAVmodel))
#' set.seed(1)
#' expr <- matrix(stats::rnorm(length(genes) * 3), nrow = length(genes),
#'                dimnames = list(genes, paste0("S", 1:3)))
#' projectExpression(expr, rav_model = miniRAVmodel)
#' @export
projectExpression <- function(expr, prior = "C2", rav_model = NULL,
                              log_input = c("auto", "yes", "no")) {
  log_input <- match.arg(log_input)
  if (is.null(rav_model)) {
    rav_model <- GenomicSuperSignature::getModel(prior, load = TRUE)
  }
  if (methods::is(expr, "SummarizedExperiment")) {
    expr <- SummarizedExperiment::assay(expr)
  }
  expr <- .to_log_scale(expr, log_input)
  out <- GenomicSuperSignature::calculateScore(expr, rav_model)
  rownames(out) <- gsub("\\.", "-", rownames(out))
  out
}

#' Put an expression matrix on log2 scale, only if it is not already
#'
#' The RAV basis is log-native and `calculateScore()` does no internal log, so
#' linear input (raw counts / RSEM / TPM) must be log2-transformed before
#' projection while already-log input (microarray, log2-CPM) must not be, to
#' avoid double-logging. `"auto"` decides by scale: data with any negative value
#' is treated as already-transformed (centred/z-scored), otherwise a 99th
#' percentile above `thresh` marks a linear scale (log2 expression essentially
#' never exceeds ~20).
#'
#' @param m A genes-x-samples numeric matrix.
#' @param log_input `"auto"`, `"yes"`, or `"no"`.
#' @param thresh Upper-quantile cutoff separating linear from log scale.
#' @return `m`, log2-transformed per `log_input`.
#' @keywords internal
#' @noRd
.to_log_scale <- function(m, log_input = "auto", thresh = 30) {
  if (identical(log_input, "no")) return(m)
  if (identical(log_input, "yes")) return(log2(m + 1))
  if (min(m, na.rm = TRUE) < 0) return(m)                 # already centred/log
  if (stats::quantile(m, 0.99, na.rm = TRUE) > thresh) log2(m + 1) else m
}

#' List ready-to-infer phenotypes
#'
#' Reads each bundled head from `inst/extdata/` to surface its training
#' metadata. The list is the contract: every row corresponds to a
#' phenotype that `inferpheno()` can be called on today.
#'
#' @return A data.frame with columns `phenotype`, `kind`, `n_train`,
#'   `metric` (which out-of-fold metric `perf` reports), `perf`, and the
#'   leave-one-SOURCE-out external-validation verdict: `transfer`
#'   (`"transfers"` / `"fails"` / `"untested"`), `transfer_perf` (the external
#'   metric value), and `transfer_note`. A head can score high on `perf`
#'   (in-cohort out-of-fold) yet `fail` external transfer when its training
#'   cohorts confound the label — see `sex`.
#' @seealso [inferpheno()]
#' @examples
#' availablePheno()
#' @export
availablePheno <- function() {
  reg <- .pheno_registry()
  rows <- lapply(names(reg), function(ph) {
    rds <- system.file("extdata", reg[[ph]]$rds, package = "InferPheno")
    if (!nzchar(rds)) {
      stop("Bundled head file missing for ", ph,
           ". Re-install InferPheno.")
    }
    m <- readRDS(rds)
    perf <- .head_perf(m)
    data.frame(phenotype     = ph,
               kind          = m@kind,
               n_train       = as.integer(m@training_meta$n_train %||% NA),
               metric        = perf$metric %||% NA_character_,
               perf          = perf$value %||% NA_real_,
               transfer      = reg[[ph]]$transfer %||% NA_character_,
               transfer_perf = reg[[ph]]$transfer_perf %||% NA_real_,
               transfer_note = reg[[ph]]$transfer_note %||% NA_character_,
               stringsAsFactors = FALSE)
  })
  do.call(rbind, rows)
}

#' RAV columns used by a bundled head
#' @param phenotype Character. One of `availablePheno()$phenotype`.
#' @return A character vector of the RAV column names the phenotype's head
#'   uses as features.
#' @examples
#' listRAVs("msi_status")
#' @export
listRAVs <- function(phenotype) {
  .load_pheno_model(phenotype)@rav_cols
}

#' @export
print.inferpheno_result <- function(x, ...) {
  # `[.data.frame` keeps this subclass but drops the attributes that carry the
  # phenotype/kind/performance metadata. A subset is no longer a self-describing
  # result, so print it as the plain data.frame it has become.
  if (is.null(attr(x, "phenotype"))) {
    return(print(structure(x, class = "data.frame"), ...))
  }
  perf <- attr(x, "head_perf")
  cat(sprintf("inferpheno() result for phenotype '%s' (kind: %s)\n",
              attr(x, "phenotype"), attr(x, "kind")))
  if (!is.null(perf) && !is.null(perf$value) && !is.na(perf$value)) {
    cat(sprintf("Head performance: %s = %.3f (n_train = %s)\n",
                perf$metric, perf$value,
                if (is.null(perf$n_train)) "NA" else perf$n_train))
  } else {
    cat("Head performance: unavailable\n")
  }
  cat("\n")
  NextMethod()
  invisible(x)
}

.predict_head <- function(m, newx) {
  if (m@kind == "binary") {
    pr <- stats::predict(m@model, newdata = newx, type = "prob")
    data.frame(class = stats::predict(m@model, newdata = newx),
               prob_positive = pr[, 2],
               prob_negative = pr[, 1],
               confidence    = pmax(pr[, 1], pr[, 2]))
  } else if (m@kind == "multiclass") {
    pr <- stats::predict(m@model, newdata = newx, type = "prob")
    cbind(data.frame(class = stats::predict(m@model, newdata = newx)),
          as.data.frame(pr),
          data.frame(confidence = apply(pr, 1, max)))
  } else if (m@kind == "regression") {
    data.frame(prediction = as.numeric(stats::predict(m@model, newdata = newx)),
               confidence = NA_real_)
  } else if (m@kind == "cox") {
    data.frame(linear_predictor = as.numeric(
                 stats::predict(m@model, newdata = newx, type = "lp")),
               confidence = NA_real_)
  } else stop("Unknown head kind: ", m@kind)
}

.head_perf <- function(m) {
  meta <- m@training_meta
  metric <- if (!is.null(meta$auc_oof) && !is.na(meta$auc_oof)) {
    "auc_oof"
  } else if (!is.null(meta$r2_oof) && !is.na(meta$r2_oof)) {
    "r2_oof"
  } else if (!is.null(meta$accuracy_oof) && !is.na(meta$accuracy_oof)) {
    "accuracy_oof"
  } else if (!is.null(meta$c_index_oof) && !is.na(meta$c_index_oof)) {
    "c_index_oof"
  } else NA_character_
  value <- if (!is.na(metric)) meta[[metric]] else NA_real_
  list(metric = metric, value = value,
       n_train = meta$n_train %||% NA_integer_)
}

# Published GSS priors resolvable through getModel(); anything else is treated
# as a custom basis that must be supplied locally (it is too large to bundle).
.known_priors <- function() c("C2", "PLIERpriors")

.resolve_rav_model <- function(prior) {
  if (prior %in% .known_priors()) {
    return(GenomicSuperSignature::getModel(prior, load = TRUE))
  }
  # An explicit override always wins: it lets developers point at the full
  # basis (e.g. to add or refit heads). Honour it before the bundled subset.
  path <- getOption("inferpheno.rav_model_path",
                    Sys.getenv("INFERPHENO_RAV_MODEL", unset = ""))
  if (nzchar(path)) {
    if (!file.exists(path)) {
      stop("inferpheno.rav_model_path / INFERPHENO_RAV_MODEL points to a ",
           "missing file: ", path, call. = FALSE)
    }
    return(readRDS(path))
  }
  # Otherwise fall back to the compact loadings subset bundled with the
  # package. `calculateScore()` accepts a plain genes-x-RAVs matrix, and the
  # subset carries every RAV used by the supported heads, so this serves all
  # built-in predictions without the full basis.
  subset_rds <- system.file("extdata", paste0(prior, "_subset.rds"),
                            package = "InferPheno")
  if (nzchar(subset_rds)) {
    return(readRDS(subset_rds))
  }
  stop("Phenotype head uses a custom RAV basis ('", prior, "') that is not a ",
       "GenomicSuperSignature prior, and no bundled subset ('", prior,
       "_subset.rds') was found.\n",
       "Point InferPheno at the full basis .rds with one of:\n",
       "  options(inferpheno.rav_model_path = '/path/to/", prior, ".rds')\n",
       "  Sys.setenv(INFERPHENO_RAV_MODEL = '/path/to/", prior, ".rds')\n",
       "or pass a pre-loaded basis via inferpheno(..., rav_model = readRDS(...)).",
       call. = FALSE)
}

.load_pheno_model <- function(phenotype) {
  reg <- .pheno_registry()
  if (!phenotype %in% names(reg)) {
    stop("Unknown phenotype: ", phenotype, ". See availablePheno().")
  }
  rds <- system.file("extdata", reg[[phenotype]]$rds, package = "InferPheno")
  if (!nzchar(rds)) stop("Bundled head file missing for ", phenotype, ".")
  readRDS(rds)
}

.pheno_registry <- function() {
  # `transfer` records the leave-one-SOURCE-out external-validation verdict
  # from the standing gate (InferPhenoBuilder external_validation.R). It is
  # deliberately surfaced by availablePheno() so a head that scores well
  # in-cohort but fails on a genuinely external source is not silently
  # advertised as trustworthy. Values: "transfers" / "fails" / "untested".
  # `transfer_perf` is the external metric value (same kind as `perf`).
  list(
    sex         = list(rds = "sex.rds",         kind = "binary",
                       transfer = "transfers", transfer_perf = 0.909,
                       transfer_note = paste0(
                         "external GTEx AUC 0.909 (re-scored 2026-09-12 for the ",
                         "K=5 asymptote_cal panel; the prior K=8 panel scored ",
                         "0.940). NOTE: an earlier build on a ",
                         "raw-linear projection FAILED (AUC 0.540) because ",
                         "stability selection on the mis-scaled basis picked ",
                         "cancer-type-proxy RAVs, not sex RAVs. Rebuilding on the ",
                         "log-native projection (2026-07-23) selected genuine ",
                         "Y-chromosome RAVs (EIF1AY/DDX3Y/KDM5D-loaded) that ",
                         "transfer to normal tissue.")),
    msi_status  = list(rds = "msi_status.rds",  kind = "binary",
                       transfer = "transfers", transfer_perf = 0.824,
                       transfer_note = "external GEO CRC (Affymetrix) AUC 0.824."),
    msi_score   = list(rds = "msi_score.rds",   kind = "regression",
                       transfer = "untested", transfer_perf = NA_real_,
                       transfer_note = "no external MANTIS cohort wired yet."),
    cms_subtype = list(rds = "cms_subtype.rds", kind = "multiclass",
                       transfer = "transfers", transfer_perf = 0.924,
                       transfer_note = paste0(
                         "external GSE37892 CRCSC gold CMS accuracy 0.924 ",
                         "(macro-F1 0.884), re-scored 2026-09-12. The previous ",
                         "K=3 head trained on 1,355 samples scored 0.653 here; ",
                         "the gain comes from both the asymptote_cal resize to ",
                         "K=10 and the TCGA sample_id fix that restored 162 ",
                         "training samples.")),
    neutrophil  = list(rds = "neutrophil.rds",  kind = "regression",
                       transfer = "transfers", transfer_perf = 0.615,
                       transfer_note = paste0(
                         "external SLE-WB clinical CBC Spearman 0.615, at the ",
                         "marker-signature measurement ceiling of 0.634 ",
                         "(re-scored 2026-09-12). The previous K=10 head scored ",
                         "0.610, so the asymptote_cal resize to K=5 costs ",
                         "in-cohort R2 but not external transfer.")),
    sle_neutrophil = list(rds = "sle_neutrophil.rds", kind = "regression",
                       transfer = "untested", transfer_perf = NA_real_,
                       transfer_note = paste0(
                         "no external cohort: SLE-WB is itself the external set ",
                         "for `neutrophil`, so the two cannot validate each other. ",
                         "NOTE: SLE-WB is longitudinal (853 samples, 156 ",
                         "participants) and this head was sized with folds ",
                         "assigned per sample, so its R2 0.662 is mildly optimistic ",
                         "(0.645 under participant-grouped folds)."))
  )
}

`%||%` <- function(a, b) if (is.null(a)) b else a
