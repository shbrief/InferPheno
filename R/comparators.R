## Wrappers for the supervised comparator tools used in Aim 1. Each returns
## predicted class + probabilities in a uniform shape so they can be
## compared via paired_delong() against RAV-RF.

#' CMSclassifier on a genes-x-samples expression matrix
#'
#' Works around a packaging bug in CMSclassifier 1.0.0 where the
#' `model` dataset name in `data/datalist` is not recognised by recent
#' R, so the internal `data(model)` call inside `listModelGenes("RF")`
#' fails with `object 'mgenes' not found`. We pre-load the exported
#' `mgenes` and `finalModel` objects into a temporary local copy of
#' `listModelGenes` before delegating.
#'
#' @param expr Genes-x-samples matrix with Entrez gene IDs as rownames.
#' @return A list with `class` and `prob`.
#' @export
run_cmsclassifier <- function(expr) {
    if (!requireNamespace("CMSclassifier", quietly = TRUE))
        stop("CMSclassifier not installed.")
    ## Patch the CMSclassifier namespace. Both `listModelGenes("RF")`
    ## and `classifyCMS.RF` start with `data(model)`, which is a no-op
    ## under recent R (datalist entry not recognised), leaving the
    ## downstream calls without `mgenes` / `finalModel`. We replace both
    ## with closures that pull the exported objects directly.
    ns <- asNamespace("CMSclassifier")
    if (!isTRUE(attr(ns$classifyCMS.RF, "rpi_patched"))) {
        patched_listModelGenes <- function(method) {
            if (method == "RF")  return(CMSclassifier::mgenes)
            if (method == "SSP") return(rownames(CMSclassifier::centroids))
            NULL
        }
        patched_classifyCMS.RF <- function(Exp, center = TRUE,
                                          minPosterior = 0.5) {
            Exp <- ns$naImpute(Exp, G = CMSclassifier::mgenes)
            if (center) Exp <- Exp - rowMeans(Exp)
            res <- as.data.frame(stats::predict(
                CMSclassifier::finalModel, t(Exp), type = "prob"))
            names(res) <- paste0(names(res), ".posteriorProb")
            res$predictedCMS <- res$nearestCMS <- apply(res, 1, function(z) {
                y <- which(z == max(z))
                paste(c("CMS1", "CMS2", "CMS3", "CMS4")[y], collapse = ",")
            })
            res$predictedCMS[which(apply(res[, 1:4], 1, max) < minPosterior)] <-
                NA
            names(res) <- paste0("RF.", names(res))
            res
        }
        attr(patched_classifyCMS.RF, "rpi_patched") <- TRUE
        for (sym in c("listModelGenes", "classifyCMS.RF")) {
            unlockBinding(sym, ns)
        }
        assign("listModelGenes", patched_listModelGenes, envir = ns)
        assign("classifyCMS.RF", patched_classifyCMS.RF, envir = ns)
        for (sym in c("listModelGenes", "classifyCMS.RF")) {
            lockBinding(sym, ns)
        }
    }
    ## classifyCMS.RF -> naImpute uses `names(Exp)` and rbind, which
    ## requires Exp to be a data.frame, not a matrix.
    expr_df <- as.data.frame(expr, stringsAsFactors = FALSE,
                             check.names = FALSE)
    rf <- CMSclassifier::classifyCMS(expr_df, method = "RF")
    cls  <- stats::setNames(rf$predictedCMS$RF, colnames(expr_df))
    ## classifyCMS scatters per-class posteriors into 4 columns; rebuild
    ## a samples-x-classes matrix.
    det  <- rf$RF.details
    cols <- paste0("RF.CMS", 1:4, ".posteriorProb")
    prob <- as.matrix(det[, cols, drop = FALSE])
    colnames(prob) <- c("CMS1", "CMS2", "CMS3", "CMS4")
    rownames(prob) <- rownames(det)
    list(class = cls, prob = prob)
}

#' CMScaller on a genes-x-samples expression matrix
#'
#' @param expr Genes-x-samples matrix.
#' @param RNAseq Logical; pass-through to [CMScaller::CMScaller()].
#' @return A list with `class` and `prob`.
#' @export
run_cmscaller <- function(expr, RNAseq = TRUE) {
    if (!requireNamespace("CMScaller", quietly = TRUE))
        stop("CMScaller not installed.")
    res <- CMScaller::CMScaller(emat = expr, RNAseq = RNAseq, doPlot = FALSE)
    prob <- res[, c("CMS1", "CMS2", "CMS3", "CMS4")]
    list(class = as.character(res$prediction), prob = prob)
}

#' XIST / Y-chromosome panel sex classifier
#'
#' @param expr Genes-x-samples log-expression matrix with HGNC symbols as
#'   rownames.
#' @return A list with `class`, `prob`, and the underlying `score`.
#' @export
run_xist_y_panel <- function(expr) {
    y_panel <- c("RPS4Y1", "EIF1AY", "DDX3Y", "USP9Y", "UTY", "NLGN4Y", "KDM5D")
    have <- intersect(c("XIST", y_panel), rownames(expr))
    if (length(have) < 2L) stop("Not enough panel genes present.")
    z <- t(scale(t(expr[have, , drop = FALSE])))
    xist <- if ("XIST" %in% have) z["XIST", ] else rep(0, ncol(expr))
    y_genes <- setdiff(have, "XIST")
    y_mean <- if (length(y_genes) > 0L)
        colMeans(z[y_genes, , drop = FALSE])
    else rep(0, ncol(expr))
    score <- y_mean - xist
    if (any(is.na(score))) score[is.na(score)] <- 0
    prob_df <- data.frame(female = stats::pnorm(-score),
                          male   = stats::pnorm( score),
                          stringsAsFactors = FALSE)
    rownames(prob_df) <- colnames(expr)
    cls <- stats::setNames(ifelse(score > 0, "male", "female"), colnames(expr))
    list(class = cls, prob = prob_df, score = score)
}

#' MMR-gene expression-panel MSI scorer
#'
#' Rough expression-based surrogate for MSIsensor-pro on DNA. High-MSI
#' tumors typically have lower MLH1/MSH2 expression (e.g. MLH1 promoter
#' methylation), so we negate the mean z-score to get an MSI-up score.
#'
#' @param expr Genes-x-samples expression matrix with HGNC symbols as
#'   rownames.
#' @return A list with `class`, `prob`, and the underlying `score`.
#' @export
run_msi_expression_panel <- function(expr) {
    sig <- c("MLH1", "MSH2", "MSH6", "PMS2", "EPCAM",
             "MLH3", "PMS1", "MUTYH")
    have <- intersect(sig, rownames(expr))
    if (length(have) == 0L) stop("No MSI panel genes present.")
    z <- t(scale(t(expr[have, , drop = FALSE])))
    score <- -colMeans(z, na.rm = TRUE)
    if (any(is.na(score))) score[is.na(score)] <- 0
    prob_df <- data.frame(mss = stats::pnorm(-score),
                          msi = stats::pnorm( score),
                          stringsAsFactors = FALSE)
    rownames(prob_df) <- colnames(expr)
    cls <- stats::setNames(ifelse(score > stats::median(score), "msi", "mss"),
                    colnames(expr))
    list(class = cls, prob = prob_df, score = score)
}
