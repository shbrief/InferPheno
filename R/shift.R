## Cross-cohort / cross-platform / cross-tissue shift harness (Aim 2).

#' Run a train-cohorts -> test-cohorts shift experiment
#'
#' @param dat A data.frame with columns `sample_id`, `cohort`, `phenotype`,
#'   plus the RAV columns named in `rav_cols`.
#' @param rav_cols Character. RAV columns to use as features.
#' @param train_cohorts Character. Values of `dat$cohort` to use for training.
#' @param test_cohorts Character. Values of `dat$cohort` to use for evaluation.
#' @param fitter A fitter function such as [fit_rf()].
#' @param seed RNG seed.
#' @param ... Extra arguments forwarded to `fitter`.
#' @return A list with the fitted `model` and a per-cohort metric
#'   `by_cohort` data.frame.
#' @export
run_shift <- function(dat, rav_cols, train_cohorts, test_cohorts,
                      fitter = fit_rf, seed = 1L, ...) {
    train_df <- dat[dat$cohort %in% train_cohorts, , drop = FALSE]
    test_df  <- dat[dat$cohort %in% test_cohorts,  , drop = FALSE]
    stopifnot(nrow(train_df) > 0, nrow(test_df) > 0)
    model <- fitter(train_df, rav_cols, seed = seed, ...)

    rows <- lapply(split(test_df, test_df$cohort), function(sub) {
        pr <- predict_inferpheno(model, sub)
        if (length(levels(factor(sub$phenotype))) == 2L) {
            pos <- levels(factor(sub$phenotype))[2]
            pos_prob <- if (is.data.frame(pr$prob) || is.matrix(pr$prob))
                pr$prob[, pos] else pr$prob
            a <- auc_ci(sub$phenotype, pos_prob, positive = pos)
            data.frame(cohort = unique(sub$cohort), n = nrow(sub),
                       auc = a$auc, lo = a$lo, hi = a$hi)
        } else {
            mm <- multiclass_metrics(sub$phenotype, pr$class,
                                     as.data.frame(pr$prob))
            data.frame(cohort = unique(sub$cohort), n = nrow(sub),
                       accuracy = mm$accuracy, macro_f1 = mm$macro_f1)
        }
    })
    list(model = model, by_cohort = do.call(rbind, rows))
}

#' Compute a shift-severity score for a (train, test) cohort pair
#'
#' Top-PC mean distance on RAV space. Useful as the x-axis of Fig 3.
#'
#' @inheritParams run_shift
#' @param npc Number of principal components to summarise distance over.
#' @return A non-negative numeric scalar.
#' @export
shift_severity <- function(dat, rav_cols, train_cohorts, test_cohorts,
                           npc = 3L) {
    mat <- as.matrix(dat[, rav_cols, drop = FALSE])
    pcs <- stats::prcomp(mat, scale. = TRUE)$x[, seq_len(npc), drop = FALSE]
    tr  <- colMeans(pcs[dat$cohort %in% train_cohorts, , drop = FALSE])
    te  <- colMeans(pcs[dat$cohort %in% test_cohorts,  , drop = FALSE])
    sqrt(sum((tr - te)^2))
}
