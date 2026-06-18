## Performance metrics with proper CIs and paired comparisons.

#' Binary AUC with DeLong 95% CI
#'
#' @param truth Vector of true class labels.
#' @param prob Numeric vector of positive-class probabilities.
#' @param positive Optional positive class label. Defaults to the second
#'   level of `as.factor(truth)`.
#' @param conf.level Confidence level for the DeLong CI.
#' @return A list with `auc`, `lo`, `hi`, and the underlying `pROC::roc`
#'   object.
#' @export
auc_ci <- function(truth, prob, positive = NULL, conf.level = 0.95) {
    truth <- as.factor(truth)
    if (is.null(positive)) positive <- levels(truth)[2]
    negative <- setdiff(levels(truth), positive)[1]
    r <- pROC::roc(truth, prob, levels = c(negative, positive),
                   direction = "<", quiet = TRUE)
    ci <- pROC::ci.auc(r, conf.level = conf.level, method = "delong")
    list(auc = as.numeric(pROC::auc(r)),
         lo  = as.numeric(ci[1]),
         hi  = as.numeric(ci[3]),
         roc = r)
}

#' Paired DeLong test for two ROC curves on the same truth vector
#'
#' @param truth Vector of true class labels.
#' @param prob_a Probabilities from method A.
#' @param prob_b Probabilities from method B.
#' @return An `htest` from [pROC::roc.test()].
#' @export
paired_delong <- function(truth, prob_a, prob_b) {
    truth <- as.factor(truth)
    r_a <- pROC::roc(truth, prob_a, quiet = TRUE)
    r_b <- pROC::roc(truth, prob_b, quiet = TRUE)
    pROC::roc.test(r_a, r_b, method = "delong", paired = TRUE)
}

#' Brier score and expected calibration error
#'
#' @param truth Vector of true labels.
#' @param prob Numeric vector of positive-class probabilities.
#' @param positive Optional positive class label.
#' @param bins Number of equal-width probability bins for ECE.
#' @return A list with `brier` and `ece`.
#' @export
brier_ece <- function(truth, prob, positive = NULL, bins = 10L) {
    truth <- as.factor(truth)
    if (is.null(positive)) positive <- levels(truth)[2]
    y <- as.integer(truth == positive)
    brier <- mean((prob - y)^2)
    edges <- seq(0, 1, length.out = bins + 1L)
    bin <- findInterval(prob, edges, all.inside = TRUE)
    ece <- 0
    for (b in unique(bin)) {
        ix <- bin == b
        if (!any(ix)) next
        ece <- ece + (sum(ix) / length(prob)) * abs(mean(prob[ix]) - mean(y[ix]))
    }
    list(brier = brier, ece = ece)
}

#' Multiclass accuracy, macro F1, and per-class AUCs
#'
#' @param truth Factor of true class labels.
#' @param class Factor or character of predicted class labels.
#' @param prob_mat A samples-x-classes matrix or data.frame of class
#'   probabilities.
#' @return A list with `accuracy`, `macro_f1`, `per_class_auc`, and
#'   `confusion`.
#' @export
multiclass_metrics <- function(truth, class, prob_mat) {
    truth <- as.factor(truth)
    class <- factor(class, levels = levels(truth))
    cm <- caret::confusionMatrix(class, truth)
    per_class_auc <- vapply(levels(truth), function(lab) {
        if (!lab %in% colnames(prob_mat)) return(NA_real_)
        pROC::auc(pROC::roc(truth == lab, prob_mat[, lab], quiet = TRUE))
    }, numeric(1))
    list(
        accuracy = unname(cm$overall["Accuracy"]),
        macro_f1 = mean(cm$byClass[, "F1"], na.rm = TRUE),
        per_class_auc = per_class_auc,
        confusion = cm$table
    )
}

#' Regression metrics for the MSI score task
#'
#' @param truth Numeric vector of true values.
#' @param pred Numeric vector of predicted values.
#' @return A list with `rmse`, `mae`, `r2`, and `spearman`.
#' @export
regression_metrics <- function(truth, pred) {
    res <- truth - pred
    list(
        rmse     = sqrt(mean(res^2)),
        mae      = mean(abs(res)),
        r2       = 1 - sum(res^2) / sum((truth - mean(truth))^2),
        spearman = stats::cor(truth, pred, method = "spearman")
    )
}
