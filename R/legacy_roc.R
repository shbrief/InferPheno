## Legacy ROC helpers retained for backward compatibility with the
## existing notebooks. For new code prefer [auc_ci()] (which also returns a
## DeLong CI) and `pROC::ggroc()` on the resulting `roc` object.

#' Manual AUC by trapezoidal integration
#'
#' Kept for backward compatibility with notebooks that pre-date
#' [auc_ci()]. New code should use [auc_ci()].
#'
#' @param pred_probs Numeric vector of predicted positive-class probabilities.
#' @param actual_labels Integer/logical vector of true labels (1 = positive).
#' @return A scalar AUC.
#' @export
calculate_auc <- function(pred_probs, actual_labels) {
    ord <- order(-pred_probs)
    pred_probs <- pred_probs[ord]
    actual_labels <- actual_labels[ord]
    n_pos <- sum(actual_labels == 1)
    n_neg <- sum(actual_labels == 0)
    tpr <- cumsum(actual_labels == 1) / n_pos
    fpr <- cumsum(actual_labels == 0) / n_neg
    sum(diff(fpr) * (tpr[-1] + tpr[-length(tpr)]) / 2)
}

#' Base-R ROC plotter
#'
#' Kept for backward compatibility. New code should call `pROC::ggroc()`
#' on the `roc` object returned by [auc_ci()].
#'
#' @inheritParams calculate_auc
#' @return Invisibly `NULL`; called for the plot side effect.
#' @export
plot_roc <- function(pred_probs, actual_labels) {
    ord <- order(-pred_probs)
    pred_probs <- pred_probs[ord]
    actual_labels <- actual_labels[ord]
    n_pos <- sum(actual_labels == 1)
    n_neg <- sum(actual_labels == 0)
    tpr <- cumsum(actual_labels == 1) / n_pos
    fpr <- cumsum(actual_labels == 0) / n_neg
    auc <- sum(diff(fpr) * (tpr[-1] + tpr[-length(tpr)]) / 2)
    plot(fpr, tpr, type = "l", col = "blue",
         xlab = "False Positive Rate", ylab = "True Positive Rate",
         main = "ROC Curve")
    graphics::abline(0, 1, lty = 2, col = "gray")
    graphics::legend("bottomright",
                     legend = sprintf("AUC = %.3f", auc), bty = "n")
}

#' pROC-based ROC plotter
#'
#' @param actual_labels Vector of true class labels.
#' @param pred_probs Numeric vector of positive-class probabilities.
#' @param plot_title Plot title.
#' @return Invisibly `NULL`; called for the plot side effect.
#' @export
plot_roc_pROC <- function(actual_labels, pred_probs, plot_title = "ROC Curve") {
    roc_obj <- pROC::roc(actual_labels, pred_probs, quiet = TRUE)
    plot(roc_obj, main = plot_title, col = "blue", lwd = 2)
    graphics::abline(0, 1, lty = 2, col = "gray")
    graphics::legend("bottomright",
                     legend = sprintf("AUC = %.3f", pROC::auc(roc_obj)),
                     bty = "n")
}
