## Unified fitters. Every aim calls these — never re-implements them.

#' Fit a random forest on RAV columns
#'
#' @param df A data.frame containing `phenotype` and the columns in `rav_cols`.
#' @param rav_cols Character. RAV column names to use as features.
#' @param ntree Number of trees.
#' @param mtry Number of variables sampled per split. Defaults to
#'   `floor(sqrt(length(rav_cols)))`.
#' @param seed RNG seed.
#' @return A `randomForest` object.
#' @export
fit_rf <- function(df, rav_cols, ntree = 1000, mtry = NULL, seed = 1L) {
    stopifnot(all(rav_cols %in% colnames(df)))
    set.seed(seed)
    x <- df[, rav_cols, drop = FALSE]
    y <- df$phenotype
    if (is.character(y)) y <- factor(y)
    if (is.null(mtry)) mtry <- max(1L, floor(sqrt(length(rav_cols))))
    randomForest::randomForest(x = x, y = y, ntree = ntree, mtry = mtry,
                               importance = TRUE)
}

#' Fit a penalized logistic regression
#'
#' Returns a `cv.glmnet` fit at `lambda.min` wrapped in an
#' `inferpheno_logit` object so [predict_inferpheno()] can dispatch on it.
#'
#' @inheritParams fit_rf
#' @param alpha Elastic-net mixing parameter; 1 = lasso, 0 = ridge.
#' @return An `inferpheno_logit` object.
#' @export
fit_logit <- function(df, rav_cols, alpha = 1, seed = 1L) {
    set.seed(seed)
    x <- as.matrix(df[, rav_cols, drop = FALSE])
    y <- df$phenotype
    if (is.character(y)) y <- factor(y)
    cv <- glmnet::cv.glmnet(x, y, family = "binomial", alpha = alpha)
    structure(list(fit = cv, rav_cols = rav_cols, alpha = alpha),
              class = "inferpheno_logit")
}

#' Fit a Cox proportional hazards model with RAV scores as covariates
#'
#' @inheritParams fit_rf
#' @param time_col Name of the survival-time column.
#' @param event_col Name of the event/censoring column.
#' @return A `coxph` object.
#' @export
fit_cox <- function(df, rav_cols, time_col = "time", event_col = "event") {
    stopifnot(all(c(time_col, event_col) %in% colnames(df)))
    form <- stats::as.formula(sprintf("survival::Surv(%s, %s) ~ %s",
                                      time_col, event_col,
                                      paste(rav_cols, collapse = " + ")))
    survival::coxph(form, data = df)
}

#' Predict probability/class on a new (RAV-projected) data.frame
#'
#' Dispatches on the class of `model` so callers don't need to know whether
#' it came from [fit_rf()] or [fit_logit()].
#'
#' @param model A fitted model from [fit_rf()] or [fit_logit()].
#' @param newdf A data.frame containing the RAV columns the model was
#'   trained on.
#' @return A list with elements `prob` and `class`.
#' @export
predict_inferpheno <- function(model, newdf) {
    rav_cols <- if (inherits(model, "inferpheno_logit")) model$rav_cols
                else rownames(model$importance)
    x <- newdf[, intersect(colnames(newdf), rav_cols), drop = FALSE]
    if (inherits(model, "inferpheno_logit")) {
        p2 <- as.numeric(stats::predict(model$fit, newx = as.matrix(x),
                                        s = "lambda.min", type = "response"))
        cls_levels <- model$fit$glmnet.fit$classnames
        prob_df <- stats::setNames(
            data.frame(1 - p2, p2),
            cls_levels)
        list(prob = prob_df,
             class = ifelse(p2 > 0.5, cls_levels[2], cls_levels[1]))
    } else {
        prob <- as.data.frame(stats::predict(model, newdata = x, type = "prob"))
        cls  <- as.character(stats::predict(model, newdata = x, type = "response"))
        list(prob = prob, class = cls)
    }
}

#' Stratified k-fold cross-validation
#'
#' Produces out-of-fold probabilities and class labels in a tidy data.frame.
#'
#' @inheritParams fit_rf
#' @param k Number of folds.
#' @param fitter A fitter function such as [fit_rf()] or [fit_logit()].
#' @param ... Extra arguments forwarded to `fitter`.
#' @return A data.frame with columns `sample_id, fold, truth, class,
#'   prob_<level>...`.
#' @export
cv_oof <- function(df, rav_cols, k = 5, fitter = fit_rf, seed = 1L, ...) {
    set.seed(seed)
    folds <- caret::createFolds(df$phenotype, k = k, list = TRUE,
                                returnTrain = FALSE)
    out <- vector("list", k)
    for (i in seq_len(k)) {
        test_ix <- folds[[i]]
        train  <- df[-test_ix, , drop = FALSE]
        test   <- df[ test_ix, , drop = FALSE]
        m      <- fitter(train, rav_cols, ...)
        pr     <- predict_inferpheno(m, test)
        prob_mat <- if (is.matrix(pr$prob) || is.data.frame(pr$prob))
            as.data.frame(pr$prob)
        else data.frame(positive = pr$prob)
        out[[i]] <- cbind(
            data.frame(sample_id = test$sample_id, fold = i,
                       truth = test$phenotype, class = pr$class,
                       stringsAsFactors = FALSE),
            prob_mat
        )
    }
    do.call(rbind, out)
}
