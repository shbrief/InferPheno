## RAV x phenotype atlas builder (Aim 3).

#' Score one RAV against one phenotype
#'
#' Uses a univariate test for screening plus an AUC summary. For binary
#' phenotypes, returns the DeLong AUC + 95% CI and the Wilcoxon p-value.
#' For multiclass phenotypes, returns eta-squared and the Kruskal-Wallis
#' p-value.
#'
#' @param rav_score Numeric vector of RAV scores for one RAV across samples.
#' @param pheno Factor or character of phenotype labels (same length as
#'   `rav_score`).
#' @param k Reserved for future stratified-CV extension.
#' @param seed RNG seed.
#' @return A one-row data.frame or `NULL` if the phenotype has fewer than 2
#'   levels.
#' @export
score_rav_phenotype <- function(rav_score, pheno, k = 5L, seed = 1L) {
    pheno <- factor(pheno)
    if (nlevels(pheno) < 2L) return(NULL)

    if (nlevels(pheno) == 2L) {
        a <- auc_ci(pheno, rav_score)
        wp <- suppressWarnings(stats::wilcox.test(rav_score ~ pheno)$p.value)
        data.frame(auc = a$auc, auc_lo = a$lo, auc_hi = a$hi,
                   wilcox_p = wp, kind = "binary")
    } else {
        kw <- suppressWarnings(stats::kruskal.test(rav_score ~ pheno)$p.value)
        eta2 <- {
            ss_b <- sum(tapply(rav_score, pheno,
                               function(z) length(z) * (mean(z) - mean(rav_score))^2))
            ss_b / sum((rav_score - mean(rav_score))^2)
        }
        data.frame(eta2 = eta2, kw_p = kw, kind = "multiclass")
    }
}

#' Fast rank-based AUC for one binary phenotype against many RAVs
#'
#' Mann-Whitney-U identity, vectorized over RAVs. Takes pre-computed
#' column ranks so it can be called cheaply inside permutation loops.
#' @keywords internal
.abs_auc_vec <- function(y_pos_lgl, ranks) {
    n_pos <- sum(y_pos_lgl)
    n_neg <- length(y_pos_lgl) - n_pos
    if (n_pos == 0L || n_neg == 0L) return(rep(NA_real_, ncol(ranks)))
    r_pos <- colSums(ranks[y_pos_lgl, , drop = FALSE])
    auc <- (r_pos - n_pos * (n_pos + 1) / 2) / (n_pos * n_neg)
    pmax(auc, 1 - auc) - 0.5
}

#' Batched abs-AUC permutation null over many label shufflings
#'
#' Builds an N x B indicator matrix of permuted positive-class labels,
#' then `crossprod(ranks, B_mat)` produces a K x B matrix of positive-
#' rank sums in a single BLAS call. Much faster than per-permutation
#' looping when both K (RAVs) and N (samples) are large.
#' @keywords internal
.abs_auc_perm_batch <- function(y_pos_lgl, ranks, n_perm) {
    n <- length(y_pos_lgl)
    n_pos <- sum(y_pos_lgl)
    n_neg <- n - n_pos
    if (n_pos == 0L || n_neg == 0L)
        return(matrix(NA_real_, ncol(ranks), n_perm))
    B <- matrix(0, n, n_perm)
    for (b in seq_len(n_perm)) B[sample.int(n, n_pos), b] <- 1
    r_pos <- crossprod(ranks, B)
    auc <- (r_pos - n_pos * (n_pos + 1) / 2) / (n_pos * n_neg)
    pmax(auc, 1 - auc) - 0.5
}

#' Vectorized eta-squared for one multiclass phenotype across many RAVs
#' @keywords internal
.eta2_vec <- function(y_fac, R) {
    grand <- colMeans(R)
    ss_tot <- colSums(sweep(R, 2L, grand, "-")^2)
    grp_means <- vapply(levels(y_fac), function(g)
        colMeans(R[y_fac == g, , drop = FALSE]),
        numeric(ncol(R)))
    grp_n <- as.integer(table(y_fac)[levels(y_fac)])
    diff2 <- sweep(grp_means, 1L, grand, "-")^2
    ss_b <- as.vector(diff2 %*% grp_n)
    ss_b / ss_tot
}

#' Build the full RAV x phenotype atlas
#'
#' Vectorized scoring with phenotype-kind-aware permutation nulls:
#' binary phenotypes use abs(AUC-0.5) as the effect, multiclass use
#' eta^2, and ordinal use abs(Spearman rho). Each kind permutes against
#' its own null.
#'
#' @param rav Samples-x-RAVs matrix.
#' @param pheno_list A named list of vectors aligned to `rownames(rav)`.
#'   Class is inferred (logical/2-level factor -> binary; ordered factor
#'   or numeric -> ordinal; >2 unordered levels -> multiclass).
#' @param n_perm Number of label-shuffling permutations for the FDR null.
#' @param seed RNG seed.
#' @return A tidy data.frame with one row per (phenotype, RAV) pair and
#'   columns `effect` (kind-appropriate), `kind`, `n`, `perm_p`, `fdr`.
#' @export
build_atlas <- function(rav, pheno_list, n_perm = 200L, seed = 1L,
                        verbose = TRUE) {
    set.seed(seed)
    rav_names <- colnames(rav)
    out_rows <- list()
    for (pn in names(pheno_list)) {
        y_raw <- pheno_list[[pn]]
        keep <- !is.na(y_raw)
        y <- y_raw[keep]; R <- rav[keep, , drop = FALSE]
        n <- length(y)
        if (verbose) message(sprintf(
            "[atlas] %s: n=%d, scoring %d RAVs x %d perms",
            pn, n, ncol(R), n_perm))
        ranks <- apply(R, 2L, rank, ties.method = "average")
        if (is.ordered(y) || is.numeric(y)) {
            kind <- "ordinal"
            y_rank <- rank(as.numeric(y), ties.method = "average")
            yc <- y_rank - mean(y_rank)
            denom_y <- sqrt(sum(yc^2))
            rc <- sweep(ranks, 2L, colMeans(ranks), "-")
            denom_r <- sqrt(colSums(rc^2))
            obs <- abs(crossprod(rc, yc)[, 1L] / (denom_r * denom_y))
            B <- matrix(0, n, n_perm)
            for (b in seq_len(n_perm)) B[, b] <- sample(yc)
            null_max <- abs(crossprod(rc, B) /
                            (denom_r %o% sqrt(colSums(B^2))))
        } else {
            y_fac <- factor(y)
            if (nlevels(y_fac) < 2L) next
            if (nlevels(y_fac) == 2L) {
                kind <- "binary"
                pos <- levels(y_fac)[2L]
                obs <- .abs_auc_vec(y == pos, ranks)
                null_max <- .abs_auc_perm_batch(y == pos, ranks, n_perm)
            } else {
                kind <- "multiclass"
                obs <- .eta2_vec(y_fac, R)
                null_max <- vapply(seq_len(n_perm), function(b)
                    .eta2_vec(y_fac[sample.int(n)], R),
                    numeric(ncol(R)))
            }
        }
        perm_p <- (rowSums(null_max >= obs, na.rm = TRUE) + 1) /
                  (ncol(null_max) + 1)
        out_rows[[pn]] <- data.frame(
            phenotype = pn, rav = rav_names, kind = kind, n = n,
            effect = unname(obs), perm_p = perm_p,
            stringsAsFactors = FALSE
        )
    }
    out <- do.call(rbind, out_rows)
    out$fdr <- stats::p.adjust(out$perm_p, method = "BH")
    out
}
