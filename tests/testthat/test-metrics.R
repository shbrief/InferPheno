test_that("auc_ci recovers a perfect classifier", {
  truth <- factor(c("neg", "neg", "neg", "pos", "pos", "pos"))
  prob  <- c(0.1, 0.2, 0.3, 0.7, 0.8, 0.9)
  r <- auc_ci(truth, prob)
  expect_equal(r$auc, 1)
  expect_lte(r$lo, 1)
  expect_gte(r$hi, 1)
})

test_that("auc_ci on a useless score is near 0.5", {
  set.seed(1L)
  truth <- factor(rep(c("neg", "pos"), each = 50))
  prob  <- runif(100)
  r <- auc_ci(truth, prob)
  expect_lt(abs(r$auc - 0.5), 0.15)
})

test_that("regression_metrics on perfect prediction has rmse 0 and r2 1", {
  truth <- c(1, 2, 3, 4, 5)
  m <- regression_metrics(truth, truth)
  expect_equal(m$rmse, 0)
  expect_equal(m$mae, 0)
  expect_equal(m$r2, 1)
  expect_equal(m$spearman, 1)
})

test_that("brier_ece bounds are respected", {
  truth <- factor(c("0", "0", "1", "1"))
  prob  <- c(0.2, 0.3, 0.7, 0.8)
  m <- brier_ece(truth, prob, positive = "1")
  expect_gte(m$brier, 0)
  expect_lte(m$brier, 1)
  expect_gte(m$ece, 0)
  expect_lte(m$ece, 1)
})
