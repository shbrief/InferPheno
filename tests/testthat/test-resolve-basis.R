test_that("bundled cohesion subset covers every head's RAV panel", {
  rds <- system.file("extdata", "cohesion_subset.rds", package = "InferPheno")
  skip_if(!nzchar(rds), "cohesion_subset.rds not installed")
  subset_mat <- readRDS(rds)
  expect_true(is.matrix(subset_mat))
  expect_type(subset_mat, "double")

  union_ravs <- unique(unlist(lapply(availablePheno()$phenotype, listRAVs)))
  expect_true(all(union_ravs %in% colnames(subset_mat)))
})

test_that("inferpheno runs from the bundled subset with no override set", {
  rds <- system.file("extdata", "cohesion_subset.rds", package = "InferPheno")
  skip_if(!nzchar(rds), "cohesion_subset.rds not installed")

  # Ensure neither override is in effect, restoring afterwards.
  old_opt <- getOption("inferpheno.rav_model_path")
  old_env <- Sys.getenv("INFERPHENO_RAV_MODEL", unset = NA)
  options(inferpheno.rav_model_path = NULL)
  Sys.unsetenv("INFERPHENO_RAV_MODEL")
  on.exit({
    options(inferpheno.rav_model_path = old_opt)
    if (!is.na(old_env)) Sys.setenv(INFERPHENO_RAV_MODEL = old_env)
  }, add = TRUE)

  subset_mat <- readRDS(rds)
  set.seed(1)
  genes <- rownames(subset_mat)
  expr <- matrix(stats::rnorm(length(genes) * 6), nrow = length(genes),
                 dimnames = list(genes, paste0("S", seq_len(6))))

  res <- inferpheno(expr, "sex")
  expect_s3_class(res, "inferpheno_result")
  expect_equal(nrow(res), 6L)
  expect_true("confidence" %in% colnames(res))
})
