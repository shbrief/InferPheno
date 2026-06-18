test_that("proj_path honors INFERPHENO_ROOT", {
  withr::with_envvar(c(INFERPHENO_ROOT = "/tmp/inferpheno_test"), {
    expect_equal(proj_path("a", "b.txt"), "/tmp/inferpheno_test/a/b.txt")
  })
})

test_that("standardize_pheno renames and drops NA phenotypes", {
  df <- data.frame(
    sid = c("s1", "s2", "s3"),
    cancer = c("BRCA", "BRCA", "OV"),
    status = c("msi", NA, "mss"),
    age = c(50, 60, 70),
    stringsAsFactors = FALSE
  )
  out <- standardize_pheno(df, "sid", "cancer", "status", age = "age")
  expect_equal(nrow(out), 2L)
  expect_equal(out$sample_id, c("s1", "s3"))
  expect_equal(out$cohort, c("BRCA", "OV"))
  expect_equal(out$age, c(50, 70))
})

test_that("join_rav_pheno aligns and inner-joins", {
  rav <- matrix(seq_len(6), nrow = 3, ncol = 2,
                dimnames = list(c("s1", "s2", "s3"), c("RAV1", "RAV2")))
  pheno <- data.frame(
    sample_id = c("s2", "s3", "s4"),
    cohort    = c("A", "B", "C"),
    phenotype = c("x", "y", "z"),
    stringsAsFactors = FALSE
  )
  joined <- join_rav_pheno(rav, pheno)
  expect_equal(nrow(joined), 2L)
  expect_equal(joined$sample_id, c("s2", "s3"))
  expect_true(all(c("RAV1", "RAV2") %in% colnames(joined)))
})

test_that("join_rav_pheno errors when there is no overlap", {
  rav <- matrix(1, nrow = 1, ncol = 1, dimnames = list("s1", "RAV1"))
  pheno <- data.frame(sample_id = "s2", cohort = "A", phenotype = "x",
                      stringsAsFactors = FALSE)
  expect_error(join_rav_pheno(rav, pheno), "No overlapping sample IDs")
})
