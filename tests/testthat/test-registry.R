test_that("availablePheno returns the expected five heads", {
  tbl <- availablePheno()
  expect_setequal(tbl$phenotype,
                  c("sex", "msi_status", "msi_score", "cms_subtype",
                    "neutrophil"))
  expect_true(all(tbl$kind %in% c("binary", "multiclass", "regression", "cox")))
  expect_true(all(c("n_train", "metric", "perf") %in% colnames(tbl)))
})

test_that("inferpheno errors on unknown phenotype", {
  expect_error(inferpheno(matrix(0, 2, 2), "not_a_phenotype"),
               "Unknown phenotype")
})
