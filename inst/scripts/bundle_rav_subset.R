## Build the compact cohesion basis bundled with the package.
##
## calculateScore() reads only the genes-x-RAVs loadings and scores each RAV
## independently, so the supported heads need only the union of their RAV
## panels (~72 of 1,760 RAVs). We ship that column-subset as a plain matrix;
## calculateScore() accepts a matrix directly and returns bit-identical scores
## for the retained RAVs. This lets inferpheno() run without the ~190 MB basis.
##
## Run from the InferPheno package root, with the full cohesion RAVmodel
## reachable via the same override inferpheno() honours:
##   INFERPHENO_RAV_MODEL=/path/to/RAVmodel_cohesion.rds \
##     Rscript inst/scripts/bundle_rav_subset.R

suppressPackageStartupMessages({
  library(methods)
  library(InferPheno)
  library(GenomicSuperSignature)
})

extdata_dir <- file.path("inst", "extdata")

## --- RAV union across the bundled heads (authoritative @rav_cols) ----------
heads <- c("sex", "msi_status", "msi_score", "cms_subtype", "neutrophil",
           "sle_neutrophil")
union_ravs <- unique(unlist(lapply(heads, function(ph) {
  readRDS(file.path(extdata_dir, paste0(ph, ".rds")))@rav_cols
})))
message("RAV union across heads: ", length(union_ravs), " RAVs")

## --- Full cohesion basis (override path) -----------------------------------
path <- getOption("inferpheno.rav_model_path",
                  Sys.getenv("INFERPHENO_RAV_MODEL", unset = ""))
if (!nzchar(path) || !file.exists(path)) {
  stop("Set INFERPHENO_RAV_MODEL (or option inferpheno.rav_model_path) to the ",
       "full cohesion RAVmodel .rds before running this script.")
}
full <- readRDS(path)
loadings <- RAVindex(full)   # genes x 1,760 matrix (== assay(full))

missing <- setdiff(union_ravs, colnames(loadings))
if (length(missing)) {
  stop("Full basis is missing RAVs the heads require: ",
       paste(missing, collapse = ", "))
}

subset_mat <- loadings[, union_ravs, drop = FALSE]
storage.mode(subset_mat) <- "double"

out <- file.path(extdata_dir, "cohesion_subset.rds")
saveRDS(subset_mat, out, compress = "xz")
message(sprintf("Wrote %s: %d genes x %d RAVs, %.1f MB",
                out, nrow(subset_mat), ncol(subset_mat),
                file.info(out)$size / 1e6))
