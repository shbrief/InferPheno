# Bundled phenotype heads

This directory holds the serialized `RAVPhenoModel` objects that ship
with the package, one per `availablePheno()` row, plus the compact
cohesion basis (`cohesion_subset.rds`) they project onto.

Each head is stored slimmed: the caret `train` object's training-frame
baggage (the captured formula/`terms` environment, `glm$model`, resampling
folds) is stripped at build time, leaving only what `predict()` needs.
This is what keeps each `.rds` ~0.1 MB instead of tens of MB. The stripping
happens inside `InferPhenoBuilder::build_pheno_head()` (via its internal
`slim_fitted_head()`), which re-checks that predictions are unchanged
before writing.

The shipped heads are fit on the **cohesion-gated de-redundant RAVmodel**
(1,760 RAVs; `rav_prior = "cohesion"`). Each head is a non-penalized
model (binary → `glm`, multiclass → `multinom`, regression → `lm`) over a
top-K-by-|β| RAV panel, refit on all training data. They are produced by
`InferPhenoBuilder::build_all_heads()` (driver:
`InferPhenoBuilder/inst/scripts/build_cohesion_heads.R`).

Panel K per head is the size chosen by the paper's survivor-ranked,
multi-metric selection (a proper scoring rule, gated by calibration slope
∈ [0.85, 1.15] and selection stability, under an a-priori sample-size
ceiling; visualization §8, Table 5). This favours a small, well-calibrated
panel over a larger one that scores higher on raw discrimination but
overfits — most visibly for `cms_subtype` (K=3). `sex` was not evaluated in
§8 (its TCGA data was archived during that run) and falls back to its §4
discrimination-only K=8.

| File            | Phenotype   | Kind       | Panel K | OOF metric                |
|-----------------|-------------|------------|--------:|---------------------------|
| sex.rds         | sex         | binary     |       8 | AUC 0.994                 |
| msi_status.rds  | msi_status  | binary     |      20 | AUC 0.951                 |
| msi_score.rds   | msi_score   | regression |      50 | R² 0.701, ρ 0.441         |
| cms_subtype.rds | cms_subtype | multiclass |       3 | acc 0.714, macro-F1 0.674 |
| neutrophil.rds  | neutrophil  | regression |      10 | R² 0.987, ρ 0.966         |

## Basis loading

`cohesion_subset.rds` is a plain genes-x-RAVs loadings matrix holding only
the union of every head's RAV panel (86 of 1,760 RAVs). Because
`calculateScore()` scores each RAV independently, this column-subset gives
bit-identical scores, and `inferpheno()` resolves it by default — no setup.
It is rebuilt by `../scripts/bundle_rav_subset.R` from the full basis.

The full cohesion basis (~190 MB) is only needed to develop/refit heads
(it spans RAVs the subset omits). Point `inferpheno()` at it to override the
bundled subset:

```r
options(inferpheno.rav_model_path = "/path/to/RAVmodel_cohesion.rds")
# or: Sys.setenv(INFERPHENO_RAV_MODEL = "...")
# or: inferpheno(expr, "sex", rav_model = readRDS("..."))
```

## Legacy heads

[`legacy/`](legacy/) holds the immediately-preceding cohesion heads, backed up
automatically by the builder before the §8 panel-size reconciliation overwrote
them (e.g. `cms_subtype` at the old K=20). Kept for A/B comparison against the
current §8-selected panels.

The earlier heads — fit on the published **C2** basis (4,764 RAVs) by
`benchmark/aim1_*.R` and shipped before the cohesion switch — are kept in
[`legacy_c2/`](legacy_c2/) for A/B comparison. They use `rav_prior = "C2"`,
which `inferpheno()` still resolves automatically via `getModel()`.
