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
top-K-by-|β| RAV panel, refit on all training data. As of 2026-09-12 all six shipped heads are sized by a single procedure,
`asymptote_cal`, at one protocol and seed. The build driver, per-head selection
tables and the validation report for that run live in
`InferPhenoPaper/paper/asymptote_cal_final/`. Earlier heads are preserved under
`legacy_2026-09-12/`. Heads are produced by
`InferPhenoBuilder::build_pheno_head()` / `build_all_heads()`.

`cohesion_subset.rds` must be regenerated whenever a shipped panel changes: it
carries exactly the union of the shipped heads' `@rav_cols` (currently 97 RAVs),
so `inferpheno()` cannot resolve a basis for a RAV it omits. Rebuild it with
`inst/scripts/bundle_rav_subset.R`.

Panel K per head is the size chosen by the paper's survivor-ranked,
multi-metric selection (a proper scoring rule, gated by calibration slope
∈ [0.85, 1.15] and selection stability, under an a-priori sample-size
ceiling; visualization §8, Table 5). This favours a small, well-calibrated
panel over a larger one that scores higher on raw discrimination but
overfits — most visibly for `cms_subtype` (K=3). `sex` was not evaluated in
§8 (its TCGA data was archived during that run) and falls back to its §4
discrimination-only K=8.

| File            | Phenotype   | Kind       | Panel K | OOF metric                | External transfer            |
|-----------------|-------------|------------|--------:|---------------------------|------------------------------|
| sex.rds         | sex         | binary     |       8 | AUC 0.994                 | **FAILS** — AUC 0.540 (GTEx) |
| msi_status.rds  | msi_status  | binary     |      20 | AUC 0.951                 | transfers — AUC 0.824 (GEO)  |
| msi_score.rds   | msi_score   | regression |      50 | R² 0.701, ρ 0.441         | untested (no MANTIS cohort)  |
| cms_subtype.rds | cms_subtype | multiclass |       3 | acc 0.714, macro-F1 0.674 | transfers — acc 0.653 (gold) |
| neutrophil.rds  | neutrophil  | regression |      10 | R² 0.987, ρ 0.966         | transfers — ρ 0.610 (SLE-WB) |

> **The `sex` head does not transfer across cohorts.** Its high OOF AUC is a
> training-cohort confound — sex was nearly identical to cancer type in the TCGA
> training cohorts (OV 100% female, PRAD 100% male), so its panel RAVs encode
> tumor biology, not sex, and collapse to near chance (AUC 0.540) on GTEx normal
> tissue. `inferpheno(expr, "sex")` warns at call time; do not use it for
> cross-cohort inference pending a de-confounded retrain. The `transfer` column
> of `availablePheno()` carries the same verdict for every head, produced by the
> standing external-validation gate
> (`InferPhenoBuilder/inst/scripts/external_validation/`).

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
