# Bundled phenotype heads

This directory holds the serialized `RAVPhenoModel` objects that ship
with the package, one per `availablePheno()` row.

The shipped heads are fit on the **cohesion-gated de-redundant RAVmodel**
(2,377 RAVs; `rav_prior = "cohesion"`). Each head is a non-penalized
model (binary → `glm`, multiclass → `multinom`, regression → `lm`) over a
minimal top-K-by-|β| RAV panel, refit on all training data. They are
produced by `aim4/build_cohesion_heads.R` in the cohesion de-redundance
tree and validated by `aim4` section 6 (out-of-fold prediction outcomes).

Panel K is the minimal sufficient size (smallest top-K within tolerance
of the asymptote), except `cms_subtype`, kept at the peak-anchored K=20
because its recommended K=5 overfits to acc 0.853.

| File            | Phenotype   | Kind       | Panel K | OOF metric                |
|-----------------|-------------|------------|--------:|---------------------------|
| sex.rds         | sex         | binary     |       5 | AUC 0.986                 |
| msi_status.rds  | msi_status  | binary     |      20 | AUC 0.974                 |
| msi_score.rds   | msi_score   | regression |      30 | R² 0.610, ρ 0.404         |
| cms_subtype.rds | cms_subtype | multiclass |      20 | acc 0.925, macro-F1 0.917 |
| neutrophil.rds  | neutrophil  | regression |       3 | R² 0.965, ρ 0.937         |

## Basis loading

The cohesion basis (~240 MB) is too large to bundle, so it is **not** a
`GenomicSuperSignature::getModel()` prior. Point `inferpheno()` at the
basis `.rds` before predicting:

```r
options(inferpheno.rav_model_path = "/path/to/RAVmodel_cohesion.rds")
# or: Sys.setenv(INFERPHENO_RAV_MODEL = "...")
# or: inferpheno(expr, "sex", rav_model = readRDS("..."))
```

## Legacy heads

The previous heads — fit on the published **C2** basis (4,764 RAVs) by
`benchmark/aim1_*.R` and shipped before the cohesion switch — are kept in
[`legacy_c2/`](legacy_c2/) for A/B comparison. They use `rav_prior = "C2"`,
which `inferpheno()` still resolves automatically via `getModel()`.
