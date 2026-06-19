# InferPheno

Infer clinical and molecular phenotypes from gene expression data using a
fixed, unsupervised Replicable Axes of Variation (RAV) basis from
[GenomicSuperSignature](https://github.com/shbrief/GenomicSuperSignature).
The package ships pre-fit supervised *heads* for sex, microsatellite
instability (MSI) status, MSI score, consensus molecular subtype (CMS), and
neutrophil infiltration. It is designed as a drop-in alternative to
task-specific classifiers when interpretability and cross-platform transfer
matter more than peak accuracy on the source platform.

## How it works

Inference is a two-stage pipeline:

1. **Project** — new expression data is scored against a *fixed*,
   pre-computed RAV basis. No fitting happens here; the basis is the same
   one used to train every head, which is what makes predictions
   comparable across platforms and cohorts.
2. **Predict** — a small supervised head (a model over a minimal top-K RAV
   panel) maps the projected scores to a phenotype.

A head is an `RAVPhenoModel`: a phenotype name, the RAV prior it was built
on, the RAV columns it uses as features, the fitted model, its `kind`
(`binary` / `multiclass` / `regression` / `cox`), and training metadata.

## Install

```r
# install.packages("devtools")
devtools::install("~/GSS/InferPheno")
```

## Quick start

```r
library(InferPheno)

availablePheno()          # what can I infer today?
res <- inferpheno(expr, "sex")
res                       # prints the head, its OOF performance, and predictions
```

`expr` is a genes-by-samples numeric matrix or `SummarizedExperiment` with
HGNC symbols as rownames. The result is an `inferpheno_result` data.frame —
one row per sample — carrying a per-sample `confidence` column and
`phenotype`, `kind`, and `head_perf` attributes describing the head used.

## Available phenotypes

The shipped heads are fit on the **cohesion-gated de-redundant RAVmodel**
(1,760 RAVs; `rav_prior = "cohesion"`). Each is a non-penalized model over a
minimal top-K-by-|β| RAV panel, refit on all training data.

| Phenotype     | Kind        | Panel K | Out-of-fold performance     |
|---------------|-------------|--------:|-----------------------------|
| `sex`         | binary      |       5 | AUC 0.990                   |
| `msi_status`  | binary      |      20 | AUC 0.951                   |
| `msi_score`   | regression  |      30 | R² 0.617, ρ 0.517           |
| `cms_subtype` | multiclass  |      20 | acc 0.925, macro-F1 0.914   |
| `neutrophil`  | regression  |       3 | R² 0.899, ρ 0.783           |

`availablePheno()` is the live contract — it reads each bundled head and is
the authoritative source for what `inferpheno()` can be called on. Use
`listRAVs(phenotype)` to see the exact RAV panel a head uses.

## Supplying the RAV basis

Heads built on a published GenomicSuperSignature prior (`C2`, `PLIERpriors`)
resolve their basis automatically via `GenomicSuperSignature::getModel()`.

The cohesion basis the shipped heads use (~190 MB) is too large to bundle
and is **not** a `getModel()` prior, so it must be supplied locally. Point
InferPheno at the basis `.rds` once per session:

```r
options(inferpheno.rav_model_path = "/path/to/RAVmodel_cohesion.rds")
# or:   Sys.setenv(INFERPHENO_RAV_MODEL = "/path/to/RAVmodel_cohesion.rds")
# or pass it directly, overriding both:
inferpheno(expr, "sex", rav_model = readRDS("/path/to/RAVmodel_cohesion.rds"))
```

Legacy heads fit on the published **C2** basis (4,764 RAVs) are kept in
`inst/extdata/legacy_c2/` for A/B comparison; they resolve automatically.

## Scope

InferPheno is the **deployment** package: it serves already-built heads and
nothing more. The development and benchmarking tooling lives in two companion
packages:

- **[InferPhenoBuilder](https://github.com/shbrief/InferPhenoBuilder)** —
  produces the heads: cohort assembly and data prep (`project_to_rav()`,
  `join_rav_pheno()`, `standardize_pheno()`, `load_rav_model()`), RAV-panel
  selection, model fitting and CV (`fit_logit()`, `fit_rf()`, `fit_cox()`,
  `cv_oof()`), the RAV × phenotype atlas (`build_atlas()`,
  `score_rav_phenotype()`), and benchmark metrics (`auc_ci()`,
  `paired_delong()`, `brier_ece()`, `multiclass_metrics()`,
  `regression_metrics()`).
- **[InferPhenoPaper](https://github.com/shbrief/InferPhenoPaper)** —
  reproduces the manuscript comparisons: the distribution-shift harness
  (`run_shift()`, `shift_severity()`), the task-specific comparators
  (`run_cmsclassifier()`, `run_cmscaller()`, `run_xist_y_panel()`,
  `run_msi_expression_panel()`), and the per-task benchmark vignettes.

See `quick-start.Rmd` for a hands-on held-out walkthrough (calling CMS on a
metastatic-CRC cohort the heads never saw) and `InferPheno.Rmd` for the
end-user inference workflow.

## Provenance

The cohesion-gated de-redundant RAV basis and the head-fitting methodology
come from
`~/GSS/GenomicSuperSignatureLibrary/refinebioRseq/RAVmodel_536_adaptive_cohesion/cohesion/`
(`METHODS.md`, `RESULTS.md`).

## License

Artistic-2.0
