# InferPheno 0.99.2 (2026-09-12)

## All heads sized by one procedure (`asymptote_cal`)

Every shipped head is now sized by the **same** rule — `asymptote_cal`: hard
calibration and EPV-ceiling gates, *K* at the discrimination asymptote, bootstrap
reproducibility demoted to an advisory — at one protocol and one seed (B = 200,
tau = 0.6, 25 x 5-fold CV with in-fold re-ranking). Previously the shipped set was
a mix of rules and vintages. `msi_status` (K = 20) and `msi_score` (K = 50) were
already consistent and are unchanged, reproducing bit-for-bit.

| Head | Before | Now | Primary metric |
|---|---|---|---|
| `sex` | K = 8 | **K = 5** | AUC 0.99993 -> 0.99938 |
| `neutrophil` | K = 10 | **K = 5** | R2 0.98679 -> 0.97847 |
| `cms_subtype` | K = 3, n = 1355 | **K = 10, n = 1517** | accuracy 0.71365 -> **0.89189** |
| `sle_neutrophil` | *not shipped* | **K = 10** (new) | R2 0.66206 |

For `sex` and `neutrophil` the new panel is a strict prefix of the old one: the
ranking is unchanged and the rule simply stops earlier. The small in-cohort loss
is the deliberate parsimony trade, not a regression.

## `cms_subtype` retrained on corrected data

The previous head was trained on 1,355 samples because the CMS loader dropped every
TCGA sample through a dotted-vs-hyphenated `sample_id` mismatch. With the IDs
corrected the training set is 1,517 across all 8 cohorts, and accuracy rises from
0.714 to 0.892 (macro-F1 0.674 -> 0.871). The panel keeps all three previous RAVs
and adds seven.

## New head: `sle_neutrophil`

Neutrophil counts from SLE whole blood (E-GEOD-65391, 853 samples), K = 10,
R2 0.662. Two caveats are recorded in `availablePheno()$transfer_note`: it has no
external cohort (SLE-WB is itself the external set for `neutrophil`, so the two
cannot validate each other), and the cohort is longitudinal — 853 samples are 156
participants — while the head was sized with folds assigned per sample, so R2 0.662
is mildly optimistic (0.645 under participant-grouped folds).

## External transfer figures withheld for two heads

`cms_subtype` and `neutrophil` changed enough that their recorded external scores
(0.653 and 0.610) no longer describe the shipped head. `transfer_perf` is now `NA`
for both, with the provenance in `transfer_note`, rather than carrying a stale
number forward. Re-scoring against the external gate is pending.

## Bundled basis regenerated

`cohesion_subset.rds` covers the union of the six shipped panels: **97 RAVs**
(was 88). Prior heads and basis are preserved under
`inst/extdata/legacy_2026-09-12/`.

# InferPheno 0.99.1 (2026-07-23)

## Log-native projection (default behaviour change)

`inferpheno()` and `projectExpression()` now project expression on **log2 scale**
by default (`log_input = "auto"`). `calculateScore()` row-centres each gene and
applies no internal log, while the RAV basis is log-native, so linear input (raw
counts / RSEM / TPM) must be log-transformed before projection. The new default
is **scale-aware**: it applies `log2(x + 1)` only when the input looks linear
(99th percentile > ~30) and a genuine expression matrix (no negative values),
leaving already-log input (microarray, log2-CPM) and centred/z-scored input
untouched — so it never double-logs. Override with `log_input = "yes"` / `"no"`.

## Heads rebuilt on the log-native projection

The RSEM-trained heads (`sex`, `msi_status`, `msi_score`) were rebuilt on the
log-native projection; `cms_subtype` and `neutrophil` were unchanged (their
training inputs are centred / already-log, so auto-detect leaves them). The
bundled `cohesion_subset` basis was regenerated to cover the new panels.

## `sex` head now TRANSFERS (verdict corrected)

The `sex` head previously carried `transfer = "fails"` (external GTEx AUC 0.540),
attributed to cohort confounding (sex ~ cancer type in TCGA training). Re-analysis
showed the failure was in large part a **projection-scale artifact**: stability
selection on the raw-linear projection picked cancer-type-proxy RAVs rather than
sex RAVs. Rebuilding on the log-native projection selects genuine Y-chromosome
RAVs (top loadings on EIF1AY / DDX3Y / KDM5D), and the head now transfers to GTEx
normal tissue at **AUC 0.940**. `availablePheno()` and the runtime warning are
updated accordingly. `msi_status` external transfer is unchanged (0.824), i.e.
scale-invariant.
