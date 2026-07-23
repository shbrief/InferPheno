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
