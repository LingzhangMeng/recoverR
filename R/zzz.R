# ============================================================================
# recoverR - zzz.R
#
# Package-level declarations.
# ============================================================================

#' @importFrom stats setNames reorder
#' @keywords internal
NULL

# ggplot2's non-standard evaluation means column names appear as bare symbols in
# aes(); R CMD check cannot see them and reports "no visible binding for global
# variable". Declaring them here is the idiomatic fix and keeps the check clean
# without suppressing anything broadly.
utils::globalVariables(c(
  "V", "value", "label", "fill", "answered", "series", "step", "time", "event",
  "shapley", "view", "V_full", "loss_when_removed", "effect", "verdict",
  "cv_r2", "layer", "mean_conf", "observed_acc", "n", "confusability",
  "confidence", "set_size", "set_size_f", "abstain", "to", "from",
  "diameter_mm", "performance", "panel_f", "score", "cluster", "mean_z",
  "Freq", "phenotype", "z", "auc", "lab", "cohort", "biomarker",
  "ci_lo", "ci_hi", "mean_auc", "c_index", "set", "os_cindex", "os_logrank_p",
  "ablation", "n_views", "delta", "gene", "r_meth_expr", "p_value",
  "null_q95", "delta_observed", "panel_genes"
))
