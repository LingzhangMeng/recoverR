# ============================================================================
# recoverR - guardrails.R
#
# Every function here exists because a real analysis silently produced a wrong
# answer without raising an error. They are exported because the failure modes
# are general: they will bite anyone doing cross-validated multi-omics work, and
# a package that only documents them in a vignette does not prevent them.
#
# Design rule: these wrappers WARN or STOP. They never silently correct a
# plausible-looking but wrong result, because a silent correction is how the
# original error survived to the manuscript in the first place.
# ============================================================================


#' Safe `cv.glmnet`
#'
#' `glmnet::cv.glmnet()` silently downgrades `type.measure = "auc"` to
#' `"deviance"` when a fold holds fewer than 10 observations, emitting only a
#' warning. `cv$cvm` then contains DEVANCE rather than AUC, so the idiom
#' `which.max(cv$cvm)` maximises *error*. At the sample sizes typical of
#' immune-oncology cohorts (n = 16-30) this downgrade is the norm, not the
#' exception.
#'
#' This wrapper requests deviance explicitly (so behaviour does not depend on n)
#' and returns no `cvm`-as-AUC ambiguity for callers to trip over.
#'
#' @param x Predictor matrix.
#' @param y Response.
#' @param family glmnet family.
#' @param alpha Elastic-net mixing parameter.
#' @param nfolds Number of folds.
#' @param type.measure Selection measure. Defaults to `"deviance"`.
#' @return The `cv.glmnet` object, with an attribute recording the measure used.
#' @export
#' @param ... Passed through to the learner.
cv_glmnet_safe <- function(x, y, family = "binomial", alpha = 0.5,
                           nfolds = 5L, type.measure = "deviance", ...) {
  if (!requireNamespace("glmnet", quietly = TRUE)) {
    stop("cv_glmnet_safe() needs the 'glmnet' package.", call. = FALSE)
  }
  if (identical(type.measure, "auc")) {
    warning("type.measure = 'auc' will be downgraded to 'deviance' by ",
            "cv.glmnet whenever a fold holds <10 observations, making cv$cvm ",
            "a deviance rather than an AUC. Use type.measure = 'deviance' and ",
            "compute AUC yourself, or call recoverR::perf_auc() on held-out ",
            "predictions.", call. = FALSE)
  }
  min_fold <- floor(length(y) / nfolds)
  if (min_fold < 10L && identical(type.measure, "auc")) {
    type.measure <- "deviance"
  }
  out <- glmnet::cv.glmnet(x, y, family = family, alpha = alpha,
                           nfolds = nfolds, type.measure = type.measure, ...)
  attr(out, "recoverR_type_measure") <- type.measure
  out
}


#' Safe AUC
#'
#' `pROC::roc.test()` errors when the ROC curve is perfect (AUC == 1), which is
#' precisely where the strongest biomarkers sit, so its p-value comes back `NA`
#' and propagates `NA` through any downstream false-discovery adjustment. The
#' Mann-Whitney (rank-sum) test is the exact nonparametric test of AUC = 0.5, is
#' tie-aware, and never degenerates.
#'
#' @param score Numeric scores.
#' @param label Two-level outcome; the second level is treated as the positive
#'   class.
#' @return A list with `auc` and `p` (Mann-Whitney, two-sided).
#' @export
auc_safe <- function(score, label) {
  ok <- is.finite(score) & !is.na(label)
  score <- score[ok]; label <- as.factor(label[ok])
  if (length(unique(label)) < 2L) return(list(auc = NA_real_, p = NA_real_))
  a <- perf_auc(as.numeric(label == levels(label)[2]), score)
  p <- tryCatch(
    stats::wilcox.test(score ~ label, exact = FALSE)$p.value,
    error = function(e) NA_real_
  )
  list(auc = a, p = p)
}


#' Guards for the structural mistakes
#'
#' @name structural_guards
#' @keywords internal
NULL

#' Assert that column names are unique across concatenated blocks
#'
#' `cbind()`ing per-layer principal components is natural and wrong: every layer
#' names its columns `PC1..PCn`, so the combined frame has duplicate names and
#' `lm()` then returns `NA` for the combined model - silently turning the
#' headline comparison into a row of `NA`s.
#'
#' @param ... Matrices or data frames.
#' @return The column-bound matrix with layer-prefixed, unique names.
#' @export
safe_cbind <- function(...) {
  blocks <- list(...)
  nms <- names(blocks)
  if (is.null(nms)) nms <- paste0("block", seq_along(blocks))
  blocks <- lapply(seq_along(blocks), function(i) {
    b <- as.matrix(blocks[[i]])
    colnames(b) <- paste0(nms[i], "_", colnames(b))
    b
  })
  out <- do.call(cbind, blocks)
  if (anyDuplicated(colnames(out))) {
    stop("duplicate column names survived prefixing; fix the input names.",
         call. = FALSE)
  }
  out
}

#' Assert that no row is entirely missing
#'
#' `if (x) next` over a vector that contains `NA` fails with "missing value where
#' TRUE/FALSE needed". The usual cause is an `all.x = TRUE` merge: samples with
#' no annotation carry `NA`, and a filter written for numeric data silently
#' becomes a filter over `NA`.
#'
#' @param x Vector.
#' @param what Label used in the error message.
#' @export
assert_no_na <- function(x, what = "x") {
  if (anyNA(x)) {
    n <- sum(is.na(x))
    stop(sprintf("%s contains %d NA value(s) (%d%%). Filter them explicitly or ",
                 what, n, round(100 * n / length(x))),
         "replace them with an intended value; do not let them reach a boolean ",
         "condition.", call. = FALSE)
  }
  invisible(TRUE)
}

#' Assert that a scale has enough discrete values
#'
#' `paste0("C", NA)` yields the literal string `"CNA"`, adding a phantom factor
#' level. A `scale_fill_manual()` with 3 colours then fails at DRAW time with
#' "Insufficient values in manual scale. 4 needed but only 3 provided" - long
#' after the mistake was made.
#'
#' @param levels Character vector of levels actually present.
#' @param values Palette or named vector.
#' @export
assert_palette_covers <- function(levels, values) {
  lv <- unique(as.character(levels))
  lv <- lv[!is.na(lv)]
  if (length(lv) > length(values)) {
    stop(sprintf("palette has %d value(s) but %d level(s) are present (%s). ",
                 length(values), length(lv), paste(lv, collapse = ", ")),
         "A common cause is a paste0() branch that produced a literal 'NA' ",
         "level.", call. = FALSE)
  }
  invisible(TRUE)
}

#' Detect an accidental literal "NA" level
#'
#' @param x Character vector.
#' @export
detect_cna_levels <- function(x) {
  hits <- unique(as.character(x)[grepl("NA$", as.character(x)) & !is.na(x)])
  if (length(hits)) {
    warning("detected suspicious level(s) ending in 'NA': ",
            paste(hits, collapse = ", "),
            ". This usually comes from paste0(prefix, NA).", call. = FALSE)
  }
  invisible(hits)
}
