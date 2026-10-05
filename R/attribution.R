# ============================================================================
# recoverR - attribution.R
#
# Out-of-sample, null-calibrated attribution of omics layers.
#
# THE PROBLEM THIS SOLVES
#   `variancePartition` decomposes variance across NAMED COVARIATES within a
#   fitted model; `MOFA2` decomposes variance across LATENT FACTORS within a
#   fitted model. Both are descriptive decompositions of a model that has
#   already seen all the data. Neither answers the question a user actually
#   asks before building a multi-omics predictor: "if I add this layer, does
#   out-of-sample performance go up, and is any apparent gain distinguishable
#   from noise?"
#
# THE ANSWER IMPLEMENTED HERE
#   Exact Shapley (LMG) decomposition of a cross-validated performance
#   functional, plus a dimension-matched permutation null.
# ============================================================================


# ----------------------------------------------------------------------------
# Performance functionals
# ----------------------------------------------------------------------------

#' Out-of-sample R-squared
#' @keywords internal
perf_r2 <- function(obs, pred) {
  ok <- is.finite(obs) & is.finite(pred)
  if (sum(ok) < 5) return(NA_real_)
  1 - sum((obs[ok] - pred[ok])^2) / sum((obs[ok] - mean(obs[ok]))^2)
}

#' Out-of-sample AUC (rank-based, tie-aware)
#' @keywords internal
perf_auc <- function(obs, pred) {
  ok <- is.finite(obs) & is.finite(pred)
  if (sum(ok) < 5) return(NA_real_)
  if (length(unique(obs[ok])) < 2) return(NA_real_)
  # Mann-Whitney U formulation: AUC = P(pred_case > pred_control) + 0.5 * P(tie).
  # Written directly rather than via pROC because pROC::roc.test() fails at
  # AUC == 1, which is exactly where the strongest markers sit.
  r <- rank(pred[ok], ties.method = "average")
  pos <- obs[ok] == max(obs[ok])
  n1 <- sum(pos); n0 <- sum(!pos)
  if (n1 == 0 || n0 == 0) return(NA_real_)
  (sum(r[pos]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

#' Out-of-sample concordance index (Harrell's C)
#' @keywords internal
perf_cindex <- function(obs, pred, time = NULL, event = NULL) {
  if (is.null(time) || is.null(event)) return(NA_real_)
  ok <- is.finite(pred) & is.finite(time) & !is.na(event)
  if (sum(ok) < 10 || length(unique(event[ok])) < 2) return(NA_real_)
  if (!requireNamespace("survival", quietly = TRUE)) return(NA_real_)
  d <- data.frame(t = time[ok], e = event[ok], v = pred[ok])
  cx <- tryCatch(survival::coxph(survival::Surv(t, e) ~ v, data = d),
                 error = function(e) NULL)
  if (is.null(cx)) return(NA_real_)
  unname(summary(cx)$concordance[1])
}

#' Resolve a performance functional by name
#'
#' @param family One of `"gaussian"`, `"binomial"`, `"survival"`.
#' @return A function `f(obs, pred, ...)` returning a single numeric score,
#'   oriented so that larger is better.
#' @keywords internal
perf_functional <- function(family = c("gaussian", "binomial", "survival")) {
  family <- match.arg(family)
  switch(family,
    gaussian = function(obs, pred, ...) perf_r2(obs, pred),
    binomial = function(obs, pred, ...) perf_auc(obs, pred),
    survival = function(obs, pred, ...) perf_cindex(obs, pred, time = time, event = event)
  )
}


# ----------------------------------------------------------------------------
# Cross-validated performance of one view subset
# ----------------------------------------------------------------------------

#' Stratified fold assignment, shared across every subset
#'
#' Sharing ONE fold vector across all view subsets is both a statistical and a
#' computational requirement. Statistically, comparing V(S) and V(S') computed
#' on different random splits confounds the subset effect with split noise;
#' with n in the tens this dominates. Computationally, it lets the fold
#' structure be computed once.
#'
#' @keywords internal
make_folds <- function(y, nfolds = 5L, seed = 1L, strata = NULL) {
  set.seed(seed)
  n <- length(y)
  fold <- integer(n)
  if (is.null(strata)) {
    fold <- rep(seq_len(nfolds), length.out = n)[sample(n)]
  } else {
    for (lv in unique(strata)) {
      idx <- which(strata == lv)
      fold[idx] <- rep(seq_len(nfolds), length.out = length(idx))[sample(length(idx))]
    }
  }
  fold
}

#' Fit on a training fold and predict a held-out fold
#'
#' Feature standardisation and hyper-parameter selection happen strictly inside
#' the training fold. This is the single most common source of optimistic bias
#' in published multi-omics benchmarks.
#'
#' @keywords internal
fit_predict <- function(Xtr, ytr, Xte, family, alpha_grid = c(0.1, 0.5, 1),
                        nfolds_inner = 5L, seed = 1L, ...) {
  if (ncol(Xtr) < 1L) return(rep(NA_real_, nrow(Xte)))
  # Drop zero-variance columns. A constant (or all-zero) column makes cv.glmnet's
  # lambda path degenerate, which surfaced as a silent NA performance rather than
  # an error -- so a layer that is merely uninformative looked like a numerical
  # failure. Removing the columns is the correct handling.
  vz <- matrixStats::colVars(Xtr, na.rm = TRUE)
  keep <- is.finite(vz) & vz > 0
  if (!any(keep)) return(rep(NA_real_, nrow(Xte)))
  Xtr <- Xtr[, keep, drop = FALSE]
  Xte <- Xte[, keep, drop = FALSE]

  if (requireNamespace("glmnet", quietly = TRUE)) {
    fam <- switch(family, gaussian = "gaussian", binomial = "binomial", survival = "cox")
    yfit <- if (family == "survival") ytr else ytr
    mu <- colMeans(Xtr, na.rm = TRUE)
    sdv <- matrixStats::colSds(Xtr, na.rm = TRUE)
    sdv[!is.finite(sdv) | sdv == 0] <- 1
    Xtr_s <- sweep(sweep(Xtr, 2, mu, "-"), 2, sdv, "/")
    Xte_s <- sweep(sweep(Xte, 2, mu, "-"), 2, sdv, "/")

    best <- NULL
    for (a in alpha_grid) {
      cv <- tryCatch({
        set.seed(seed)
        glmnet::cv.glmnet(Xtr_s, yfit, family = fam, alpha = a,
                          nfolds = nfolds_inner,
                          # NB: cv.glmnet SILENTLY downgrades type.measure from
                          # "auc" to "deviance" when folds hold <10 observations.
                          # Requesting deviance explicitly keeps the behaviour
                          # deterministic instead of depending on n.
                          type.measure = "deviance")
      }, error = function(e) NULL)
      if (is.null(cv)) next
      if (is.null(best) || min(cv$cvm, na.rm = TRUE) < best$cvm) {
        best <- list(cvm = min(cv$cvm, na.rm = TRUE), lambda = cv$lambda.min, alpha = a)
      }
    }
    if (is.null(best)) {
      m <- tryCatch(glmnet::glmnet(Xtr_s, yfit, family = fam, alpha = 0.5), error = function(e) NULL)
      if (is.null(m)) return(rep(NA_real_, nrow(Xte)))
      # A fit given no `lambda` carries the WHOLE path, so predict() on it returns an
      # n_test x length(m$lambda) MATRIX rather than a vector. as.numeric() then flattens it and
      # cv_perf()'s `oof[te] <- ...` RECYCLES the first n_test values -- i.e. it scores the whole
      # held-out fold with the largest-lambda (near-intercept-only) end of the path, an
      # UNDECLARED model, with nothing but a "number of items to replace is not a multiple of
      # replacement length" warning. Measured 2026-09-27 at n = 16: one fold returned 500 values
      # for 5 held-out samples and moved that arm from 0.6000 to 0.6545.
      # The fallback therefore DECLARES one model: the most regularised point of the same path,
      # equivalent to the column the recycling silently used (conservative, never optimistic),
      # and the length is asserted so no future edit can recycle again.
      lam <- max(m$lambda)
      p <- as.numeric(stats::predict(m, Xte_s, s = lam))
      if (length(p) != nrow(Xte)) {
        stop("fit_predict(): the no-lambda fallback returned ", length(p),
             " predictions for ", nrow(Xte), " held-out rows", call. = FALSE)
      }
      return(p)
    }
    m <- glmnet::glmnet(Xtr_s, yfit, family = fam, alpha = best$alpha, lambda = best$lambda)
    return(as.numeric(stats::predict(m, Xte_s)))
  }

  # Base-R fallback: never silently skip, but also never require glmnet.
  if (family == "binomial") {
    d <- data.frame(y = ytr, Xtr, check.names = TRUE)
    m <- tryCatch(stats::glm(y ~ ., data = d, family = stats::binomial()), error = function(e) NULL)
    if (is.null(m)) return(rep(NA_real_, nrow(Xte)))
    return(as.numeric(stats::predict(m, newdata = data.frame(Xte, check.names = TRUE), type = "response")))
  }
  d <- data.frame(y = ytr, Xtr, check.names = TRUE)
  m <- tryCatch(stats::lm(y ~ ., data = d), error = function(e) NULL)
  if (is.null(m)) return(rep(NA_real_, nrow(Xte)))
  as.numeric(stats::predict(m, newdata = data.frame(Xte, check.names = TRUE)))
}

#' Cross-validated performance of a view subset
#'
#' Returns the pooled out-of-fold performance of the subset `S`, i.e. the
#' empirical estimate of `V(S)` used by the Shapley decomposition.
#'
#' @param X_list Named list of view matrices (samples x features).
#' @param y Response: numeric, a two-level factor, or a `Surv` object.
#' @param subset Character vector of view names.
#' @param family `"gaussian"`, `"binomial"` or `"survival"`.
#' @param fold Optional shared fold vector (see [make_folds()]).
#' @return A list with `perf`, the pooled out-of-fold predictions `oof`, and the
#'   resolved `family`.
#' @export
#' @param nfolds Number of cross-validation folds.
#' @param seed RNG seed.
#' @param ... Passed through to the learner.
cv_perf <- function(X_list, y, subset = names(X_list),
                    family = c("binomial", "gaussian", "survival"),
                    fold = NULL, nfolds = 5L, seed = 1L, ...) {
  family <- match.arg(family)
  # Canonical column order: without this, the column ordering of the concatenated
  # matrix depends on how the caller listed the views, and glmnet's internal
  # lambda path is not bit-reproducible across differing column counts. Sorting
  # makes V(S) depend only on the SET of views, never on their listing order.
  subset <- sort(intersect(subset, names(X_list)))
  n <- nrow(X_list[[1]])

  if (family == "survival") {
    time <- as.numeric(y[, 1]); event <- as.integer(y[, 2])
    strata <- event
  } else {
    time <- event <- NULL
    strata <- if (is.factor(y) || is.character(y)) as.character(y) else NULL
  }
  if (is.null(fold)) fold <- make_folds(seq_len(n), nfolds = nfolds, seed = seed, strata = strata)

  oof <- rep(NA_real_, n)
  for (k in sort(unique(fold))) {
    tr <- fold != k; te <- fold == k
    if (sum(tr) < 10L) next
    if (family == "survival") {
      if (length(unique(event[tr])) < 2L) next
      ytr <- y[tr]
    } else if (family == "binomial") {
      if (length(unique(y[tr])) < 2L) next
      ytr <- y[tr]
    } else {
      ytr <- y[tr]
    }
    Xtr <- do.call(cbind, lapply(subset, function(v) X_list[[v]][tr, , drop = FALSE]))
    Xte <- do.call(cbind, lapply(subset, function(v) X_list[[v]][te, , drop = FALSE]))
    p <- fit_predict(Xtr, ytr, Xte, family, seed = seed, nfolds_inner = nfolds, ...)
    # Defence in depth: `oof[te] <- <vector of the wrong length>` RECYCLES in R, which is exactly
    # how an undeclared model came to score a whole held-out fold with only a warning (see
    # fit_predict). The assertion costs nothing and cannot be silently re-broken.
    if (length(p) != sum(te)) {
      stop("cv_perf(): the learner returned ", length(p), " predictions for ", sum(te),
           " held-out rows (fold ", k, ", view subset {", paste(subset, collapse = ","), "})",
           call. = FALSE)
    }
    oof[te] <- p
  }

  f <- switch(family,
    gaussian = function(o, p) perf_r2(o, p),
    binomial = function(o, p) perf_auc(o, p),
    survival = function(o, p) perf_cindex(o, p, time = time, event = event)
  )
  obs <- if (family == "survival") rep(NA_real_, n) else as.numeric(y)
  list(perf = f(obs, oof), oof = oof, family = family, fold = fold,
       n_views = length(subset), subset = subset)
}


# ----------------------------------------------------------------------------
# Exact Shapley (LMG) decomposition
# ----------------------------------------------------------------------------

#' Exact Shapley decomposition of out-of-sample performance
#'
#' For `K` views there are `2^K` subsets. The Shapley value of view `k` is
#'
#' \deqn{\phi_k = \sum_{S \subseteq V \setminus \{k\}}
#'   \frac{|S|!\,(K-|S|-1)!}{K!}\,\bigl[V(S \cup \{k\}) - V(S)\bigr]}
#'
#' which is the unique attribution that is efficient
#' (\eqn{\sum_k \phi_k = V(V) - V(\emptyset)}), symmetric, and order-independent
#' (the LMG / Shapley-value regression decomposition of Lindeman, Merenda and
#' Gold). Because \eqn{V} here is an out-of-sample functional, \eqn{\phi_k} is a
#' **predictive** attribution rather than the within-model variance share
#' reported by `variancePartition` or `MOFA2`.
#'
#' `K <= 6` by default: the cost is exactly `2^K` cross-validated fits, and the
#' `V(S)` values are cached and reused across all `K` views, so the whole
#' decomposition costs one pass over the subset lattice rather than `K` passes.
#'
#' @param X_list Named list of view matrices.
#' @param y Response.
#' @param family Performance functional family.
#' @param max_views Refuse to run above this many views (cost is exponential).
#' @return A `data.frame` of Shapley values with `V(alone)`, `V(full)` and the
#'   unique contribution `V(full) - V(full \ k)`.
#' @export
#' @param nfolds Number of cross-validation folds.
#' @param seed RNG seed.
#' @param ... Passed through to the learner.
shapley_attribution <- function(X_list, y,
                                family = c("binomial", "gaussian", "survival"),
                                nfolds = 5L, seed = 1L, max_views = 6L, ...) {
  family <- match.arg(family)
  views <- names(X_list)
  K <- length(views)
  if (K > max_views) {
    stop("shapley_attribution() enumerates 2^K subsets; K = ", K,
         " exceeds max_views = ", max_views,
         ". Raise max_views deliberately or drop views.", call. = FALSE)
  }

  strata <- if (family == "survival") as.integer(y[, 2]) else
    if (is.factor(y) || is.character(y)) as.character(y) else NULL
  fold <- make_folds(seq_len(nrow(X_list[[1]])), nfolds = nfolds, seed = seed, strata = strata)

  # --- cache V(S) over the subset lattice ------------------------------------
  # PITFALL: subset keys must be canonical in VIEW ORDER, not sorted
  # alphabetically. combn() emits combinations in `views` order, while c(s, v)
  # inside the Shapley loop appends the target view last; sorting both would not
  # reconcile them, and a missing key silently fails as "subscript out of
  # bounds" rather than as a wrong number. Rank by position in `views`.
  vrank <- stats::setNames(seq_along(views), views)
  key <- function(s) {
    s <- s[!is.na(s)]
    # NB: an empty-string list name cannot be assigned in R
    # (subsets[[""]] <- x errors), so the empty set gets a sentinel key.
    if (!length(s)) return("<empty>")
    paste(s[order(vrank[s])], collapse = "+")
  }

  subsets <- list()
  for (m in 0:K) {
    cmb <- if (m == 0L) list(character(0)) else utils::combn(views, m, simplify = FALSE)
    for (s in cmb) subsets[[key(s)]] <- s
  }
  if (length(subsets) != 2^K) {
    stop("subset lattice is incomplete: ", length(subsets), " of ", 2^K,
         " subsets keyed; view names must be unique.", call. = FALSE)
  }
  V <- vapply(subsets, function(s) {
    if (length(s) == 0L) return(0)              # V(emptyset) := 0 by definition
    cv_perf(X_list, y, subset = s, family = family, fold = fold,
            nfolds = nfolds, seed = seed, ...)$perf
  }, numeric(1))

  # --- Shapley values --------------------------------------------------------
  fact <- factorial(0:K)
  phi <- setNames(numeric(K), views)
  for (i in seq_along(views)) {
    v <- views[i]
    others <- setdiff(views, v)
    acc <- 0
    for (m in 0:(K - 1L)) {
      for (s in utils::combn(others, m, simplify = FALSE)) {
        w <- fact[m + 1L] * fact[K - m] / fact[K + 1L]
        acc <- acc + w * (V[[key(c(s, v))]] - V[[key(s)]])
      }
    }
    phi[i] <- acc
  }

  V_full <- V[[key(views)]]
  alone <- vapply(views, function(v) V[[key(v)]], numeric(1))
  uniq  <- vapply(views, function(v) V_full - V[[key(setdiff(views, v))]], numeric(1))

  data.frame(
    view = views,
    shapley = as.numeric(phi),
    V_alone = as.numeric(alone),
    unique_contribution = as.numeric(uniq),
    V_full = V_full,
    V_empty = 0,
    n_views = K,
    stringsAsFactors = FALSE
  )[order(-phi), ]
}


# ----------------------------------------------------------------------------
# Dimension-matched permutation null
# ----------------------------------------------------------------------------

#' Permutation null for one view's marginal contribution
#'
#' Rows of `X_list[[view]]` are permuted, which re-pairs each sample with
#' another sample's profile for that view ONLY. This destroys the
#' view-to-response association while preserving (a) the view's internal
#' covariance structure and marginal distributions, and (b) every other view's
#' alignment with the response. The null hypothesis is therefore precisely
#' "this view carries no sample-specific information about the response,
#' given the others" - which is the hypothesis a user cares about.
#'
#' The null is dimension-matched by construction: the permuted view has exactly
#' the same number of features as the real one, so no correction for
#' dimensionality is needed (unlike comparing against a random-noise view).
#'
#' @param B Number of permutations.
#' @return A data.frame with the observed `delta`, the null mean/sd, and a
#'   one-sided p-value of the form (1 + number of null values at least as
#'   extreme) / (B + 1).
#' @export
#' @param X_list Named list of view matrices (samples x features), rows aligned across views.
#' @param y Response: a two-level factor (binomial), a numeric vector (gaussian), or a `Surv` object (survival).
#' @param view Name of the view to permute.
#' @param family Performance functional family.
#' @param nfolds Number of cross-validation folds.
#' @param seed RNG seed.
#' @param ncores Number of worker processes for the permutation null. Seeds are pre-generated, so results are numerically identical to a serial run.
#' @param ... Passed through to the learner.
permutation_null <- function(X_list, y, view, family = c("binomial", "gaussian", "survival"),
                             B = 199L, nfolds = 5L, seed = 1L, ncores = 1L, ...) {
  family <- match.arg(family)
  if (!view %in% names(X_list)) stop("unknown view: ", view, call. = FALSE)
  views <- names(X_list)
  others <- setdiff(views, view)

  strata <- if (family == "survival") as.integer(y[, 2]) else
    if (is.factor(y) || is.character(y)) as.character(y) else NULL
  fold <- make_folds(seq_len(nrow(X_list[[1]])), nfolds = nfolds, seed = seed, strata = strata)

  delta_of <- function(Xl) {
    v_full <- cv_perf(Xl, y, views, family = family, fold = fold, nfolds = nfolds, seed = seed, ...)$perf
    v_drop <- cv_perf(Xl, y, others, family = family, fold = fold, nfolds = nfolds, seed = seed, ...)$perf
    v_full - v_drop
  }

  obs <- delta_of(X_list)
  n <- nrow(X_list[[1]])

  # The null is the computational bottleneck: each replicate costs two full
  # cross-validated evaluations, so B = 199 over 4 views is ~1600 fits. Seeds are
  # pre-generated so that parallel and serial runs give IDENTICAL draws --
  # otherwise a parallel run would silently return different numbers from the
  # serial one, which is worse than being slow.
  seeds <- sample.int(.Machine$integer.max, B)
  one <- function(b) {
    set.seed(seeds[b])
    Xp <- X_list
    Xp[[view]] <- Xp[[view]][sample(n), , drop = FALSE]
    tryCatch(delta_of(Xp), error = function(e) NA_real_)
  }

  nul <- if (ncores > 1L && .Platform$OS.type != "windows" &&
             requireNamespace("parallel", quietly = TRUE)) {
    unlist(parallel::mclapply(seq_len(B), one, mc.cores = min(ncores, B),
                              mc.preschedule = TRUE), use.names = FALSE)
  } else {
    vapply(seq_len(B), one, numeric(1))
  }
  nul <- nul[is.finite(nul)]

  data.frame(
    view = view,
    delta_observed = obs,
    null_mean = mean(nul),
    null_sd = stats::sd(nul),
    null_q95 = as.numeric(stats::quantile(nul, 0.95, na.rm = TRUE)),
    p_value = (1 + sum(nul >= obs)) / (1 + length(nul)),
    B = length(nul),
    stringsAsFactors = FALSE
  )
}


# ----------------------------------------------------------------------------
# User-facing entry point
# ----------------------------------------------------------------------------

#' Attribute out-of-sample performance across omics layers
#'
#' @param X_list Named list of view matrices (samples x features, rows aligned).
#' @param y Response.
#' @param family Performance functional.
#' @param n_perm Number of permutations for the marginal null (0 to skip).
#' @param verbose Print progress.
#' @return An object of class `recoverR_attribution` with the Shapley
#'   decomposition, the per-view permutation null, and the ablation table.
#' @export
#' @param nfolds Number of cross-validation folds.
#' @param seed RNG seed.
#' @param max_views Refuse to run above this many views; the cost is 2^K cross-validated fits.
#' @param ncores Number of worker processes for the permutation null. Seeds are pre-generated, so results are numerically identical to a serial run.
#' @param ... Passed through to the learner.
attrib_views <- function(X_list, y,
                         family = c("binomial", "gaussian", "survival"),
                         nfolds = 5L, seed = 1L, n_perm = 199L,
                         max_views = 6L, ncores = 1L, verbose = TRUE, ...) {
  family <- match.arg(family)
  stopifnot(is.list(X_list), length(X_list) >= 1L, !is.null(names(X_list)))
  n <- vapply(X_list, nrow, integer(1))
  if (length(unique(n)) != 1L) {
    stop("all views must have the same number of rows; got ",
         paste(unique(n), collapse = ", "), call. = FALSE)
  }

  if (verbose) message("Shapley decomposition over ", length(X_list), " views (",
                       length(X_list)^2, " subsets max)...")
  sh <- shapley_attribution(X_list, y, family = family, nfolds = nfolds,
                            seed = seed, max_views = max_views, ...)

  pn <- NULL
  if (n_perm > 0L) {
    if (verbose) message("Permutation null (B = ", n_perm, ") per view...")
    pn <- do.call(rbind, lapply(names(X_list), function(v) {
      tryCatch(permutation_null(X_list, y, v, family = family, B = n_perm,
                                nfolds = nfolds, seed = seed, ncores = ncores, ...),
               error = function(e) NULL)
    }))
  }

  ab <- tryCatch(ablate_views(X_list, y, family = family, nfolds = nfolds, seed = seed, ...),
                 error = function(e) NULL)

  structure(list(shapley = sh, null = pn, ablation = ab,
                 family = family, n_views = length(X_list)),
            class = "recoverR_attribution")
}

#' Forward and leave-one-view-out ablation
#'
#' Forward selection orders views by their marginal gain given what has already
#' been added; leave-one-out removes each view from the full model. Reading the
#' two together distinguishes a view that is useful on its own from one that is
#' useful only in combination, and exposes views that are actively harmful.
#'
#' @export
#' @param X_list Named list of view matrices (samples x features), rows aligned across views.
#' @param y Response: a two-level factor (binomial), a numeric vector (gaussian), or a `Surv` object (survival).
#' @param family Performance functional family.
#' @param nfolds Number of cross-validation folds.
#' @param seed RNG seed.
#' @param ... Passed through to the learner.
ablate_views <- function(X_list, y, family = c("binomial", "gaussian", "survival"),
                         nfolds = 5L, seed = 1L, ...) {
  family <- match.arg(family)
  views <- names(X_list)
  strata <- if (family == "survival") as.integer(y[, 2]) else
    if (is.factor(y) || is.character(y)) as.character(y) else NULL
  fold <- make_folds(seq_len(nrow(X_list[[1]])), nfolds = nfolds, seed = seed, strata = strata)

  vp <- function(s) {
    if (length(s) == 0L) return(0)
    cv_perf(X_list, y, s, family = family, fold = fold, nfolds = nfolds, seed = seed, ...)$perf
  }

  V_full <- vp(views)
  loo <- data.frame(
    view = views,
    V_without = vapply(views, function(v) vp(setdiff(views, v)), numeric(1)),
    stringsAsFactors = FALSE
  )
  loo$loss_when_removed <- V_full - loo$V_without

  chosen <- character(0); remaining <- views; fwd <- list(); cur <- 0
  repeat {
    if (length(remaining) == 0L) break
    gains <- vapply(remaining, function(v) {
      tryCatch(vp(c(chosen, v)) - cur, error = function(e) NA_real_)
    }, numeric(1))
    if (all(!is.finite(gains))) break
    best <- remaining[which.max(gains)]
    cur <- vp(c(chosen, best))
    chosen <- c(chosen, best)
    fwd[[length(fwd) + 1L]] <- data.frame(step = length(chosen), view = best,
                                          gain = max(gains, na.rm = TRUE), V = cur)
    remaining <- setdiff(remaining, best)
  }
  forward <- if (length(fwd)) do.call(rbind, fwd) else NULL

  list(V_full = V_full, leave_one_out = loo[order(-loo$loss_when_removed), ],
       forward = forward)
}

#' @export
print.recoverR_attribution <- function(x, ...) {
  cat("<recoverR_attribution>", x$n_views, "views |", x$family, "endpoint\n\n")
  cat("Shapley decomposition of out-of-sample performance:\n")
  s <- x$shapley
  s$shapley <- round(s$shapley, 4); s$V_alone <- round(s$V_alone, 4)
  s$unique_contribution <- round(s$unique_contribution, 4)
  print(s[, c("view", "shapley", "V_alone", "unique_contribution")], row.names = FALSE)
  if (!is.null(x$null)) {
    cat("\nDimension-matched permutation null:\n")
    n <- x$null
    for (i in seq_len(nrow(n))) {
      cat(sprintf("  %-14s delta = %+.4f   p = %.4f  (B = %d)\n",
                  n$view[i], n$delta_observed[i], n$p_value[i], n$B[i]))
    }
  }
  if (!is.null(x$ablation)) {
    cat("\nV(full) =", round(x$ablation$V_full, 4), "\n")
    cat("Leave-one-view-out loss:\n")
    l <- x$ablation$leave_one_out
    l$loss_when_removed <- round(l$loss_when_removed, 4)
    l$V_without <- round(l$V_without, 4)
    print(l, row.names = FALSE)
  }
  invisible(x)
}
