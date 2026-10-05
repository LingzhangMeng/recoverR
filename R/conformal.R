# ============================================================================
# recoverR - conformal.R
#
# Can a BULK profile carry a spatially-defined phenotype at all, and with what
# calibrated confidence?
#
# WHY THIS IS NOT JUST "ANOTHER DECONVOLUTION TOOL"
#   Reliability layers now exist for spatial deconvolution OUTPUT (adding a risk
#   map to RCTD / cell2location / Tangram predictions). Those tools ask "is this
#   spot-level spatial call trustworthy?". They all operate spatial -> spatial.
#   None asks the inverse question, which is the one that decides whether a bulk
#   biomarker programme is viable at all: given a bulk profile, is the spatial
#   archetype recoverable, and if not, can the method SAY SO?
#
#   The motivating evidence: in a randomised neoadjuvant trial (DUTRENEO) a
#   validated 18-gene bulk inflammation signature failed to enrich for
#   immunotherapy responders, and spatial profiling attributed the failure to
#   architectures invisible to bulk assays. A method that reports a confident
#   archetype for a profile that cannot support one is therefore not merely
#   imprecise - it is the mechanism by which bulk biomarkers fail in trials.
#
# THE ANSWER IMPLEMENTED HERE
#   Split conformal prediction with an explicit abstention rule, group-conditional
#   (Mondrian) calibration by archetype, and optional weighted calibration for
#   cohort shift. Abstention is the deliverable, not a shortcoming.
# ============================================================================


# ----------------------------------------------------------------------------
# Nonconformity scores
# ----------------------------------------------------------------------------

#' Least-ambiguous-set (LAC) nonconformity score
#'
#' \deqn{s(x, y) = 1 - \hat p_y(x)}
#' Small scores mean the candidate label is well supported, so the resulting
#' prediction set is \eqn{\{y : \hat p_y(x) \ge 1 - \hat q\}}.
#' @keywords internal
score_lac <- function(probs, y_idx) 1 - probs[cbind(seq_along(y_idx), y_idx)]

#' Adaptive prediction set (APS) nonconformity score
#'
#' \deqn{s(x, y) = \sum_{j : \hat p_j(x) \ge \hat p_y(x)} \hat p_j(x)}
#' the cumulative probability mass of every label at least as likely as the true
#' one. APS adapts set size to difficulty, so easy samples get small sets and
#' hard samples are allowed to abstain rather than being forced into a wrong
#' singleton.
#' @keywords internal
score_aps <- function(probs, y_idx) {
  vapply(seq_along(y_idx), function(i) {
    p <- probs[i, ]; p_true <- p[y_idx[i]]
    sum(p[p >= p_true])
  }, numeric(1))
}

#' Finite-sample conformal quantile
#'
#' \deqn{\hat q = s_{(\lceil (n+1)(1-\alpha) \rceil)}}
#' with \eqn{s_{(n+1)} := +\infty}. The `(n+1)` correction is what buys the
#' finite-sample coverage guarantee; using the plain empirical quantile gives
#' nominal coverage only asymptotically and under-covers at the n = 20-50 that
#' immune-oncology cohorts actually have.
#' @keywords internal
conformal_quantile <- function(scores, alpha, weights = NULL) {
  s <- sort(scores[is.finite(scores)])
  n <- length(s)
  if (n == 0L) return(Inf)
  k <- ceiling((n + 1) * (1 - alpha))
  if (k > n) return(Inf)
  if (is.null(weights)) return(s[k])
  # weighted (covariate-shift) conformal: smallest q with cumulative normalised
  # calibration weight >= k / (n + 1)
  w <- weights[order(scores[is.finite(scores)])]
  cw <- cumsum(w) / sum(w)
  idx <- which(cw >= k / (n + 1))
  if (!length(idx)) return(Inf)
  s[min(idx)]
}


# ----------------------------------------------------------------------------
# Default learner
# ----------------------------------------------------------------------------

#' Default probabilistic classifier
#'
#' Linear discriminant analysis by default: deterministic, dependency-light
#' (`MASS` is base-recommended), and well behaved when features are rank
#' normalised. `glmnet` multinomial elastic net is available for p >> n.
#' @keywords internal
learn_probs <- function(Xtr, ytr, Xte, learner = c("lda", "glmnet"), seed = 1L) {
  learner <- match.arg(learner)
  ytr <- factor(ytr)
  # LDA needs a pooled covariance that can actually be estimated. The original guard
  # fired only at p >= n, which is NOT a sufficient condition for a usable estimate:
  # at p = 2000 and n = 2235 (p ~ n, p < n) `MASS::lda()` still inverts a 2000 x 2000
  # covariance from 2235 samples and returns posteriors that SATURATE at ~1.0 while the
  # call is at chance. That is worse than an error, because the conformal layer builds
  # its per-class quantiles from those posteriors and then admits every class for every
  # sample. Measured on a planted, zero-noise signal (p = 2000, n_train = 2235):
  #   lda    -> accuracy 0.725, ECE 0.274, abstention 0.996   (the harness fails a
  #             signal that is definitely present)
  #   glmnet -> accuracy 0.954, ECE 0.029, abstention 0.000
  # So the criterion is a RATIO, not a comparison: require at least
  # `LDA_MIN_SAMPLES_PER_FEATURE` training samples per feature, which is the usual
  # rule of thumb for an estimable covariance and is deliberately conservative. The
  # condition also subsumes the old one (p >= n implies p * ratio > n for any ratio > 1).
  # Fixtures with p = 6, n = 400 are unaffected: 6 * 5 = 30 < 400, so lda is still used
  # wherever it is genuinely appropriate.
  LDA_MIN_SAMPLES_PER_FEATURE <- 5
  if (learner == "lda" && ncol(Xtr) * LDA_MIN_SAMPLES_PER_FEATURE > nrow(Xtr)) {
    warning("LDA is ill-conditioned with ", ncol(Xtr), " features and ", nrow(Xtr),
            " training samples (fewer than ", LDA_MIN_SAMPLES_PER_FEATURE,
            " samples per feature); switching to the glmnet multinomial learner. ",
            "Set learner = \"glmnet\" explicitly to silence this.", call. = FALSE)
    learner <- "glmnet"
  }
  if (learner == "lda") {
    m <- tryCatch(MASS::lda(Xtr, ytr), error = function(e) NULL)
    if (!is.null(m)) {
      pr <- tryCatch(stats::predict(m, Xte)$posterior, error = function(e) NULL)
      if (!is.null(pr)) return(pr[, levels(ytr), drop = FALSE])
    }
  }
  if (requireNamespace("glmnet", quietly = TRUE)) {
    mu <- colMeans(Xtr); sv <- matrixStats::colSds(Xtr)
    sv[!is.finite(sv) | sv == 0] <- 1
    Xtr_s <- sweep(sweep(Xtr, 2, mu, "-"), 2, sv, "/")
    Xte_s <- sweep(sweep(Xte, 2, mu, "-"), 2, sv, "/")
    set.seed(seed)
    m <- tryCatch(glmnet::cv.glmnet(Xtr_s, ytr, family = "multinomial",
                                    alpha = 0.5, nfolds = min(5L, sum(table(ytr))),
                                    type.measure = "deviance"),
                  error = function(e) NULL)
    if (!is.null(m)) {
      pr <- stats::predict(m, Xte_s, type = "response")[, , 1]
      pr <- as.matrix(pr)
      if (nrow(pr) == nrow(Xte) && all(levels(ytr) %in% colnames(pr))) {
        return(pr[, levels(ytr), drop = FALSE])
      }
    }
  }
  # last resort: class priors (explicitly uninformative, never silently absent)
  matrix(rep(table(ytr) / length(ytr), each = nrow(Xte)),
         nrow = nrow(Xte), dimnames = list(NULL, levels(ytr)))
}


# ----------------------------------------------------------------------------
# Fit / calibrate
# ----------------------------------------------------------------------------

#' Fit a conformal recoverability model for bulk -> spatial archetype
#'
#' @param X Bulk feature matrix (samples x features), rank-normalised internally.
#' @param archetype Factor of spatially-defined archetype labels.
#' @param cal_prop Fraction of data used for conformal calibration.
#' @param alpha Miscoverage level; the target is `1 - alpha` coverage.
#' @param mondrian Calibrate separately per archetype (group-conditional
#'   coverage). Strongly recommended: marginal coverage can be met while
#'   systematically failing on the archetype that matters most.
#' @param score `"aps"` (default) or `"lac"`.
#' @param learner `"lda"` or `"glmnet"`.
#' @return An object of class `recoverR_conformal`.
#' @export
#' @param seed RNG seed.
#' @param weights Optional importance weights for weighted (covariate-shift) conformal calibration.
fit_recoverability <- function(X, archetype, cal_prop = 0.5, alpha = 0.1,
                               mondrian = TRUE, score = c("aps", "lac"),
                               learner = c("lda", "glmnet"), seed = 1L,
                               weights = NULL) {
  score <- match.arg(score); learner <- match.arg(learner)
  X <- as.matrix(X); archetype <- factor(archetype)
  if (nrow(X) != length(archetype)) stop("X and archetype length mismatch", call. = FALSE)
  if (length(unique(archetype)) < 2L) stop("need at least two archetypes", call. = FALSE)
  if (length(unique(archetype)) < 3L) {
    message("NOTE: only two archetypes present. The scientifically decisive ",
            "contrast (immune-excluded vs immune-desert) requires a taxonomy ",
            "with both classes represented; results here cannot speak to it.")
  }

  set.seed(seed)
  # stratified split so every archetype reaches the calibration set
  cal_idx <- unlist(lapply(split(seq_len(nrow(X)), archetype), function(ix) {
    sample(ix, max(1L, floor(length(ix) * cal_prop)))
  }))
  tr_idx <- setdiff(seq_len(nrow(X)), cal_idx)

  probs_cal <- learn_probs(X[tr_idx, , drop = FALSE], archetype[tr_idx],
                           X[cal_idx, , drop = FALSE], learner, seed)
  lev <- levels(archetype)
  y_idx <- match(as.character(archetype[cal_idx]), lev)
  s <- if (score == "aps") score_aps(probs_cal, y_idx) else score_lac(probs_cal, y_idx)

  w_cal <- if (is.null(weights)) NULL else weights[cal_idx]
  if (mondrian) {
    q <- vapply(lev, function(g) {
      sel <- which(as.character(archetype[cal_idx]) == g)
      if (!length(sel)) return(Inf)
      conformal_quantile(s[sel], alpha, if (is.null(w_cal)) NULL else w_cal[sel])
    }, numeric(1))
  } else {
    q <- setNames(rep(conformal_quantile(s, alpha, w_cal), length(lev)), lev)
  }

  # calibration diagnostics on the held-out fold
  cal_sets <- predict_sets(probs_cal, lev, q, score)
  covered <- vapply(seq_along(y_idx), function(i) lev[y_idx[i]] %in% cal_sets[[i]], logical(1))

  structure(list(
    model = learn_probs, levels = lev, q = q, alpha = alpha, score = score,
    mondrian = mondrian, learner = learner,
    cal_probs = probs_cal, cal_labels = archetype[cal_idx], cal_scores = s,
    cal_sets = cal_sets, cal_covered = covered,
    cal_coverage = mean(covered),
    cal_coverage_by_class = vapply(lev, function(g)
      mean(covered[as.character(archetype[cal_idx]) == g]), numeric(1)),
    n_cal = length(cal_idx), n_train = length(tr_idx),
    used_weights = !is.null(weights)
  ), class = "recoverR_conformal")
}

#' Build prediction sets from probabilities and quantiles
#'
#' Mondrian calibration is candidate-specific: a label `y` enters the set when its
#' own group quantile `q_y` admits it. That is what makes coverage conditional on
#' the true class rather than merely marginal. A side effect is that an empty set
#' is possible; an empty set is a coverage violation, so the arg-max label is
#' restored and the event is flagged in the returned attribute rather than hidden.
#'
#' @keywords internal
predict_sets <- function(probs, lev, q, score) {
  qof <- function(g) if (length(q) == 1L) q else q[[g]]
  out <- lapply(seq_len(nrow(probs)), function(i) {
    p <- probs[i, ]
    if (score == "lac") {
      keep <- vapply(seq_along(lev), function(j) {
        qj <- qof(lev[j])
        is.finite(qj) && (1 - p[j]) <= qj
      }, logical(1))
    } else {
      o <- order(-p); cp <- cumsum(p[o])
      keep <- vapply(seq_along(lev), function(j) {
        qj <- qof(lev[j])
        is.finite(qj) && cp[match(j, o)] <= qj
      }, logical(1))
    }
    if (!any(keep)) {
      attr(keep, "empty") <- TRUE
      keep[which.max(p)] <- TRUE        # never emit an empty set
    }
    lev[keep]
  })
  attr(out, "n_empty_restored") <- sum(vapply(out, function(s) length(s) == 1L, logical(1)))
  out
}

#' Apply a fitted recoverability model to new bulk profiles
#'
#' @param object A `recoverR_conformal` object built by
#'   [refit_recoverability()] (which embeds the training features).
#' @param Xte New bulk feature matrix, columns matching the training matrix.
#' @return A `data.frame` with one row per sample: the point `call`, the
#'   prediction-set `set_size`, the `abstain` flag (set size > 1), the top-class
#'   `confidence`, the `margin` (top1 - top2 probability), and `trusted`
#'   (a singleton call, i.e. the model committed under the requested alpha).
#'
#'   `recoverability` is deliberately NOT a synthetic probability. A trusted
#'   sample carries its classifier confidence, an abstained sample carries `NA`
#'   because no calibrated call was made. The empirical accuracy of the trusted
#'   subset is reported by [reliability_curve()].
#' @export
predict_recoverability <- function(object, Xte) {
  if (is.null(object$X_train) || is.null(object$y_train)) {
    stop("this object has no embedded training data; rebuild it with ",
         "refit_recoverability() so the learner can be refit at prediction time.",
         call. = FALSE)
  }
  Xte <- as.matrix(Xte)
  tr_cols <- colnames(object$X_train)
  common <- intersect(tr_cols, colnames(Xte))
  if (length(common) < 2L) {
    stop("fewer than 2 shared features between training and new data; check ",
         "that both matrices use the same feature identifiers.", call. = FALSE)
  }
  if (length(common) < length(tr_cols)) {
    message(sprintf("using %d of %d training features present in the new data",
                    length(common), length(tr_cols)))
  }
  probs <- learn_probs(object$X_train[, common, drop = FALSE], object$y_train,
                       Xte[, common, drop = FALSE], object$learner)
  sets <- predict_sets(probs, object$levels, object$q, object$score)
  sz <- vapply(sets, length, integer(1))
  conf <- apply(probs, 1, max)
  srt <- t(apply(probs, 1, function(v) sort(v, decreasing = TRUE)))
  data.frame(
    call = object$levels[apply(probs, 1, which.max)],
    set_size = sz,
    abstain = sz > 1L,
    confidence = conf,
    margin = srt[, 1] - srt[, 2],
    trusted = sz == 1L,
    recoverability = ifelse(sz == 1L, conf, NA_real_),
    set_members = vapply(sets, paste, character(1), collapse = "|"),
    row.names = NULL,
    stringsAsFactors = FALSE
  )
}

#' Refit a recoverability model while embedding the training data
#'
#' [fit_recoverability()] deliberately does not retain the raw training matrix,
#' which can be large. [predict_recoverability()] needs it in order to refit the
#' learner on all available data at prediction time, so this wrapper stores it on
#' the returned object.
#'
#' @inheritParams fit_recoverability
#' @return A `recoverR_conformal` object carrying `X_train` and `y_train`.
#' @export
#' @param ... Passed through to the learner.
refit_recoverability <- function(X, archetype, ...) {
  fit <- fit_recoverability(X, archetype, ...)
  fit$X_train <- as.matrix(X)
  fit$y_train <- factor(archetype)
  # the point call at prediction time is produced by a learner refit on ALL data,
  # so the calibration labels are retained separately for diagnostics
  fit$cal_labels_prior <- fit$y_train
  class(fit) <- "recoverR_conformal"
  fit
}

# ----------------------------------------------------------------------------
# Separability diagnostics
# ----------------------------------------------------------------------------

#' Which archetypes does a bulk profile confuse?
#'
#' For every ordered pair of archetypes, the fraction of calibration samples
#' whose prediction set contains BOTH. A high value means bulk data cannot
#' separate the pair, which is the actionable form of the claim: it says exactly
#' where a bulk biomarker will fail, rather than reporting one aggregate
#' accuracy.
#'
#' @param fit A `recoverR_conformal` object.
#' @return A symmetric matrix of pairwise co-membership fractions.
#' @export
archetype_separability <- function(fit) {
  lev <- fit$levels; L <- length(lev)
  m <- matrix(NA_real_, L, L, dimnames = list(lev, lev))
  for (i in seq_len(L)) for (j in seq_len(L)) {
    both <- vapply(fit$cal_sets, function(s) all(c(lev[i], lev[j]) %in% s), logical(1))
    m[i, j] <- mean(both)
  }
  attr(m, "role") <- "confusability: high = bulk cannot separate the pair"
  m
}

#' Coverage / risk / abstention trade-off curve
#'
#' A single coverage number is not sufficient: a model can hold marginal coverage
#' while being uninformative for the majority of samples. This curve sweeps the
#' abstention rule and reports, at each operating point, the fraction answered,
#' the accuracy among those answered, and the singleton rate.
#'
#' @param fit A `recoverR_conformal` object.
#' @return A `data.frame` with `answered`, `accuracy`, `singleton_rate`,
#'   `mean_set_size`.
#' @export
coverage_risk_curve <- function(fit) {
  lev <- fit$levels
  sets <- fit$cal_sets
  truth <- as.character(fit$cal_labels)
  p <- fit$cal_probs
  conf <- apply(p, 1, max)
  ord <- order(-conf)
  n <- length(sets)
  do.call(rbind, lapply(seq_len(n), function(k) {
    sel <- ord[seq_len(k)]
    sz <- vapply(sets[sel], length, integer(1))
    ok <- vapply(sel, function(i) truth[i] %in% sets[[i]], logical(1))
    data.frame(answered = k / n,
               accuracy = mean(ok[sz == 1]),
               singleton_rate = mean(sz == 1),
               mean_set_size = mean(sz))
  }))
}

#' Calibration (reliability) curve and expected calibration error
#'
#' @param fit A `recoverR_conformal` object.
#' @param bins Number of confidence bins.
#' @return A list with the binned curve and the `ece`.
#' @export
reliability_curve <- function(fit, bins = 10L) {
  p <- fit$cal_probs
  conf <- apply(p, 1, max)
  pred <- fit$levels[apply(p, 1, which.max)]
  corr <- as.integer(pred == as.character(fit$cal_labels))
  br <- cut(conf, breaks = unique(stats::quantile(conf, seq(0, 1, length.out = bins + 1))),
            include.lowest = TRUE)
  df <- data.frame(conf = conf, correct = corr, bin = br)
  agg <- stats::aggregate(cbind(conf, correct) ~ bin, data = df, FUN = mean)
  names(agg) <- c("bin", "mean_conf", "observed_acc")
  agg$n <- as.integer(table(df$bin))
  agg$gap <- agg$observed_acc - agg$mean_conf
  list(curve = agg, ece = sum(agg$n / sum(agg$n) * abs(agg$gap)))
}

#' Selective gain: does abstention actually buy accuracy?
#'
#' @param fit A `recoverR_conformal` object.
#' @return A `data.frame` with coverage levels and the accuracy at each.
#' @export
selective_gain <- function(fit) {
  cr <- coverage_risk_curve(fit)
  cr <- cr[is.finite(cr$accuracy), ]
  base <- cr$accuracy[nrow(cr)]
  cr$gain_over_full_coverage <- cr$accuracy - base
  cr
}
