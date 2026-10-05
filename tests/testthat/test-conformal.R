# Tests for the conformal recoverability layer.
# The decisive property is the finite-sample coverage guarantee: on exchangeable
# data the empirical coverage must sit at or just above 1 - alpha. If this test
# fails, every other number the layer reports is meaningless.

make_bulk_archetype <- function(n = 400, seed = 1, sep = 1.6, n_classes = 3) {
  set.seed(seed)
  lab <- rep(letters[seq_len(n_classes)], length.out = n)
  cen <- matrix(0, n_classes, 6)
  # separate classes along distinct axes so the problem is learnable
  for (i in seq_len(n_classes)) cen[i, i] <- sep
  X <- cen[match(lab, letters[seq_len(n_classes)]), , drop = FALSE] +
    matrix(stats::rnorm(n * 6), n, 6)
  colnames(X) <- paste0("f", seq_len(6))
  list(X = X, archetype = factor(lab))
}

test_that("split conformal attains its finite-sample coverage target", {
  d <- make_bulk_archetype(n = 500, seed = 3, sep = 1.6)
  fit <- fit_recoverability(d$X, d$archetype, alpha = 0.10, mondrian = FALSE,
                            cal_prop = 0.5, seed = 1L)
  expect_gte(fit$cal_coverage, 1 - fit$alpha - 0.06)
})

test_that("Mondrian calibration gives class-conditional coverage", {
  d <- make_bulk_archetype(n = 600, seed = 4, sep = 1.4)
  fit <- fit_recoverability(d$X, d$archetype, alpha = 0.10, mondrian = TRUE,
                            cal_prop = 0.5, seed = 1L)
  expect_true(all(fit$cal_coverage_by_class >= 1 - fit$alpha - 0.12))
})

test_that("a smaller alpha produces larger prediction sets", {
  d <- make_bulk_archetype(n = 400, seed = 6, sep = 1.2)
  loose <- fit_recoverability(d$X, d$archetype, alpha = 0.20, seed = 1L)
  tight <- fit_recoverability(d$X, d$archetype, alpha = 0.05, seed = 1L)
  m <- function(f) mean(vapply(f$cal_sets, length, integer(1)))
  expect_gt(m(tight), m(loose))
})

test_that("abstention rises as the classes become indistinguishable", {
  easy <- make_bulk_archetype(n = 400, seed = 8, sep = 2.5)
  hard <- make_bulk_archetype(n = 400, seed = 8, sep = 0.0)
  f_easy <- fit_recoverability(easy$X, easy$archetype, alpha = 0.10, seed = 1L)
  f_hard <- fit_recoverability(hard$X, hard$archetype, alpha = 0.10, seed = 1L)
  abst <- function(f) mean(vapply(f$cal_sets, length, integer(1)) > 1)
  expect_gt(abst(f_hard), abst(f_easy))
})

test_that("separability is symmetric with unit diagonal-ish structure", {
  d <- make_bulk_archetype(n = 400, seed = 10, sep = 0.8)
  fit <- fit_recoverability(d$X, d$archetype, alpha = 0.10, seed = 1L)
  m <- archetype_separability(fit)
  expect_equal(dim(m), c(3L, 3L))
  expect_true(isSymmetric(unname(m), tol = 1e-8))
  expect_true(all(m >= 0 & m <= 1))
})

test_that("coverage-risk and reliability curves are well formed", {
  d <- make_bulk_archetype(n = 300, seed = 12, sep = 1.5)
  fit <- fit_recoverability(d$X, d$archetype, alpha = 0.10, seed = 1L)
  cr <- coverage_risk_curve(fit)
  expect_true(all(c("answered", "accuracy", "singleton_rate") %in% names(cr)))
  expect_true(all(cr$answered >= 0 & cr$answered <= 1))
  rc <- reliability_curve(fit)
  expect_true(is.finite(rc$ece) && rc$ece >= 0 && rc$ece <= 1)
})

test_that("prediction sets are never empty", {
  d <- make_bulk_archetype(n = 300, seed = 14, sep = 0.3)
  fit <- fit_recoverability(d$X, d$archetype, alpha = 0.05, mondrian = TRUE, seed = 1L)
  sizes <- vapply(fit$cal_sets, length, integer(1))
  expect_true(all(sizes >= 1L))
})

test_that("predict_recoverability refuses an object without training data", {
  d <- make_bulk_archetype(n = 200, seed = 16)
  fit <- fit_recoverability(d$X, d$archetype, seed = 1L)
  expect_error(predict_recoverability(fit, d$X), "no embedded training data")
})

test_that("refit + predict returns the documented columns and abstains sensibly", {
  d <- make_bulk_archetype(n = 400, seed = 18, sep = 1.2)
  fit <- refit_recoverability(d$X, d$archetype, alpha = 0.10, seed = 1L)
  pr <- predict_recoverability(fit, d$X)
  expect_true(all(c("call", "set_size", "abstain", "confidence", "trusted",
                    "recoverability") %in% names(pr)))
  expect_equal(nrow(pr), nrow(d$X))
  expect_true(all(pr$set_size >= 1L))
  # recoverability is NA exactly for abstained samples by design
  expect_true(all(is.na(pr$recoverability[pr$abstain])))
})

# ---------------------------------------------------------------------------
# The learner guard's CRITERION (added 0.1.2 after the positive control that found it).
#
# The original guard fired only at p >= n. That is not sufficient: at p = 2000,
# n = 2235 (p < n) MASS::lda() inverts a 2000 x 2000 pooled covariance and returns
# SATURATED posteriors (confidence ~1.0 at chance accuracy), and the conformal layer
# then builds its per-class quantiles from them and admits every class for every
# sample. Measured on a real region-scale matrix with a planted, zero-noise signal:
# lda -> accuracy 0.725, ECE 0.274, abstention 0.996; glmnet -> 0.954 / 0.029 / 0.000.
# These two tests pin the criterion in BOTH directions: it must fire just below n,
# and it must NOT fire where lda is genuinely appropriate.
# ---------------------------------------------------------------------------

test_that("the learner guard fires when p is close to n, not only at p >= n", {
  set.seed(2)
  n <- 120; p <- 100                 # p < n, but only 1.2 samples per feature
  X <- matrix(stats::rnorm(n * p), n, p, dimnames = list(NULL, paste0("p", seq_len(p))))
  y <- factor(rep(c("a", "b"), each = n / 2))
  expect_warning(recoverR:::learn_probs(X, y, X, learner = "lda"),
                 "ill-conditioned")
})

test_that("the learner guard does NOT fire where lda is well conditioned (p << n)", {
  d <- make_bulk_archetype(n = 400, seed = 20)      # p = 6, n = 400 -> 30 < 400
  expect_silent(recoverR:::learn_probs(d$X, d$archetype, d$X, learner = "lda"))
})

test_that("the default learner recovers a PLANTED signal at p close to n", {
  # The positive control, in miniature: a label built from the first 10 features with
  # no noise must come back with high accuracy. Before 0.1.2 this returned ~0.73 under
  # the default learner; the point of the guard is that the default is safe.
  set.seed(3)
  n <- 240; p <- 240                 # p == n, the worst case the old guard caught...
  X <- matrix(stats::rnorm(n * p), n, p, dimnames = list(NULL, paste0("p", seq_len(p))))
  s <- rowMeans(scale(X[, 1:10, drop = FALSE]))
  y <- factor(ifelse(s > stats::median(s), "hi", "lo"))
  prob <- suppressWarnings(recoverR:::learn_probs(X, y, X))
  acc <- mean(colnames(prob)[apply(prob, 1, which.max)] == as.character(y))
  expect_gt(acc, 0.90)
})
