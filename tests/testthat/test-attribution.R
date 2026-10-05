# Tests for the attribution algorithms.
# These verify the AXIOMS, not just that the code runs -- the efficiency axiom in
# particular is what makes the decomposition interpretable, and it is easy to
# break silently by mis-keying the subset lattice.

make_views <- function(n = 80, K = 3, p = 6, seed = 1, signal_view = 1) {
  set.seed(seed)
  X <- lapply(seq_len(K), function(k) {
    m <- matrix(stats::rnorm(n * p), n, p)
    if (k == signal_view) m[, 1] <- m[, 1] + 1.5 * ((seq_len(n) %% 2) * 2 - 1)
    m
  })
  names(X) <- paste0("view", seq_len(K))
  y <- factor(ifelse(X[[signal_view]][, 1] > stats::median(X[[signal_view]][, 1]), "b", "a"),
              levels = c("a", "b"))
  list(X = X, y = y)
}

test_that("Shapley values satisfy the efficiency axiom", {
  d <- make_views(n = 70, K = 3, seed = 11)
  sh <- shapley_attribution(d$X, d$y, family = "binomial", nfolds = 3L, seed = 1L)
  # sum of Shapley values must equal V(V) - V(empty), with V(empty) := 0
  expect_equal(sum(sh$shapley), unique(sh$V_full), tolerance = 1e-8)
})

test_that("the permutation null separates real signal from noise", {
  set.seed(4)
  n <- 240
  signal <- rnorm(n)
  # BOTH views are 4 features wide, but every column of `real` carries the signal
  # (signal + independent noise). A "real" view that is 3/4 pure noise makes the
  # two views statistically comparable and the test meaningless -- which is how an
  # earlier version of this test was written.
  X <- list(
    real  = matrix(signal, n, 4) + matrix(stats::rnorm(n * 4, sd = 0.6), n, 4),
    noise = matrix(stats::rnorm(n * 4), n, 4)
  )
  y <- factor(ifelse(signal > median(signal), "b", "a"), levels = c("a", "b"))

  sh <- shapley_attribution(X, y, family = "binomial", nfolds = 3L, seed = 1L)
  expect_true(all(is.finite(sh$shapley)))
  # a noise view can carry a non-zero Shapley value: random columns change the
  # elastic-net fit and can act like noise-injection regularisation. What must
  # hold is a dominance gap, not a zero.
  expect_gt(sh$shapley[sh$view == "real"], 2.5 * sh$shapley[sh$view == "noise"])

  pr <- permutation_null(X, y, view = "real",  family = "binomial", B = 49L, nfolds = 3L, seed = 1L)
  pn <- permutation_null(X, y, view = "noise", family = "binomial", B = 49L, nfolds = 3L, seed = 1L)
  expect_true(all(is.finite(c(pr$p_value, pn$p_value))))
  expect_true(pr$p_value >= 0 && pr$p_value <= 1)

  # THE DISCRIMINATOR: for a real signal the observed marginal contribution sits
  # far ABOVE its own permutation null. For a noise view the permutation is
  # distribution-preserving, so the observed value sits essentially AT the null
  # mean and the p-value is uniform under H0 by construction. Asserting
  # "p > 0.05" for one seed would therefore fail ~5% of the time by design; the
  # gap-to-null is the stable, correct statistic.
  gap_real  <- pr$delta_observed - pr$null_mean
  gap_noise <- abs(pn$delta_observed - pn$null_mean)
  expect_gt(gap_real, 5 * gap_noise)
  expect_gt(gap_real, 0)
})

test_that("Shapley attribution does not depend on the view listing order", {
  d <- make_views(n = 70, K = 3, seed = 5)
  s1 <- shapley_attribution(d$X, d$y, family = "binomial", nfolds = 3L, seed = 1L)
  s2 <- shapley_attribution(rev(d$X), d$y, family = "binomial", nfolds = 3L, seed = 1L)
  m <- merge(s1[, c("view", "shapley")], s2[, c("view", "shapley")], by = "view")
  # exact order-independence holds mathematically; the residual is float-level
  # noise from glmnet's internal lambda path, so the tolerance is stated rather
  # than pretended away. The ranking must be identical.
  expect_equal(m$shapley.x, m$shapley.y, tolerance = 1e-2)
  expect_identical(order(-m$shapley.x), order(-m$shapley.y))
})

test_that("subset lattice of the wrong size is refused rather than mis-keyed", {
  d <- make_views(n = 40, K = 2, seed = 7)
  names(d$X) <- c("dup", "dup")     # duplicate names cannot form a lattice
  expect_error(shapley_attribution(d$X, d$y, family = "binomial", nfolds = 2L),
               "lattice|unique|subscript")
})

test_that("shapley_attribution refuses to run above max_views", {
  d <- make_views(n = 60, K = 3, seed = 9)
  expect_error(shapley_attribution(d$X, d$y, family = "binomial", max_views = 2L),
               "2\\^K|max_views")
})

test_that("permutation null is valid and sensitive", {
  d <- make_views(n = 60, K = 2, seed = 13, signal_view = 1)
  pn <- permutation_null(d$X, d$y, view = "view1", family = "binomial",
                         B = 29L, nfolds = 3L, seed = 1L)
  expect_true(pn$p_value >= 0 && pn$p_value <= 1)
  expect_gt(pn$delta_observed, pn$null_mean)          # the signal view beats its null
})

test_that("leave-one-out flags an actively harmful view", {
  set.seed(21)
  n <- 90
  sig <- rnorm(n)
  y <- factor(ifelse(sig > median(sig), "b", "a"), levels = c("a", "b"))
  X <- list(
    informative = matrix(sig, ncol = 1),
    pure_noise  = matrix(stats::rnorm(n * 20), n, 20)
  )
  ab <- ablate_views(X, y, family = "binomial", nfolds = 3L, seed = 1L)
  expect_gt(ab$leave_one_out$loss_when_removed[ab$leave_one_out$view == "informative"], 0)
})

test_that("performance functionals agree with reference implementations", {
  set.seed(2)
  y <- rep(c(0, 1), each = 20)
  s <- c(rnorm(20, -1), rnorm(20, 1))
  a <- recoverR:::perf_auc(y, s)
  if (requireNamespace("pROC", quietly = TRUE)) {
    expect_equal(a, as.numeric(pROC::auc(pROC::roc(y, s, quiet = TRUE))), tolerance = 1e-9)
  }
  expect_equal(recoverR:::perf_r2(1:8, 1:8), 1)   # n >= 5 is enforced
})

test_that("perf_auc survives a perfect classifier (where pROC::roc.test fails)", {
  y <- rep(c(0, 1), each = 10)
  s <- c(rep(0, 10), rep(1, 10))
  expect_equal(recoverR:::perf_auc(y, s), 1)
  r <- auc_safe(s, factor(y))
  expect_equal(r$auc, 1)
  expect_true(is.finite(r$p))          # pROC::roc.test() would give NA here
})

test_that("fit_predict returns ONE prediction per held-out row when the inner CV cannot be formed", {
  # Regression for the defect measured on 2026-09-27: when EVERY inner cv.glmnet() errored,
  # fit_predict() fell back to glmnet() with NO lambda, so predict() returned the whole
  # n_test x 100 path as a matrix; as.numeric() flattened it and cv_perf()'s `oof[te] <- ...`
  # RECYCLED the first n_test values into the fold (500 values for 5 held-out rows at n = 16),
  # scoring the fold with an undeclared model behind a single warning.
  # The trigger is deterministic: 2 minority-class samples with 5 inner folds makes every
  # cv.glmnet() error ("one multinomial or binomial class has 1 or 0 observations") while a
  # plain glmnet() fit on the same training set still succeeds -- exactly the regime the
  # defect fired in.
  set.seed(3)
  n <- 12L; p <- 4L
  Xtr <- matrix(stats::rnorm(n * p), n, p)
  ytr <- factor(c("b", "b", rep("a", n - 2L)), levels = c("a", "b"))
  Xte <- matrix(stats::rnorm(7L * p), 7L, p)

  # The trigger is deliberately degenerate, so glmnet warns repeatedly ("one multinomial or
  # binomial class has fewer than 8 observations; dangerous ground") - expected, and suppressed
  # here so `R CMD check`'s test output stays clean and a future reader does not read those
  # warnings as evidence of a problem. What is being tested is the LENGTH of the return value.
  pr <- suppressWarnings(fit_predict(Xtr, ytr, Xte, "binomial", nfolds_inner = 5L))
  expect_length(pr, 7L)                # the old code returned 7 * 100 = 700 values here
  expect_true(all(is.finite(pr)))
})

test_that("cv_perf STOPs on a wrong-length prediction vector instead of recycling it", {
  # Defence in depth for the same defect: `oof[te] <- <wrong length>` recycles silently in R.
  # Assert the invariant directly, so a future learner that returns a matrix or a path again
  # fails loudly at the assignment site rather than quietly scoring the fold.
  local_mocked_bindings(fit_predict = function(...) rep(0.5, 1000L))
  set.seed(5)
  n <- 30L
  X <- list(v1 = matrix(stats::rnorm(n * 3L), n, 3L))
  y <- factor(rep(c("a", "b"), length.out = n), levels = c("a", "b"))
  expect_error(cv_perf(X, y, family = "binomial", nfolds = 3L), "predictions for")
})
