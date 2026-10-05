## ============================================================================
## make_tutorial.R — regenerate every figure and every printed table quoted in
##                   TUTORIAL.md and README.md.
##
##   Rscript tutorial/make_tutorial.R          # from the package root
##
## What it does
##   * builds small SYNTHETIC designs with a known planted signal (no study data
##     is used or needed — the tutorial must be reproducible by a stranger);
##   * runs the package's three modules on them (attribution, conformal
##     recoverability, sampling adequacy) plus the guardrails;
##   * writes every figure to man/figures/ as a VECTOR PDF **and** a JPEG, by the
##     package's own save_plot() — the JPEG is what Markdown embeds, the PDF is
##     the editable/publication copy. Both are asserted to exist afterwards.
##
## Output is deterministic: every stochastic step is seeded, and the seeds are
## printed below so a reader can reproduce a given panel exactly.
## ============================================================================

suppressPackageStartupMessages(library(recoverR))

FIGDIR <- file.path("man", "figures")
dir.create(FIGDIR, showWarnings = FALSE, recursive = TRUE)

DPI <- 300L
figs <- character(0)

## write a figure and assert the PDF+JPEG pair really landed (a half-written
## pair is the classic silent failure — the package's own save_plot() warns,
## this makes it fatal for the tutorial build)
fig <- function(plot, name, width, height) {
  stem <- file.path(FIGDIR, name)
  save_plot(plot, stem, width = width, height = height, dpi = DPI, quiet = TRUE)
  pair <- paste0(stem, c(".pdf", ".jpeg"))
  missing <- pair[!file.exists(pair)]
  if (length(missing)) stop("figure pair incomplete: ", paste(missing, collapse = ", "))
  if (!is_vector_pdf(paste0(stem, ".pdf"))) stop("not vector: ", stem, ".pdf")
  figs <<- c(figs, basename(pair))
  cat(sprintf("  figure ok (pdf+jpeg, vector verified): %s\n", basename(stem)))
}

banner <- function(x) cat("\n", strrep("=", 78), "\n## ", x, "\n", strrep("=", 78), "\n", sep = "")
tab <- function(x, row.names = FALSE, digits = 4) {
  x <- as.data.frame(x)
  num <- vapply(x, is.numeric, logical(1))
  x[num] <- lapply(x[num], function(v) ifelse(is.na(v), NA, signif(v, digits)))
  print(x, row.names = row.names)
}

cat("recoverR", as.character(packageVersion("recoverR")),
    "| R", paste(R.version$major, R.version$minor, sep = "."), "\n")

## ============================================================================
## PART 1 — Layer attribution: does adding a view improve OUT-OF-SAMPLE
##          prediction, and is an added view actively harmful?
##
## Design: 120 samples, four views of 30 features. v1 carries the outcome
## (signal 1.0), v2 carries a weaker version of it (0.6), v3 is a NOISY COPY of
## v1 (signal 0) — a redundant layer, the case a bare variance share gets wrong
## — and v4 is pure noise.
## ============================================================================
banner("PART 1 — attribution: 4 synthetic views, planted signal in v1/v2, v3 redundant")
set.seed(20261005)
n <- 120L
y <- rnorm(n)
make_view <- function(p, signal, sd_noise = 1) {
  m <- matrix(rnorm(n * p), n, p)
  if (signal > 0) m[, 1] <- signal * y + rnorm(n)
  m + matrix(rnorm(n * p, sd = sd_noise), n, p)
}
X <- list(
  v1_signal    = make_view(30, 1.0),   # planted
  v2_weaker    = make_view(30, 0.6),   # planted, weaker
  v3_redundant = {                     # a NOISY COPY of a planted view: it
    z <- make_view(30, 1.0)            # shares v1's signal but carries none
    z + matrix(rnorm(n * 30, sd = 1.5), n, 30)   # of its own, i.e. a redundant
  },                                   # layer — the case a bare variance share
  v4_noise     = make_view(30, 0.0)    # gets wrong.  v4 is pure noise.
)
cat("\nview dimensions (samples x features):\n"); print(sapply(X, dim))

res <- attrib_views(X, y, family = "gaussian", nfolds = 5L, seed = 1L,
                    n_perm = 49L, ncores = 4L, verbose = FALSE)
cat("\n-- res (Shapley + unique contribution) --\n"); print(res)

cat("\n-- dimension-matched permutation null --\n")
null_tab <- as.data.frame(res$null)[, c("view", "delta_observed", "null_mean", "null_q95", "p_value")]
tab(null_tab)

cat("\n-- leave-one-out (out-of-sample performance without each view) --\n")
tab(res$ablation$leave_one_out)

cat("\n-- forward selection (views added greedily) --\n")
tab(res$ablation$forward)
cat(sprintf("\nV(full) = %.4f\n", res$ablation$V_full))

fig(plot_shapley(res),             "tut-01-shapley",        6.5, 4.0)
fig(plot_attribution_plate(res),   "tut-02-plate",         10.0, 7.5)   # wide enough that
fig(plot_leave_one_out(res),       "tut-03-leave-one-out",  6.5, 3.5)   # the two panel titles
fig(plot_forward_selection(res),   "tut-04-forward",        6.5, 3.5)

## ============================================================================
## PART 2 — Bulk -> spatial recoverability with abstention.
##
## The predictor sees only the AGGREGATED profile of a region; the label is the
## spatial archetype. Nothing about the label enters the feature vector, so the
## question is not circular: it measures what aggregation destroys.
## ============================================================================
banner("PART 2 — conformal recoverability: 300 regions, 3 archetypes, 2 informative features")
set.seed(20261005)
nr <- 300L
lab <- sample(c("inflamed", "excluded", "desert"), nr, replace = TRUE)
shift <- c(inflamed = 1.2, excluded = -0.2, desert = -1.4)[lab]
Xr <- matrix(rnorm(nr * 40), nr, 40)
Xr[, 1] <- shift + rnorm(nr, sd = 0.9)
Xr[, 2] <- shift + rnorm(nr, sd = 1.1)
colnames(Xr) <- sprintf("g%02d", seq_len(ncol(Xr)))   # features MUST be named
cat("\narchetype counts:\n"); print(table(lab))

fit <- refit_recoverability(Xr, factor(lab), alpha = 0.10, mondrian = TRUE,
                            score = "aps", cal_prop = 0.5, seed = 1)
cat(sprintf("\ncalibration coverage: %.3f (target %.2f)\n", fit$cal_coverage, 0.90))
cat("coverage by class:\n"); print(round(fit$cal_coverage_by_class, 3))

pred <- predict_recoverability(fit, Xr)
cat(sprintf("\nabstention: %d of %d regions (%.1f%%)\n",
            sum(pred$abstain), nrow(pred), 100 * mean(pred$abstain)))
cat("\nfirst six prediction sets:\n")
tab(as.data.frame(pred)[1:6, ])

cat("\n-- pairwise confusability matrix (high = this pair cannot be told apart) --\n")
sep <- archetype_separability(fit)
print(round(sep, 3))
cat("\noff-diagonal range: ",
    sprintf("%.3f to %.3f\n", min(sep[upper.tri(sep)]), max(sep[upper.tri(sep)])))

rc <- reliability_curve(fit)
cat(sprintf("\nreliability ECE: %.3f\n", rc$ece))
cr <- coverage_risk_curve(fit)
cat(sprintf("risk-coverage curve: %d operating points; answered %.1f%% to %.1f%%\n",
            nrow(cr), 100 * min(cr$answered), 100 * max(cr$answered)))
sg <- selective_gain(fit)
tab(sg)

fig(plot_recoverability(pred),   "tut-05-recoverability", 6.5, 4.5)
fig(plot_reliability(fit),       "tut-06-reliability",    5.5, 5.0)
fig(plot_coverage_risk(fit),     "tut-07-coverage-risk",  6.0, 4.0)
fig(plot_separability(fit),      "tut-08-separability",   6.5, 5.5)

## ============================================================================
## PART 3 — Sampling adequacy: how much tissue, and how many genes, does the
##          question need?
##
## The class scale must be LARGE relative to the largest region sampled,
## otherwise a big region averages the two archetypes together and accuracy
## FALLS with diameter instead of saturating (the Hill fit then correctly
## rejects: "a parameter is pinned at its bound"). So this design uses a
## 40x40 lattice at pitch 50 with field wavelengths of 6000/7000 units
## against sampled diameters of 200-1200 — that ratio is what makes the
## plateau identifiable.
## ============================================================================
banner("PART 3 — sampling adequacy: 1600 synthetic spots, 40x40 lattice, pitch 50")
set.seed(11)
nside <- 40L; pitch <- 50
n_spot <- nside * nside
spots  <- sprintf("s%05d", seq_len(n_spot))
coords <- data.frame(x = rep(seq_len(nside), each = nside) * pitch,
                     y = rep(seq_len(nside), times = nside) * pitch,
                     barcode = spots, row.names = spots)
ext    <- nside * pitch
field  <- sin(2 * pi * (coords$x - ext / 2) / 6000) + 0.7 * cos(2 * pi * (coords$y - ext / 2) / 7000)
labels <- setNames(ifelse(field > 0, "inflamed", "desert"), spots)
cnt <- matrix(rpois(n_spot * 200, 5), nrow = n_spot,
              dimnames = list(spots, paste0("GENE", seq_len(200))))
cnt[,   1:25] <- cnt[,   1:25] + pmax(0,  field) * 6
cnt[, 26:50] <- cnt[, 26:50] + pmax(0, -field) * 6
cat("\nspot labels (the ground truth the aggregate must recover):\n"); print(table(labels))

sa <- sampling_adequacy(cnt, coords, labels = labels,
                        diameters = c(200, 300, 500, 800, 1200),
                        panel_sizes = c(50, 200), n_regions = 60L,
                        target = 0.95, seed = 1L)
cat("\n-- the surface: recoverable signal vs region diameter --\n")
tab(sa$surface[, c("diameter_mm", "panel_genes", "performance")])
cat("\n-- adequacy (per panel size) --\n")
tab(sa$adequacy[, c("panel_genes", "plateau", "half_saturation_mm", "adequacy_mm", "extrapolated")])

## The flag, demonstrated: the SAME endpoint swept over diameters that are all
## far too coarse for the class scale. Nothing here is a bug - the tool is
## refusing to report a number it cannot identify.
cat("\n-- the same endpoint over too-coarse diameters (the extrapolated case) --\n")
sa_coarse <- tryCatch(
  sampling_adequacy(cnt, coords, labels = labels,
                    diameters = c(400, 800, 1600, 2000, 2400),
                    panel_sizes = c(50, 200), n_regions = 60L, target = 0.95, seed = 1L),
  error = function(e) { cat("  (sweep refused): ", conditionMessage(e), "\n", sep = ""); NULL })
if (!is.null(sa_coarse)) {
  tab(sa_coarse$surface[, c("diameter_mm", "panel_genes", "performance")])
  tab(sa_coarse$adequacy[, c("panel_genes", "plateau", "adequacy_mm", "extrapolated")])
}

fig(plot_sampling_curve(sa$surface, adequacy = 0.95), "tut-09-sampling-curve", 7.5, 4.5)

## ============================================================================
## GUARDRAILS — the failures that return a plausible number instead of an error
## ============================================================================
banner("GUARDRAILS")
msg <- tryCatch({ assert_no_na(matrix(c(1, NA, 3))); "NO ERROR (bad)" },
                error = function(e) conditionMessage(e))
cat("assert_no_na() on a matrix with an NA:\n  ", msg, "\n", sep = "")
tmp <- tempfile(fileext = ".pdf")
grDevices::pdf(tmp); plot(1:10, 1:10); invisible(grDevices::dev.off())
cat(sprintf("is_vector_pdf() on a freshly drawn PDF: %s\n", is_vector_pdf(tmp)))

## ============================================================================
banner("DONE")
cat("figures written to ", FIGDIR, " (", length(figs), " files):\n", sep = "")
cat(paste0("  ", figs, collapse = "\n"), "\n")
cat("\nRe-run: Rscript tutorial/make_tutorial.R\n")
