# ============================================================================
# recoverR - sampling.R
#
# Sampling adequacy: turning a tissue/panel constraint into a design rule.
#
# THE PROBLEM THIS SOLVES
#   Spatially predictive signal is known empirically to collapse when a specimen
#   is small (roughly sub-3 mm regions) and when the gene panel is short (tens of
#   genes). That constraint is usually reported once, for one study, as an
#   afterthought. Nobody ships it as a design tool, so the same mistake is
#   repeated: a trial is planned around a biomarker that the planned sampling
#   cannot support.
#
# THE ANSWER IMPLEMENTED HERE
#   Simulate the sampling geometry directly. Collapse spatial measurements into
#   pseudo-bulk at a grid of region diameters and panel sizes, score the endpoint
#   at each vertex, and fit a saturating response surface. The adequacy threshold
#   is then read off that surface as an explicit quantity.
# ============================================================================


#' Collapse spatial spots into pseudo-bulk samples over circular regions
#'
#' A region of diameter `d` centred at a random valid location yields one
#' pseudo-bulk profile by summing counts over the spots it contains. Sweeping `d`
#' therefore generates the paired (bulk, spatial) data that a recoverability
#' model needs, without any new specimens.
#'
#' @param counts Matrix of counts, spots x genes, rownames matched to `coords`.
#' @param coords Data frame with numeric `x` and `y` (spot coordinates in the
#'   same units as `diameter`).
#' @param diameter Region diameter in coordinate units.
#' @param n_regions Number of regions to sample.
#' @param min_spots Discard regions with fewer than this many spots.
#' @param aggregate `"sum"` or `"mean"` over spots.
#' @param seed RNG seed.
#' @return A list with `bulk` (regions x genes), `spot_index` (list of spot
#'   rownames per region) and `n_spots`.
#' @export
pseudobulk_regions <- function(counts, coords, diameter, n_regions = 200L,
                               min_spots = 3L, aggregate = c("sum", "mean"),
                               seed = 1L) {
  aggregate <- match.arg(aggregate)
  stopifnot(all(c("x", "y") %in% names(coords)))
  coords <- as.data.frame(coords)
  if (!is.null(rownames(coords)) && is.null(rownames(counts))) {
    rownames(counts) <- rownames(coords)
  }
  xy <- as.matrix(coords[, c("x", "y")])
  r <- diameter / 2
  set.seed(seed)

  idx <- vector("list", n_regions)
  made <- 0L
  attempts <- 0L
  while (made < n_regions && attempts < n_regions * 50L) {
    attempts <- attempts + 1L
    ctr <- xy[sample(nrow(xy), 1L), ]
    dd <- sqrt(rowSums(sweep(xy, 2, ctr)^2))
    within <- which(dd <= r)
    if (length(within) < min_spots) next
    made <- made + 1L
    idx[[made]] <- within
  }
  idx <- idx[seq_len(made)]
  if (!length(idx)) {
    stop("no region of diameter ", diameter, " contained >= ", min_spots,
         " spots; the requested diameter is smaller than the spot spacing.",
         call. = FALSE)
  }

  # Keep the matrix SPARSE. A Visium section is ~20k genes x ~3k spots; densifying
  # it costs ~500 MB per sample and 16 samples would not fit comfortably. The row
  # subsets here are small, so a sparse colSums is both faster and lighter.
  sparse_in <- inherits(counts, "sparseMatrix")
  if (sparse_in && !requireNamespace("Matrix", quietly = TRUE)) {
    stop("sparse input needs the 'Matrix' package.", call. = FALSE)
  }
  bulk <- t(vapply(idx, function(ix) {
    sub <- counts[ix, , drop = FALSE]
    if (aggregate == "sum") {
      if (inherits(sub, "sparseMatrix")) Matrix::colSums(sub) else colSums(sub)
    } else {
      if (inherits(sub, "sparseMatrix")) Matrix::colMeans(sub) else colMeans(sub)
    }
  }, numeric(ncol(counts))))
  colnames(bulk) <- colnames(counts)
  rownames(bulk) <- paste0("region_", seq_along(idx))
  list(bulk = bulk, spot_index = lapply(idx, function(i) rownames(counts)[i]),
       n_spots = lengths(idx), diameter = diameter)
}


#' Hill saturation model
#'
#' \deqn{P(d) = P_\infty \frac{d^{h}}{K^{h} + d^{h}}}
#' fitted by nonlinear least squares. `K` is the half-saturation diameter and `h`
#' the cooperativity, so the model accommodates a genuine threshold (large `h`)
#' as well as a gradual rise (`h ~ 1`), which a single exponential cannot.
#'
#' The fit is VALIDATED against physical bounds before being returned. An
#' unconstrained NLS on a nearly flat or non-monotone curve can wander to a
#' "solution" with a half-saturation of tens of millions of millimetres, which
#' then propagates into the reported adequacy diameter as a number that looks
#' like a result but is nonsense. A rejected fit returns `NULL` and the caller
#' reports `NA`, so the failure is visible rather than plausible.
#'
#' @param d Diameters.
#' @param p Performance at each diameter.
#' @param start Optional starting values.
#' @param perf_range Admissible range for the plateau. Defaults to `c(0, 1)`,
#'   correct for accuracy/AUC; widen for unbounded scores.
#' @keywords internal
fit_hill <- function(d, p, start = NULL, perf_range = c(0, 1)) {
  ok <- is.finite(d) & is.finite(p) & d > 0
  d <- d[ok]; p <- p[ok]
  if (length(d) < 4L) return(NULL)
  # a curve with no variation carries no saturation information
  if (length(unique(round(p, 6))) < 3L) return(NULL)
  # RESCALE the diameter axis before fitting. With raw diameters in the hundreds
  # to thousands, `d^h` overflows to Inf as the optimiser explores large h
  # (4000^50 ~ 1e180), and nls dies with "Missing value or an infinity produced
  # when evaluating the model" -- so the whole fit silently failed on perfectly
  # well-behaved curves. Working in units of the maximum diameter keeps the base
  # at or below 1 and the model numerically stable; K is converted back after.
  dref <- max(d)
  ds <- d / dref
  if (is.null(start)) {
    start <- list(Pinf = max(p, na.rm = TRUE), Ks = stats::median(ds), h = 2)
  } else if (!is.null(start$K)) {
    start <- list(Pinf = start$Pinf, Ks = start$K / dref, h = start$h)
  }
  # BOUND the fit. An unbounded Hill fit is numerically unstable from both sides:
  # with raw diameters d^h overflows, and after rescaling both Ks^h and ds^h
  # UNDERFLOW to zero for large h, giving 0/0 ("Missing value or an infinity
  # produced when evaluating the model"). Both failures were observed. Bounding h
  # to [0.5, 10] removes the instability and is also the honest choice: a
  # cooperativity above ~10 is not interpretable for a sampling-saturation curve.
  # `algorithm = "port"` is what actually enforces the bounds.
  lo <- c(Pinf = max(1e-6, perf_range[1] + 1e-6), Ks = 1e-4, h = 0.5)
  hi <- c(Pinf = perf_range[2], Ks = 1, h = 10)
  start$Pinf <- min(max(start$Pinf, lo["Pinf"] * 1.01), hi["Pinf"] * 0.99)
  start$Ks   <- min(max(start$Ks, lo["Ks"] * 1.01), hi["Ks"] * 0.99)
  start$h    <- min(max(start$h, lo["h"] * 1.01), hi["h"] * 0.99)
  m <- tryCatch(
    stats::nls(p ~ Pinf * ds^h / (Ks^h + ds^h), start = start,
               algorithm = "port", lower = lo, upper = hi,
               control = stats::nls.control(maxiter = 500, warnOnly = TRUE)),
    error = function(e) NULL)
  if (is.null(m)) return(NULL)
  cf <- stats::coef(m)
  names(cf)[names(cf) == "Ks"] <- "K"
  cf["K"] <- cf["K"] * dref          # back to the original diameter units
  # --- bound check: reject rather than report an unphysical fit --------------
  bad <- c(
    if (!is.finite(cf["Pinf"]) || cf["Pinf"] <= perf_range[1] ||
        cf["Pinf"] > perf_range[2]) "plateau outside the admissible range",
    if (!is.finite(cf["K"]) || cf["K"] <= 0 || cf["K"] >= max(d) * 100)
      "half-saturation far outside the sampled diameter range",
    if (!is.finite(cf["h"]) || cf["h"] <= 0 || cf["h"] > 50)
      "cooperativity out of range"
  )
  if (length(bad)) {
    message("  Hill fit rejected: ", paste(bad, collapse = "; "),
            " (Pinf = ", signif(cf["Pinf"], 4), ", K = ", signif(cf["K"], 4),
            ", h = ", signif(cf["h"], 4), ")")
    return(NULL)
  }
  dmax <- max(d)
  # Reject a fit whose optimiser ran INTO a bound: h at 0.5 or 10, or Ks at 1e-4
  # or 1, is not an estimated cooperativity/affinity, it is the wall. A pure-noise
  # curve does exactly this (measured 2026-09-23: Pinf 0.733, h 0.500 pinned at the
  # lower bound, K = 1.5), and reporting it as an adequacy threshold is the same
  # class of silent-nonsense output as feeding log2(x+1) to DESeq2.
  tol <- 1e-6
  at_bound <- function(x, lo_v, hi_v) is.finite(x) && (abs(x - lo_v) < tol || abs(x - hi_v) < tol)
  if (at_bound(unname(cf["h"]), lo["h"], hi["h"]) ||
      at_bound(unname(cf["K"] / dref), lo["Ks"], hi["Ks"])) {
    message("  Hill fit rejected: a parameter is pinned at its bound - the optimiser ",
            "ran to the wall, so the curve is not identified",
            " (Pinf = ", signif(cf["Pinf"], 4), ", K = ", signif(cf["K"], 4),
            ", h = ", signif(cf["h"], 4), ")")
    return(NULL)
  }
  list(coef = cf, model = m, d_range = c(min(d), dmax),
       # A Hill curve fitted to diameters that never reach its plateau cannot
       # identify K or h, and d* then falls OUTSIDE the sampled range -- an
       # extrapolation, not a measurement. Flagged so a caller can report "the
       # required diameter exceeds what we measured" instead of quoting a number.
       d_star = function(target = 0.95) {
         v <- unname(cf["K"] * (target / (1 - target))^(1 / cf["h"]))
         attr(v, "extrapolated") <- v > dmax
         v
       })
}


#' Sampling-adequacy analysis
#'
#' Sweeps region diameter (and optionally panel size), scoring the endpoint at
#' each vertex, then fits the saturation model and reports the adequacy
#' threshold: the diameter at which the recoverable signal reaches `target` of
#' its plateau.
#'
#' @param counts Spots x genes count matrix.
#' @param coords Spot coordinates.
#' @param labels Named vector of the spatially-defined label per spot (the
#'   archetype ground truth the region belongs to).
#' @param diameters Vector of region diameters to sweep.
#' @param panel_sizes Optional vector of gene-panel sizes (top-variance genes).
#' @param n_regions Regions per vertex.
#' @param scorer Function `(feature_matrix, labels)` returning a scalar to
#'   maximise. Default: leave-one-out nearest-centroid accuracy.
#' @param target Fraction of plateau defining adequacy.
#' @param seed RNG seed.
#' @return A list with the `surface` data.frame, the fitted models per panel
#'   size, and the adequacy thresholds.
#' @export
sampling_adequacy <- function(counts, coords, labels, diameters = c(0.5, 1, 2, 3, 4, 5, 6),
                              panel_sizes = NULL, n_regions = 120L,
                              scorer = NULL, target = 0.95, seed = 1L) {
  # Keep sparse: a Visium section is ~20k genes x ~3k spots and densifying it
  # costs ~500 MB before any work happens.
  sparse_in <- inherits(counts, "sparseMatrix")
  if (!sparse_in) counts <- as.matrix(counts)
  # PITFALL (hit 2026-09-17): as.character() is as.vector(x, "character"), which
  # STRIPS ATTRIBUTES INCLUDING NAMES. Coercing first therefore destroyed the
  # barcode->label mapping, the code fell through to rownames(coords) (integer
  # indices), nothing matched, every region was labelled NA, and the sweep died
  # with a misleading "check coordinates and labels" message. Preserve names
  # explicitly around the coercion.
  nm <- names(labels)
  labels <- as.character(labels)
  if (!is.null(nm)) {
    names(labels) <- nm
  } else if (!is.null(rownames(coords))) {
    names(labels) <- rownames(coords)[seq_along(labels)]
  } else {
    stop("labels have no names and coords has no row names, so spots cannot be ",
         "matched to labels. Pass a NAMED label vector keyed by barcode.",
         call. = FALSE)
  }
  if (is.null(scorer)) scorer <- scorer_nearest_centroid

  panel_sizes <- panel_sizes %||% min(200L, ncol(counts))
  rows <- list()
  skipped <- character(0)
  for (dn in diameters) {
    # A diameter smaller than the spot pitch cannot contain the minimum number of
    # spots. That is a property of the assay geometry, not an error: skip the
    # vertex and report it, rather than aborting the whole sweep. Dying here gave
    # a misleading "check coordinates and labels" message.
    pb <- tryCatch(pseudobulk_regions(counts, coords, dn, n_regions = n_regions,
                                      min_spots = max(3L, 1L), seed = seed),
                   error = function(e) { skipped <<- c(skipped, as.character(dn)); NULL })
    if (is.null(pb)) next
    ## the region label is the majority label of its spots
    lab <- vapply(pb$spot_index, function(s) {
      tb <- table(labels[match(s, names(labels))])
      if (!length(tb)) return(NA_character_)
      names(tb)[which.max(tb)]
    }, character(1))
    keep <- !is.na(lab)
    if (sum(keep) < 8L) next
    B <- pb$bulk[keep, , drop = FALSE]; L <- lab[keep]
    for (g in panel_sizes) {
      gg <- min(g, ncol(B))
      v <- matrixStats::colVars(B)
      sel <- order(v, decreasing = TRUE)[seq_len(gg)]
      rows[[length(rows) + 1L]] <- data.frame(
        diameter_mm = dn, panel_genes = gg,
        performance = suppressWarnings(scorer(B[, sel, drop = FALSE], L)),
        n_regions = nrow(B), stringsAsFactors = FALSE)
    }
  }
  if (!length(rows)) {
    stop("no usable sampling vertices. Diameters skipped as infeasible: ",
         if (length(skipped)) paste(skipped, collapse = ", ") else "none",
         ". The smallest requested diameter may be below the spot pitch; check ",
         "the spacing of `coords` and the units of `diameters`.", call. = FALSE)
  }
  if (length(skipped)) {
    message("NOTE: skipped ", length(skipped),
            " diameter(s) below the spot pitch: ", paste(skipped, collapse = ", "))
  }
  surface <- do.call(rbind, rows)

  fits <- lapply(split(surface, surface$panel_genes), function(df) {
    fit_hill(df$diameter_mm, df$performance)
  })
  thresh <- do.call(rbind, lapply(names(fits), function(g) {
    f <- fits[[g]]
    data.frame(panel_genes = as.integer(g),
               plateau = if (is.null(f)) NA_real_ else unname(f$coef["Pinf"]),
               half_saturation_mm = if (is.null(f)) NA_real_ else unname(f$coef["K"]),
               adequacy_mm = if (is.null(f)) NA_real_ else f$d_star(target),
               stringsAsFactors = FALSE)
  }))
  # flag any panel size whose adequacy diameter lies beyond the sampled range
  thresh$extrapolated <- vapply(seq_len(nrow(thresh)), function(i) {
    f <- fits[[as.character(thresh$panel_genes[i])]]
    if (is.null(f) || !is.finite(thresh$adequacy_mm[i])) return(NA)
    isTRUE(attr(f$d_star(target), "extrapolated"))
  }, logical(1))
  structure(list(surface = surface, fits = fits, adequacy = thresh,
                 target = target, diameters = diameters, panel_sizes = panel_sizes),
            class = "recoverR_adequacy")
}

#' Default scorer: leave-one-out nearest-centroid accuracy
#'
#' Deliberately simple and deterministic. A heavier learner would confound the
#' sampling question with model capacity, which is not what is being measured
#' here.
#' @keywords internal
scorer_nearest_centroid <- function(X, labels) {
  X <- as.matrix(X); labels <- factor(labels)
  if (nlevels(labels) < 2L || nrow(X) < 4L) return(NA_real_)
  ok <- 0L
  for (i in seq_len(nrow(X))) {
    tr <- setdiff(seq_len(nrow(X)), i)
    cen <- do.call(rbind, lapply(levels(labels), function(g) {
      ix <- tr[labels[tr] == g]
      if (!length(ix)) return(rep(NA_real_, ncol(X)))
      colMeans(X[ix, , drop = FALSE])
    }))
    rownames(cen) <- levels(labels)
    if (any(!is.finite(cen))) next
    d <- sqrt(rowSums((cen - matrix(X[i, ], nrow = nlevels(labels),
                                   ncol = ncol(X), byrow = TRUE))^2))
    if (levels(labels)[which.min(d)] == as.character(labels[i])) ok <- ok + 1L
  }
  ok / nrow(X)
}

#' @export
print.recoverR_adequacy <- function(x, ...) {
  cat("<recoverR_adequacy>  target =", x$target, "of plateau\n\n")
  cat("Sampling surface (region diameter x panel size):\n")
  s <- x$surface
  s$performance <- round(s$performance, 3)
  print(utils::head(s[order(s$panel_genes, s$diameter_mm), ], 24), row.names = FALSE)
  cat("\nAdequacy thresholds:\n")
  a <- x$adequacy; a$plateau <- round(a$plateau, 3)
  a$half_saturation_mm <- round(a$half_saturation_mm, 2)
  a$adequacy_mm <- round(a$adequacy_mm, 2)
  print(a, row.names = FALSE)
  cat("\nInterpretation: a specimen or trial sampling scheme must reach the\n")
  cat("adequacy diameter at the planned panel size, or the spatial endpoint\n")
  cat("cannot be recovered from that sample.\n")
  invisible(x)
}
