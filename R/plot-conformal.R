# ============================================================================
# recoverR - plot-conformal.R
#
# Figures for calibrated recoverability. The scientific message these figures
# must carry is not "how accurate is the classifier" but "where does the model
# know that it does not know". Abstention is therefore drawn as a first-class
# outcome, not hidden inside an accuracy number.
# ============================================================================


#' Reliability curve with expected calibration error
#'
#' Confidence is only meaningful if it tracks accuracy. The diagonal is the
#' identity line; deviation above it is over-confidence, which is the failure
#' that makes a bulk biomarker look better than it is.
#'
#' @param fit A `recoverR_conformal` object.
#' @param base_size Base font size.
#' @export
plot_reliability <- function(fit, base_size = 9) {
  rc <- reliability_curve(fit)
  d <- rc$curve
  ggplot2::ggplot(d, ggplot2::aes(mean_conf, observed_acc)) +
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed",
                         colour = "grey60", linewidth = 0.35) +
    ggplot2::geom_linerange(ggplot2::aes(ymin = pmin(observed_acc, mean_conf),
                                         ymax = pmax(observed_acc, mean_conf)),
                            colour = unname(rr_pal()["neutral2"]), linewidth = 0.5) +
    ggplot2::geom_point(ggplot2::aes(size = n), colour = unname(rr_pal()["signal1"])) +
    ggplot2::geom_line(colour = unname(rr_pal()["signal1"]), linewidth = 0.5, alpha = 0.7) +
    ggplot2::scale_size_continuous(range = c(1, 3.2), guide = "none") +
    ggplot2::labs(
      title = "Calibration of the recoverability model",
      subtitle = sprintf("expected calibration error = %.3f; dashed = perfect calibration",
                         rc$ece),
      x = "mean predicted confidence", y = "observed accuracy",
      caption = "point area is the number of calibration samples in the bin") +
    ggplot2::coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
    theme_audit(base_size)
}


#' Coverage / accuracy / abstention trade-off
#'
#' A single coverage figure can be met while the model is uninformative for most
#' samples, so the trade-off is shown directly: as the model is allowed to
#' abstain (moving right to left), what accuracy does it buy, and how many
#' samples does it refuse?
#'
#' @param fit A `recoverR_conformal` object.
#' @param alpha_target Draw a marker at the fitted coverage target.
#' @export
#' @param base_size Base font size in points.
plot_coverage_risk <- function(fit, alpha_target = TRUE, base_size = 9) {
  cr <- coverage_risk_curve(fit)
  cr <- cr[is.finite(cr$accuracy), ]
  long <- rbind(
    data.frame(answered = cr$answered, value = cr$accuracy, series = "accuracy when answered"),
    data.frame(answered = cr$answered, value = cr$singleton_rate, series = "singleton (committed) rate")
  )
  pal <- c("accuracy when answered" = unname(rr_pal()["signal1"]),
           "singleton (committed) rate" = unname(rr_pal()["neutral3"]))
  p <- ggplot2::ggplot(long, ggplot2::aes(answered, value, colour = series)) +
    ggplot2::geom_line(linewidth = 0.7) +
    ggplot2::scale_colour_manual(values = pal) +
    ggplot2::labs(title = "Selective prediction: what abstention buys",
                  subtitle = "left = model allowed to abstain; right = forced to answer everything",
                  x = "fraction of samples answered", y = "proportion", colour = NULL) +
    ggplot2::coord_cartesian(ylim = c(0, 1)) +
    theme_audit(base_size) + ggplot2::theme(legend.position = "bottom")
  if (alpha_target) {
    p <- p + ggplot2::geom_vline(xintercept = 1 - fit$alpha, linetype = "dotted",
                                 colour = unname(rr_pal()["accent_warn"]), linewidth = 0.4)
  }
  p
}


#' Pairwise archetype separability
#'
#' The decisive diagnostic. For each ordered pair of archetypes, the fraction of
#' calibration samples whose conformal prediction set contained BOTH. High values
#' identify pairs that bulk data cannot separate, which is a far more actionable
#' statement than an aggregate accuracy: it says which clinical distinction the
#' assay will not deliver.
#'
#' @param fit A `recoverR_conformal` object.
#' @param base_size Base font size.
#' @export
plot_separability <- function(fit, base_size = 9) {
  m <- archetype_separability(fit)
  d <- expand.grid(from = rownames(m), to = colnames(m), stringsAsFactors = FALSE)
  d$confusability <- as.vector(m)
  d$from <- factor(d$from, levels = rownames(m))
  d$to <- factor(d$to, levels = colnames(m))
  ggplot2::ggplot(d, ggplot2::aes(to, from, fill = confusability)) +
    ggplot2::geom_tile(colour = "white", linewidth = 0.7) +
    ggplot2::geom_text(ggplot2::aes(label = sprintf("%.2f", confusability)),
                       size = base_size * 0.32,
                       colour = ifelse(d$confusability > 0.5, "white", "grey15")) +
    rr_scale_seq() +
    rr_guide_colourbar("fill") +
    ggplot2::labs(
      title = "Which archetypes does a bulk profile fail to separate?",
      subtitle = "fraction of prediction sets containing BOTH labels; high = bulk cannot distinguish the pair",
      x = NULL, y = NULL, fill = "confusability") +
    ggplot2::coord_equal() +
    theme_audit(base_size) +
    ggplot2::theme(axis.line = ggplot2::element_blank(),
                   axis.ticks = ggplot2::element_blank())
}


#' Per-sample recoverability and abstention
#'
#' Each sample is placed by classifier confidence against the size of its
#' conformal prediction set. Committed (singleton) calls sit on the bottom row;
#' abstentions stack above. Points are jittered because confidence values are
#' discrete at small n.
#'
#' @param pred The `data.frame` returned by [predict_recoverability()].
#' @param base_size Base font size.
#' @export
plot_recoverability <- function(pred, base_size = 9) {
  stopifnot(all(c("confidence", "set_size") %in% names(pred)))
  pred$set_size_f <- factor(pred$set_size)
  pal <- c(unname(rr_pal()["neutral2"]), unname(rr_pal()["signal1"]),
           unname(rr_pal()["accent_warn"]), unname(rr_pal()["accent_harm"]),
           unname(rr_pal()["neutral3"]))
  n_abst <- sum(pred$abstain, na.rm = TRUE)
  ggplot2::ggplot(pred, ggplot2::aes(set_size_f, confidence)) +
    ggplot2::geom_jitter(ggplot2::aes(colour = abstain), width = 0.16, height = 0.006,
                         size = 1.1, alpha = 0.75) +
    ggplot2::scale_colour_manual(
      values = c(`FALSE` = unname(rr_pal()["signal1"]),
                 `TRUE` = unname(rr_pal()["accent_harm"])),
      labels = c(`FALSE` = "committed (singleton)", `TRUE` = "abstained"),
      name = NULL) +
    ggplot2::labs(
      title = "Calibrated recoverability per bulk sample",
      subtitle = sprintf("%d of %d samples abstained (%.0f%%) at alpha = %.2f",
                         n_abst, nrow(pred), 100 * n_abst / nrow(pred), pred$alpha[1]),
      x = "conformal prediction-set size", y = "classifier confidence",
      caption = "abstained samples receive no recoverability score by design") +
    ggplot2::coord_cartesian(ylim = c(0, 1)) +
    theme_audit(base_size) + ggplot2::theme(legend.position = "bottom")
}


#' Sampling-adequacy curve
#'
#' Signal retained against region diameter and panel size. Reproduces, as a
#' continuous design rule, the empirical observation that spatially predictive
#' signal collapses below roughly 3 mm of tissue and a few dozen genes.
#'
#' @param df Data frame with `diameter_mm`, `performance`, and optionally
#'   `panel_genes` (for one curve per panel size).
#' @param adequacy Target fraction of the plateau performance.
#' @export
#' @param base_size Base font size in points.
plot_sampling_curve <- function(df, adequacy = 0.95, base_size = 9) {
  stopifnot(all(c("diameter_mm", "performance") %in% names(df)))
  df$panel_f <- if ("panel_genes" %in% names(df)) factor(df$panel_genes) else factor("all genes")
  p <- ggplot2::ggplot(df, ggplot2::aes(diameter_mm, performance, colour = panel_f)) +
    ggplot2::geom_line(linewidth = 0.7) +
    ggplot2::geom_point(size = 1.3) +
    ggplot2::scale_colour_manual(values = unname(rr_pal()[c("signal4", "signal3", "signal2", "signal1",
                                                            "neutral3")])[seq_len(nlevels(df$panel_f))],
                                 name = "panel (genes)") +
    ggplot2::labs(
      title = "Sampling adequacy for a spatial endpoint",
      subtitle = sprintf("dotted rule marks %.0f%% of the plateau performance", 100 * adequacy),
      x = "region diameter (mm)", y = "recoverable signal") +
    theme_audit(base_size) + ggplot2::theme(legend.position = "bottom")
  plateau <- max(df$performance, na.rm = TRUE)
  p + ggplot2::geom_hline(yintercept = adequacy * plateau, linetype = "dotted",
                          colour = unname(rr_pal()["accent_warn"]), linewidth = 0.4)
}
