# ============================================================================
# recoverR - plot-attribution.R
#
# Figures for the layer-attribution results. Every function returns a ggplot so
# panels compose with patchwork; nothing writes to disk here (see save_plot()).
#
# Colour semantics used throughout: the SIGNAL family marks a layer that
# demonstrably contributes; the NEUTRAL family marks one that is not
# distinguishable from the permutation null. The accent family is reserved for
# harm (a layer that makes out-of-sample performance worse).
# ============================================================================


#' Shapley attribution of out-of-sample performance
#'
#' The hero panel: each view's exact Shapley value, with layers that cannot be
#' distinguished from the dimension-matched permutation null shown in neutral
#' grey. Values carry direct labels, so no legend round-trip is needed to read
#' magnitudes.
#'
#' @param x A `recoverR_attribution` object.
#' @param title,subtitle Plot titles; sensible defaults are supplied.
#' @param base_size Base font size.
#' @return A ggplot object.
#' @export
plot_shapley <- function(x, title = NULL, subtitle = NULL, base_size = 9) {
  d <- x$shapley
  d$verdict <- "not distinguishable from null"
  if (!is.null(x$null)) {
    m <- x$null[, c("view", "p_value", "delta_observed")]
    d <- merge(d, m, by = "view", all.x = TRUE)
    d$verdict <- ifelse(!is.na(d$p_value) & d$p_value < 0.05,
                        "contributes", "not distinguishable from null")
  } else {
    d$p_value <- NA_real_
  }
  d$verdict <- factor(d$verdict, levels = c("contributes", "not distinguishable from null"))
  pal <- c("contributes" = unname(rr_pal()["signal1"]),
           "not distinguishable from null" = unname(rr_pal()["neutral2"]
           ))
  p <- ggplot2::ggplot(d, ggplot2::aes(shapley, reorder(view, shapley), fill = verdict)) +
    ggplot2::geom_vline(xintercept = 0, linewidth = 0.35, colour = "grey55") +
    ggplot2::geom_col(width = 0.68) +
    ggplot2::geom_text(ggplot2::aes(label = sprintf("%+.4f", shapley)),
                       hjust = ifelse(d$shapley >= 0, -0.12, 1.12),
                       size = base_size * 0.30, colour = "grey20") +
    ggplot2::scale_fill_manual(values = pal, drop = FALSE) +
    ggplot2::labs(
      title = title %||% "Exact Shapley attribution of out-of-sample performance",
      subtitle = subtitle %||% sprintf(
        "sum of Shapley values = V(all views) = %.4f; grey = p >= 0.05 vs dimension-matched permutation null",
        unique(d$V_full)[1]),
      x = "Shapley value (contribution to out-of-sample performance)",
      y = NULL, fill = NULL) +
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0.18, 0.18))) +
    theme_audit(base_size) +
    ggplot2::theme(legend.position = "bottom")
  p
}

#' Leave-one-view-out loss
#'
#' Diverging bars: positive means removing the view hurt (it was helping),
#' negative means removing it HELPED (the view was actively harmful). The second
#' case is the one that matters most and is the one a variance decomposition
#' cannot reveal.
#'
#' @param x A `recoverR_attribution` object.
#' @export
#' @param base_size Base font size in points.
plot_leave_one_out <- function(x, base_size = 9) {
  d <- x$ablation$leave_one_out
  d$effect <- ifelse(d$loss_when_removed >= 0, "helps", "actively harmful")
  pal <- c("helps" = unname(rr_pal()["signal1"]),
           "actively harmful" = unname(rr_pal()["accent_harm"]))
  ggplot2::ggplot(d, ggplot2::aes(loss_when_removed, reorder(view, loss_when_removed),
                                  fill = effect)) +
    ggplot2::geom_vline(xintercept = 0, linewidth = 0.35, colour = "grey55") +
    ggplot2::geom_col(width = 0.68) +
    ggplot2::geom_text(ggplot2::aes(label = sprintf("%+.4f", loss_when_removed)),
                       hjust = ifelse(d$loss_when_removed >= 0, -0.12, 1.12),
                       size = base_size * 0.30, colour = "grey20") +
    ggplot2::scale_fill_manual(values = pal, drop = FALSE) +
    ggplot2::labs(title = "Leave-one-view-out loss in out-of-sample performance",
                  subtitle = "negative = removing the view improved performance",
                  x = "V(full) - V(full without view)", y = NULL, fill = NULL) +
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0.2, 0.2))) +
    theme_audit(base_size) + ggplot2::theme(legend.position = "bottom")
}

#' Forward-selection staircase
#'
#' Shows the marginal gain of each added view given the views already chosen. A
#' staircase that flattens immediately is the honest picture of a redundant view
#' set.
#'
#' @param x A `recoverR_attribution` object.
#' @export
#' @param base_size Base font size in points.
plot_forward_selection <- function(x, base_size = 9) {
  f <- x$ablation$forward
  if (is.null(f)) stop("no forward-selection table available", call. = FALSE)
  f$step <- factor(f$step)
  f$label <- sprintf("%s (+%.3f)", f$view, f$gain)
  ggplot2::ggplot(f, ggplot2::aes(step, V, group = 1)) +
    ggplot2::geom_line(linewidth = 0.6, colour = "grey55") +
    ggplot2::geom_point(size = 1.8, colour = unname(rr_pal()["signal1"])) +
    ggplot2::geom_text(ggplot2::aes(label = label), vjust = -0.9,
                       size = base_size * 0.30, colour = "grey20") +
    ggplot2::labs(title = "Forward view selection",
                  subtitle = "marginal gain at each step, given what was already added",
                  x = "step", y = "out-of-sample performance") +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.08, 0.18))) +
    theme_audit(base_size)
}

#' Compose the attribution figure set into one plate
#'
#' @param x A `recoverR_attribution` object.
#' @param base_size Base font size.
#' @export
plot_attribution_plate <- function(x, base_size = 9) {
  if (!requireNamespace("patchwork", quietly = TRUE)) {
    stop("plot_attribution_plate() needs the 'patchwork' package.", call. = FALSE)
  }
  a <- plot_shapley(x, base_size = base_size)
  b <- plot_leave_one_out(x, base_size = base_size)
  c <- tryCatch(plot_forward_selection(x, base_size = base_size), error = function(e) NULL)
  if (is.null(c)) return(a / b)
  (a | b) / c + patchwork::plot_layout(heights = c(1.15, 0.85))
}


#' Ablation waterfall across an arbitrary set of view subsets
#'
#' For results produced outside [attrib_views()], e.g. comparing integration
#' strategies rather than feature subsets. Ordered so the reading is immediate.
#'
#' @param df Data frame with columns `label`, `value`, and optionally `n_views`.
#' @param value_lab Axis label.
#' @param highlight Optional label to emphasise as the hero bar.
#' @export
#' @param base_size Base font size in points.
plot_ablation_bars <- function(df, value_lab = "out-of-sample performance",
                               highlight = NULL, base_size = 9) {
  stopifnot(all(c("label", "value") %in% names(df)))
  df$hero <- if (!is.null(highlight)) df$label == highlight else FALSE
  df$fill <- ifelse(df$hero, "hero", "context")
  pal <- c(hero = unname(rr_pal()["signal1"]), context = unname(rr_pal()["neutral2"]))
  ggplot2::ggplot(df, ggplot2::aes(value, reorder(label, value), fill = fill)) +
    ggplot2::geom_col(width = 0.7) +
    ggplot2::geom_text(ggplot2::aes(label = sprintf("%.3f", value)),
                       hjust = -0.12, size = base_size * 0.30, colour = "grey20") +
    ggplot2::scale_fill_manual(values = pal, guide = "none") +
    ggplot2::labs(title = "Ablation", x = value_lab, y = NULL) +
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0, 0.16))) +
    theme_audit(base_size)
}


#' Variance partitioning across omics layers under cross-validation
#'
#' Ordered categories use an ordered (sequential) scale, because the layers are
#' not interchangeable categories -- they differ in how much they are expected to
#' carry.
#'
#' @param df Data frame with columns `layer`, `cv_r2`.
#' @param baseline Optional reference value drawn as a dashed rule.
#' @export
#' @param base_size Base font size in points.
plot_variance_partition <- function(df, baseline = NULL, base_size = 9) {
  stopifnot(all(c("layer", "cv_r2") %in% names(df)))
  p <- ggplot2::ggplot(df, ggplot2::aes(cv_r2, reorder(layer, cv_r2), fill = cv_r2)) +
    ggplot2::geom_col(width = 0.7) +
    ggplot2::geom_text(ggplot2::aes(label = sprintf("%.3f", cv_r2)),
                       hjust = -0.12, size = base_size * 0.30, colour = "grey20") +
    rr_scale_seq() +
    ggplot2::labs(title = "Variance partitioning under cross-validation",
                  x = "cross-validated R-squared", y = NULL, fill = "CV R2") +
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0, 0.16))) +
    theme_audit(base_size)
  if (!is.null(baseline)) {
    p <- p + ggplot2::geom_vline(xintercept = baseline, linetype = "dashed",
                                 colour = unname(rr_pal()["accent_harm"]), linewidth = 0.4)
  }
  p + rr_guide_colourbar("fill")
}
