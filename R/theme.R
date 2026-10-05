# ============================================================================
# recoverR - theme.R
#
# Publication figure system.
#
# Design rules applied here:
#   * one restrained palette per figure: a neutral family for context/baselines,
#     a signal family for the quantity of interest, an accent family reserved for
#     directional cues (gain / harm). Deliberately NOT a rainbow: with 4-6
#     categories, maximal hue separation reads as decoration and hides which
#     series is the hero.
#   * ordered categories get an ordered (sequential) scale, never an unrelated
#     categorical one.
#   * no grid lines, thin axis rules, direct labels where a legend would force
#     eye travel.
#   * a colourbar legend is drawn as vector rectangles. ggplot2 >= 3.5 rasterises
#     continuous colourbars by default, which embeds a 1-px bitmap strip in the
#     PDF and silently breaks "this is a vector figure". That is a house rule
#     here, not a preference.
# ============================================================================


#' Font stack for figures
#'
#' Returns a font family available on the current system, preferring the
#' sans-serif stack used by most journals. Never hard-codes a font that may be
#' missing, because a missing `base_family` fails at draw time.
#' @export
rr_base_family <- function() {
  cands <- c("Arial", "Helvetica", "Nimbus Sans", "DejaVu Sans", "sans")
  if (requireNamespace("systemfonts", quietly = TRUE)) {
    have <- tryCatch(systemfonts::system_fonts()$family, error = function(e) character(0))
    hit <- cands[cands %in% have]
    if (length(hit)) return(hit[1])
  }
  ""   # device default
}

#' recoverR publication theme
#'
#' @param base_size Base font size in points. Journal single-column figures
#'   usually sit near 7; screen inspection is more comfortable near 11.
#' @param base_family Font family; defaults to [rr_base_family()].
#' @param grid Retain a faint grid? Off by default.
#' @export
theme_audit <- function(base_size = 9, base_family = rr_base_family(), grid = FALSE) {
  th <- ggplot2::theme_classic(base_size = base_size, base_family = base_family) +
    ggplot2::theme(
      axis.line        = ggplot2::element_line(linewidth = 0.35, colour = "grey20"),
      axis.ticks       = ggplot2::element_line(linewidth = 0.35, colour = "grey20"),
      axis.text        = ggplot2::element_text(colour = "grey15"),
      axis.title       = ggplot2::element_text(colour = "grey15"),
      plot.title       = ggplot2::element_text(face = "bold", size = base_size * 1.05,
                                               colour = "grey10", hjust = 0),
      plot.subtitle    = ggplot2::element_text(size = base_size * 0.88, colour = "grey35",
                                               margin = ggplot2::margin(b = 4)),
      plot.caption     = ggplot2::element_text(size = base_size * 0.78, colour = "grey45",
                                               hjust = 0),
      plot.tag         = ggplot2::element_text(face = "bold", size = base_size * 1.25,
                                               colour = "grey10"),
      plot.title.position = "plot",
      legend.title     = ggplot2::element_text(size = base_size * 0.9),
      legend.text      = ggplot2::element_text(size = base_size * 0.85),
      legend.key.size  = ggplot2::unit(0.85, "lines"),
      legend.background = ggplot2::element_blank(),
      legend.key       = ggplot2::element_blank(),
      strip.background = ggplot2::element_blank(),
      strip.text       = ggplot2::element_text(face = "bold", size = base_size * 0.92),
      panel.grid       = ggplot2::element_blank(),
      plot.background  = ggplot2::element_rect(fill = "white", colour = NA),
      panel.background = ggplot2::element_rect(fill = "white", colour = NA)
    )
  if (grid) {
    th <- th + ggplot2::theme(
      panel.grid.major = ggplot2::element_line(linewidth = 0.15, colour = "grey92"),
      panel.grid.minor = ggplot2::element_blank())
  }
  th
}


# ----------------------------------------------------------------------------
# Palettes
# ----------------------------------------------------------------------------

#' recoverR discrete palette
#'
#' A restrained qualitative palette: greys carry context and baselines, blues
#' carry the quantity of interest, and the warm accents are reserved for
#' directional meaning (harm vs gain). Ordered factors should use [rr_scale_seq()]
#' instead.
#' @export
rr_pal <- function() {
  c(neutral1 = "#9AA0A6", neutral2 = "#C9CDD1", neutral3 = "#5F6368",
    signal1  = "#1F4E79", signal2  = "#3C78A8", signal3  = "#7FB2D4",
    signal4  = "#BBD7EA",
    accent_warn = "#D98324", accent_harm = "#B23A48", accent_good = "#3F7D5B")
}

#' Sequential scale for ordered categories
#' @export
#' @param direction 1 for dark-to-light, -1 to reverse.
#' @param ... Passed through to the learner.
rr_scale_seq <- function(direction = 1, ...) {
  cols <- c("#0B2E4F", "#1F4E79", "#3C78A8", "#7FB2D4", "#C8E0EF")
  if (direction < 0) cols <- rev(cols)
  ggplot2::scale_fill_gradientn(colours = cols, ...)
}

#' Diverging scale for signed quantities
#' @export
#' @param midpoint Value mapped to the neutral midpoint of the diverging scale.
#' @param ... Passed through to the learner.
rr_scale_div <- function(midpoint = 0, ...) {
  ggplot2::scale_fill_gradient2(
    low = "#1F4E79", mid = "#F5F5F5", high = "#B23A48", midpoint = midpoint, ...)
}

#' Colourbar guide that stays vector
#'
#' ggplot2 >= 3.5 renders continuous colourbars as a raster gradient, embedding a
#' 1-px bitmap in the exported PDF. `raster = FALSE` is deprecated in those
#' versions and is silently ignored, so the supported spelling is used here.
#' @export
#' @param aesthetic Name of the aesthetic whose guide must stay vector, e.g. `"fill"`.
rr_guide_colourbar <- function(aesthetic = "fill") {
  g <- ggplot2::guide_colourbar(display = "rectangles")
  do.call(ggplot2::guides, stats::setNames(list(g), aesthetic))
}

#' Apply the recoverR theme and the vector colourbar guard to any plot
#' @export
#' @param p A ggplot (or patchwork) object.
#' @param base_size Base font size in points.
rr_finish <- function(p, base_size = 9) {
  p + theme_audit(base_size) + rr_guide_colourbar("fill") +
    ggplot2::theme(legend.position = "right")
}


# ----------------------------------------------------------------------------
# Shared plot helpers
# ----------------------------------------------------------------------------

#' A caption carrying n and the uncertainty definition
#'
#' Statistics and error-bar definitions belong to the figure, not to a caption
#' added later, so the helper makes it hard to omit them.
#' @export
#' @param n Number of observations.
#' @param what Label for the counted unit, e.g. `"samples"`.
#' @param extra Optional extra clause appended to the caption.
rr_caption_n <- function(n, what = "samples", extra = NULL) {
  paste0("n = ", n, " ", what,
         if (!is.null(extra)) paste0("; ", extra) else "")
}
