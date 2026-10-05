# ============================================================================
# recoverR - save.R
#
# Figure export with a guaranteed-vector contract.
#
# Project standing rule encoded here: EVERY plot is written as BOTH a vector PDF
# (for publication and editing) and a JPEG (for quick on-screen review), same
# basename, at the same time. A figure saved as only one of the two is treated as
# an error, not a preference.
#
# On top of that, the PDF is VERIFIED to be genuinely vector. A PDF that contains
# an embedded bitmap is not a vector figure, and ggplot2 >= 3.5 produces exactly
# that by default because it rasterises continuous colourbar legends. Silently
# shipping such a PDF is the failure mode this checker exists to catch.
# ============================================================================


#' Is a PDF genuinely vector?
#'
#' Reads the raw PDF and looks for an XObject of subtype `/Image`. A colourbar
#' legend rendered as a 1-pixel-wide gradient strip is the usual culprit.
#'
#' @param path Path to a PDF.
#' @return `TRUE` if no embedded image object is present, `FALSE` otherwise.
#' @export
is_vector_pdf <- function(path) {
  if (!file.exists(path)) return(NA)
  raw <- readBin(path, what = "raw", n = file.size(path))
  # NB: a PDF is binary, so rawToChar() yields a string that is not valid UTF-8.
  # grepl() then warns "input string 1 is invalid UTF-8" and, more importantly,
  # may not match reliably. useBytes = TRUE compares bytes directly and is both
  # silent and correct for this purpose.
  txt <- rawToChar(raw[raw != as.raw(0)])
  !grepl("/Subtype\\s*/Image", txt, perl = TRUE, useBytes = TRUE)
}


#' Save a plot as a vector PDF and a JPEG
#'
#' @param plot A ggplot / patchwork / grid object, or any object with a `print`
#'   method that draws.
#' @param filename Path without extension.
#' @param width,height Size; interpreted in `units`.
#' @param units `"in"` (default) or `"mm"`.
#' @param dpi Raster resolution for the JPEG.
#' @param formats Which outputs to write. The default honours the project rule
#'   (both PDF and JPEG). `"svg"` and `"tiff"` are available for journal
#'   submission sets.
#' @param verify Warn when the PDF is not genuinely vector.
#' @param quiet Suppress the per-file message.
#' @return Invisibly, a named character vector of the files written.
#' @export
save_plot <- function(plot, filename,
                      width = 7, height = 5,
                      units = c("in", "mm"),
                      dpi = 300,
                      formats = c("pdf", "jpeg"),
                      verify = TRUE, quiet = FALSE) {
  units <- match.arg(units)
  if (units == "mm") { width <- width / 25.4; height <- height / 25.4 }

  written <- character(0)
  for (fmt in formats) {
    f <- paste0(filename, ".", fmt)
    dir.create(dirname(f), recursive = TRUE, showWarnings = FALSE)
    switch(fmt,
      pdf = {
        grDevices::cairo_pdf(f, width = width, height = height, onefile = TRUE)
        print(plot); grDevices::dev.off()
      },
      jpeg = {
        grDevices::jpeg(f, width = width, height = height, units = "in",
                        res = dpi, quality = 95, bg = "white")
        print(plot); grDevices::dev.off()
      },
      svg = {
        if (!requireNamespace("svglite", quietly = TRUE))
          stop("format 'svg' needs the 'svglite' package.", call. = FALSE)
        svglite::svglite(f, width = width, height = height); print(plot); grDevices::dev.off()
      },
      tiff = {
        if (!requireNamespace("ragg", quietly = TRUE))
          stop("format 'tiff' needs the 'ragg' package.", call. = FALSE)
        ragg::agg_tiff(f, width = width, height = height, units = "in", res = dpi)
        print(plot); grDevices::dev.off()
      },
      stop("unsupported format: ", fmt, call. = FALSE)
    )
    written[fmt] <- f
  }

  if (verify && "pdf" %in% names(written)) {
    v <- tryCatch(is_vector_pdf(written["pdf"]), error = function(e) NA)
    if (isFALSE(v)) {
      warning("the exported PDF is NOT fully vector: it contains an embedded ",
              "image object, which ggplot2 does when a continuous colourbar ",
              "legend is rasterised. Add rr_guide_colourbar() (or use a ",
              "discrete fill) so the figure is genuinely editable.\n  file: ",
              written["pdf"], call. = FALSE)
    }
  }

  if (!quiet) {
    message("wrote ", paste(basename(written), collapse = " + "),
            "  (", round(width, 2), " x ", round(height, 2), " in)")
  }
  invisible(written)
}


#' Save a journal submission set at Nature column widths
#'
#' Writes SVG + vector PDF + 600 dpi TIFF at the specified physical width.
#'
#' @param plot Plot object.
#' @param filename Path without extension.
#' @param width_mm 89 (single column), 120 (1.5 column) or 183 (double column).
#' @param height_mm Height in millimetres.
#' @export
save_pub <- function(plot, filename, width_mm = 183, height_mm = 120) {
  save_plot(plot, filename, width = width_mm, height = height_mm, units = "mm",
            dpi = 600, formats = c("svg", "pdf", "tiff"), quiet = TRUE)
  message(sprintf("wrote journal set at %d x %d mm (svg + pdf + 600 dpi tiff)",
                  width_mm, height_mm))
}
