# ============================================================================
# recoverR - utils.R
# ============================================================================

#' Null-coalescing operator
#'
#' Provided locally so the package does not depend on the R version that
#' introduced it in base.
#' @name grapes-or-or-grapes
#' @keywords internal
`%||%` <- function(a, b) if (is.null(a)) b else a

#' Require an optional package with an actionable message
#' @keywords internal
need_pkg <- function(pkg, why) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop("This step needs the '", pkg, "' package (", why, ").\n",
         "Install with: install.packages('", pkg, "') or BiocManager::install('",
         pkg, "').", call. = FALSE)
  }
  invisible(TRUE)
}
