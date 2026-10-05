# ============================================================================
# recoverR - spatial-io.R
#
# Reading 10x Visium output and deriving spatially-defined archetype labels.
#
# WHY THIS IS IN THE PACKAGE
#   Everything downstream of a bulk->spatial claim depends on the ground truth
#   being defined the SAME way every time. If archetype labelling is left to each
#   analysis script, two papers defining "immune-excluded" differently will
#   produce incompatible results and neither will replicate. The labelling rule
#   lives here, explicitly and inspectably.
#
# THE LABELLING LOGIC
#   The operational definition in the literature is spatial, not compositional:
#     inflamed  - immune cells present AND infiltrating the tumour compartment
#     excluded  - immune cells present but held in stroma / at the margin
#     desert    - immune cells sparse throughout
#   A single bulk average cannot express the middle case, which is the entire
#   point of the package. The implementation therefore scores every spot for a
#   tumour programme and an immune programme, assigns each spot to a tumour or
#   stromal compartment, and labels the REGION from where the immune signal sits.
# ============================================================================


# ----------------------------------------------------------------------------
# Programme definitions
# ----------------------------------------------------------------------------

#' Marker programmes used for spatial compartment and archetype scoring
#'
#' Deliberately small, canonical gene sets. A large deconvolution reference would
#' confound the spatial question with reference quality; the object of study here
#' is where the signal sits, not how many cell types can be enumerated.
#' @export
rr_spatial_programmes <- function() {
  list(
    tumour = c("EPCAM", "KRT7", "KRT8", "KRT18", "KRT19", "NKX2-1", "NAPSA",
               "SFTPB", "SFTPC", "MUC1", "CEACAM5", "TACSTD2"),
    immune = c("CD3D", "CD3E", "CD2", "CD8A", "CD8B", "GZMA", "GZMB", "PRF1",
               "NKG7", "CXCL9", "CXCL10", "IFNG", "STAT1", "CCL5"),
    stroma = c("COL1A1", "COL1A2", "COL3A1", "FAP", "ACTA2", "DCN", "PDGFRB",
               "POSTN", "LUM", "THY1")
  )
}

#' Score a spot x gene matrix for each programme
#'
#' Scores are mean z-scores across spots within a sample, so they are relative to
#' that section rather than absolute. Absolute comparison across sections is not
#' meaningful for Visium (different capture efficiency and cellularity).
#' @keywords internal
score_programmes <- function(counts, programmes, log_norm = TRUE) {
  X <- as.matrix(counts)
  if (log_norm) {
    lib <- colSums(X)
    lib[lib == 0] <- 1
    X <- log1p(sweep(X, 2, lib / stats::median(lib[lib > 0]), "/"))
  }
  out <- matrix(NA_real_, nrow = ncol(X), ncol = length(programmes),
                dimnames = list(colnames(X), names(programmes)))
  present <- list()
  for (nm in names(programmes)) {
    g <- intersect(programmes[[nm]], rownames(X))
    present[[nm]] <- g
    if (length(g) < 2L) next
    sub <- X[g, , drop = FALSE]
    z <- t(scale(t(sub)))
    z[!is.finite(z)] <- 0
    out[, nm] <- colMeans(z, na.rm = TRUE)
  }
  attr(out, "genes_present") <- present
  out
}


# ----------------------------------------------------------------------------
# Visium reader
# ----------------------------------------------------------------------------

#' Read one 10x Visium sample
#'
#' Reads a `filtered_feature_bc_matrix.h5` and its matching `spatial/` folder.
#' Handles both `tissue_positions_list.csv` (Space Ranger < 2.0, no header) and
#' `tissue_positions.csv` (>= 2.0, with header) - the two are NOT
#' interchangeable and silently reading the wrong one loses the barcode column.
#'
#' @param sample_dir Directory containing the `.h5` and a `spatial/` subfolder,
#'   or the directory that *contains* a single such sample.
#' @param sample_id Optional label; defaults to the directory name.
#' @param h5 Optional explicit path to the `.h5`.
#' @return A list with `counts` (genes x spots, sparse), `features`, `barcodes`,
#'   `coords` (barcode, x, y, in_tissue, array_row, array_col) and `scalefactors`.
#' @export
read_visium <- function(sample_dir, sample_id = NULL, h5 = NULL) {
  if (!dir.exists(sample_dir)) stop("no such directory: ", sample_dir, call. = FALSE)

  if (is.null(h5)) {
    cand <- list.files(sample_dir, pattern = "filtered_feature_bc_matrix\\.h5$",
                       full.names = TRUE, recursive = TRUE)
    if (!length(cand)) {
      cand <- list.files(sample_dir, pattern = "\\.h5$", full.names = TRUE, recursive = TRUE)
    }
    if (!length(cand)) {
      stop("no .h5 count matrix found under ", sample_dir,
           ".\nIf you only have the *_spatial.tar.gz, extract it first; the counts ",
           "and the spatial folder are separate downloads.", call. = FALSE)
    }
    if (length(cand) > 1L) {
      stop("multiple .h5 files found under ", sample_dir, ": ",
           paste(basename(cand), collapse = ", "), ". Pass `h5=` explicitly.",
           call. = FALSE)
    }
    h5 <- cand[1]
  }

  if (!requireNamespace("rhdf5", quietly = TRUE)) {
    stop("read_visium() needs the 'rhdf5' package (Bioconductor).", call. = FALSE)
  }
  m <- rhdf5::h5read(h5, "/matrix")
  shape <- as.integer(m$shape)                 # genes x spots (CSC)
  counts <- Matrix::sparseMatrix(
    i = as.integer(m$indices) + 1L,            # 0-based -> 1-based
    p = as.integer(m$indptr),
    x = as.numeric(m$data),
    dims = shape, index1 = TRUE,
    dimnames = list(NULL, NULL)
  )
  genes <- if (!is.null(m$features) && !is.null(m$features$name)) m$features$name else m$genes
  rownames(counts) <- as.character(genes)
  colnames(counts) <- as.character(m$barcodes)

  ## ---- spatial coordinates ------------------------------------------------
  spdir <- file.path(dirname(h5), "spatial")
  if (!dir.exists(spdir)) {
    spdir <- list.files(sample_dir, pattern = "spatial$", full.names = TRUE,
                        recursive = TRUE)
    spdir <- if (length(spdir)) spdir[1] else spdir
  }
  coords <- NULL; scalefactors <- NULL
  if (length(spdir) && dir.exists(spdir)) {
    pos_new <- file.path(spdir, "tissue_positions.csv")
    pos_old <- file.path(spdir, "tissue_positions_list.csv")
    if (file.exists(pos_new)) {
      coords <- utils::read.csv(pos_new, stringsAsFactors = FALSE)
      names(coords)[1] <- "barcode"
    } else if (file.exists(pos_old)) {
      coords <- utils::read.csv(pos_old, header = FALSE, stringsAsFactors = FALSE)
      names(coords) <- c("barcode", "in_tissue", "array_row", "array_col",
                         "pxl_row_in_fullres", "pxl_col_in_fullres")
    } else {
      warning("a spatial/ folder exists at ", spdir,
              " but contains neither tissue_positions.csv nor ",
              "tissue_positions_list.csv.", call. = FALSE)
    }
    if (!is.null(coords)) {
      nm <- names(coords)
      coords$x <- coords[[grep("pxl_col_in_fullres|pxl_col", nm)[1]]]
      coords$y <- coords[[grep("pxl_row_in_fullres|pxl_row", nm)[1]]]
      coords <- coords[, c("barcode", "in_tissue", "array_row", "array_col", "x", "y")]
      ok <- coords$barcode %in% colnames(counts)
      if (sum(ok) == 0L) {
        stop("no barcode overlap between the .h5 and tissue_positions. The two ",
             "files are probably from different samples.", call. = FALSE)
      }
      coords <- coords[ok, , drop = FALSE]
      # make the row names meaningful: downstream helpers fall back to
      # rownames(coords) when a label vector is unnamed, so integer indices there
      # would silently match nothing.
      rownames(coords) <- coords$barcode
      counts <- counts[, coords$barcode, drop = FALSE]
    }
    sf <- file.path(spdir, "scalefactors_json.json")
    if (file.exists(sf)) {
      scalefactors <- tryCatch(jsonlite::fromJSON(sf), error = function(e) NULL)
    }
  }

  # Coordinates are REQUIRED by every downstream spatial function. Returning
  # silently without them would let a caller run a sampling-adequacy analysis on
  # a sample with no geometry and get a plausible-looking wrong answer, so their
  # absence always warns -- including when the spatial/ folder is missing entirely
  # (the common case when only the .h5 was downloaded and the separate
  # *_spatial.tar.gz was not).
  if (is.null(coords)) {
    warning("no spatial coordinates for sample '",
            sample_id %||% basename(normalizePath(sample_dir, mustWork = FALSE)),
            "': neither a spatial/ folder nor a tissue_positions file was found.\n",
            "  Coordinates are REQUIRED for sampling adequacy and for region-based ",
            "recoverability.\n  Check that the *_spatial.tar.gz was downloaded AND ",
            "extracted next to the .h5.", call. = FALSE)
  }

  structure(list(
    sample_id = sample_id %||% basename(normalizePath(sample_dir, mustWork = FALSE)),
    counts = counts, coords = coords, scalefactors = scalefactors, h5 = h5
  ), class = "recoverR_visium")
}

#' @export
print.recoverR_visium <- function(x, ...) {
  cat("<recoverR_visium>", x$sample_id, "\n")
  cat("  counts :", nrow(x$counts), "genes x", ncol(x$counts), "spots\n")
  if (is.null(x$coords)) {
    cat("  coords : MISSING\n")
  } else {
    cat("  coords :", nrow(x$coords), "spots;",
        sum(x$coords$in_tissue == 1), "in tissue;",
        sprintf("span %.0f x %.0f px", diff(range(x$coords$x)), diff(range(x$coords$y))), "\n")
  }
  invisible(x)
}


# ----------------------------------------------------------------------------
# Archetype labelling
# ----------------------------------------------------------------------------

#' Assign each spot to a tumour or stromal compartment
#'
#' A spot is tumour-compartment when its tumour programme score exceeds its
#' stromal score and is positive. The margin is reported so the sensitivity of
#' every downstream archetype call to this one threshold is auditable.
#' @param scores Programme scores (spots x programmes).
#' @param margin Extra margin the tumour score must exceed the stromal score by.
#' @export
spot_compartment <- function(scores, margin = 0) {
  tu <- scores[, "tumour"]; st <- scores[, "stroma"]
  out <- ifelse(!is.finite(tu) | !is.finite(st), NA_character_,
         ifelse(tu > st + margin & tu > 0, "tumour", "stroma"))
  out
}

#' Apply the three-way archetype rule to per-region fractions
#'
#' **The single implementation of the rule.** Two routes feed it: the Visium route
#' ([derive_spatial_archetypes()]), where the fractions are derived from spot-level programme
#' scores, and the imaging route ([read_roi_table()]), where they come from a segmentation's
#' per-ROI counts. The rule is the scientific object, so it exists once - if the two routes
#' each carried their own copy they would be free to drift, and the drift would be invisible
#' until someone compared two published archetype calls.
#'
#' The call is made from **where** the immune signal sits, not how much of it there is:
#' * `f_immune < desert_cut` -> `desert`
#' * else `f_tumour_immune >= inflamed_cut` -> `inflamed`
#' * else -> `excluded`
#'
#' @param f_immune Fraction of the region that is immune-high.
#' @param f_tumour_immune Fraction of the TUMOUR compartment that is immune-high.
#' @param inflamed_cut Tumour-compartment fraction at or above which the call is `inflamed`.
#' @param desert_cut Total immune fraction below which the call is `desert`.
#' @return A character vector of `"desert"` / `"inflamed"` / `"excluded"`, or `NA` where either
#'   input is not finite - a region with **no tumour compartment is UNDEFINED rather than
#'   forced into a class**, which is a property of the tissue and not a defect.
#' @keywords internal
call_archetype <- function(f_immune, f_tumour_immune, inflamed_cut = 0.25, desert_cut = 0.05) {
  ifelse(!is.finite(f_immune) | !is.finite(f_tumour_immune), NA_character_,
  ifelse(f_immune < desert_cut, "desert",
  ifelse(f_tumour_immune >= inflamed_cut, "inflamed", "excluded")))
}

#' Derive a region-level immune archetype from spot-level spatial data
#'
#' The three-way call is made from WHERE the immune signal sits, not how much of
#' it there is:
#'
#' \deqn{f_{\text{imm}} = \frac{n_{\text{immune-high}}}{n_{\text{total}}}, \qquad
#'       f_{\text{tum-imm}} = \frac{n_{\text{immune-high} \,\cap\, \text{tumour}}}{n_{\text{tumour}}}}
#'
#' * `f_imm < desert_cut` -> **desert**
#' * else `f_tum_imm >= inflamed_cut` -> **inflamed**
#' * else -> **excluded**
#'
#' @param scores Programme scores (spots x programmes), as produced by
#'   [rr_spatial_programmes()] scoring.
#' @param immune_cut Immune z-score above which a spot counts as immune-high.
#' @param inflamed_cut Fraction of tumour spots that must be immune-high.
#' @param desert_cut Fraction of all spots that must be immune-high to avoid desert.
#' @param compartment Precomputed compartments, or `NULL` to compute.
#' @return A `data.frame` of per-region fractions and the call, plus the
#'   intermediate quantities so the decision is never opaque.
#' @export
derive_spatial_archetypes <- function(scores, immune_cut = 0.5, inflamed_cut = 0.25,
                                      desert_cut = 0.05, compartment = NULL) {
  if (is.null(compartment)) compartment <- spot_compartment(scores)
  imm_high <- is.finite(scores[, "immune"]) & scores[, "immune"] > immune_cut
  is_tum <- compartment == "tumour"

  n <- length(imm_high)
  n_tum <- sum(is_tum, na.rm = TRUE)
  f_imm <- mean(imm_high, na.rm = TRUE)
  f_tum_imm <- if (n_tum > 0) sum(imm_high & is_tum, na.rm = TRUE) / n_tum else NA_real_
  f_str_imm <- if (sum(!is_tum, na.rm = TRUE) > 0)
    sum(imm_high & !is_tum, na.rm = TRUE) / sum(!is_tum, na.rm = TRUE) else NA_real_

  ## the rule itself lives in call_archetype(), shared with the imaging route
  call <- call_archetype(f_imm, f_tum_imm, inflamed_cut = inflamed_cut, desert_cut = desert_cut)

  data.frame(n_spots = n, n_tumour = n_tum,
             f_immune = f_imm, f_tumour_immune = f_tum_imm,
             f_stroma_immune = f_str_imm, archetype = call,
             stringsAsFactors = FALSE)
}

#' Shape an imaging study's per-ROI table into the archetype rule's inputs
#'
#' The package's imaging arm (multiplex immunofluorescence) measures **the same three
#' programmes the Visium route infers** - panCK for the tumour compartment, CD3/CD8 for immune
#' cells, alpha-SMA/FAP for stroma - but directly, at single-cell resolution, over centimetres
#' rather than a 55 um spot grid. A segmentation reports this as **counts per region of
#' interest**: how many cells sit in the tumour and stromal compartments, and how many of those
#' are marker-positive.
#'
#' This turns that table into exactly the columns the Visium route produces
#' (`n_spots`, `n_tumour`, `f_immune`, `f_tumour_immune`, `f_stroma_immune`, `archetype`,
#' `region`, `sample_id`) and applies **the same rule** ([call_archetype()]), so the existing
#' harness - [fit_recoverability()] -> [predict_recoverability()] -> [reliability_curve()] ->
#' [archetype_separability()] -> [sampling_adequacy()] - runs on imaging data **unchanged**.
#' That equivalence is the point: an IF cohort and a Visium cohort then answer the same
#' question in the same language.
#'
#' The counts must be **counts, not fractions or areas** (a fraction in a count column is
#' caught below rather than silently producing a fraction-of-a-fraction). Any further columns
#' are passed through, so a ROI's size can travel with it - which is what a scale question
#' (e.g. "is the circuit coherent at 3 mm?") needs.
#'
#' @param x A `data.frame`, or a path to a `.csv` / `.tsv` written by the imaging pipeline.
#' @param roi,sample,n_tumour,n_stroma,n_tumour_immune,n_stroma_immune Column names holding
#'   the ROI id, the sample id, and the four compartment counts.
#' @param inflamed_cut,desert_cut Thresholds, passed to [call_archetype()].
#' @return A `data.frame` with one row per ROI, in the Visium region-table's column names plus
#'   any extra columns of `x`.
#' @export
read_roi_table <- function(x, roi = "roi_id", sample = "sample_id",
                           n_tumour = "n_tumour", n_stroma = "n_stroma",
                           n_tumour_immune = "n_tumour_immune",
                           n_stroma_immune = "n_stroma_immune",
                           inflamed_cut = 0.25, desert_cut = 0.05) {
  if (is.character(x) && length(x) == 1L && !is.data.frame(x)) {
    if (!file.exists(x)) stop("file not found: ", x, call. = FALSE)
    sep <- if (grepl("\\.(tsv|txt)$", x, ignore.case = TRUE)) "\t" else ","
    x <- utils::read.table(x, sep = sep, header = TRUE, stringsAsFactors = FALSE,
                           check.names = FALSE)
  }
  if (!is.data.frame(x))
    stop("`x` must be a data.frame or a path to a .csv/.tsv of per-ROI counts.", call. = FALSE)

  need <- c(roi, sample, n_tumour, n_stroma, n_tumour_immune, n_stroma_immune)
  miss <- setdiff(need, names(x))
  if (length(miss))
    stop("missing column(s): ", paste(miss, collapse = ", "), ".\nFound: ",
         paste(names(x), collapse = ", "),
         "\nPass the real names through the `roi`/`sample`/`n_*` arguments - do not rename ",
         "silently.", call. = FALSE)

  as_count <- function(v, nm) {
    out <- suppressWarnings(as.numeric(v))
    bad <- !is.na(out) & (!is.finite(out) | out < 0)
    if (any(bad))
      stop("column '", nm, "' has ", sum(bad), " negative or non-finite value(s).", call. = FALSE)
    if (any(!is.na(v) & is.na(out)))
      stop("column '", nm, "' has non-numeric values; it must hold COUNTS, not labels.",
           call. = FALSE)
    out
  }
  nT  <- as_count(x[[n_tumour]], n_tumour)
  nS  <- as_count(x[[n_stroma]], n_stroma)
  nTI <- as_count(x[[n_tumour_immune]], n_tumour_immune)
  nSI <- as_count(x[[n_stroma_immune]], n_stroma_immune)

  ## impossible counts: the commonest way this gets fed the wrong thing is a table whose
  ## "immune" columns are FRACTIONS, which then exceed the denominator. Fail loudly.
  over <- which((!is.na(nTI) & !is.na(nT) & nTI > nT) | (!is.na(nSI) & !is.na(nS) & nSI > nS))
  if (length(over))
    stop(sum(over), " ROI(s) have more immune cells than compartment cells (first: ",
         x[[roi]][over[1]], "). Immune counts must be SUBSETS of their compartment; a table of ",
         "fractions or areas will trip this - pass counts.", call. = FALSE)

  n_tot <- nT + nS
  f_imm      <- ifelse(n_tot > 0, (nTI + nSI) / n_tot, NA_real_)
  f_tum_imm  <- ifelse(!is.na(nT) & nT > 0, nTI / nT, NA_real_)
  f_str_imm  <- ifelse(!is.na(nS) & nS > 0, nSI / nS, NA_real_)
  call       <- call_archetype(f_imm, f_tum_imm, inflamed_cut = inflamed_cut,
                               desert_cut = desert_cut)

  ## the Visium region-table's column names, so the harness needs no special case
  out <- data.frame(n_spots = n_tot, n_tumour = nT,
                    f_immune = f_imm, f_tumour_immune = f_tum_imm,
                    f_stroma_immune = f_str_imm, archetype = call,
                    region = as.character(x[[roi]]), sample_id = as.character(x[[sample]]),
                    stringsAsFactors = FALSE)
  extra <- setdiff(names(x), need)
  for (e in extra) out[[e]] <- x[[e]]
  attr(out, "n_undefined_no_tumour") <- sum(!is.na(n_tot) & n_tot > 0 & (is.na(nT) | nT == 0))
  attr(out, "n_no_cells")            <- sum(is.na(n_tot) | n_tot == 0)
  attr(out, "cutpoints")             <- c(inflamed_cut = inflamed_cut, desert_cut = desert_cut)
  out
}

#' Compare an original section with the sub-regions drawn from it
#'
#' A section-level call and a region-level call are different objects. Reporting
#' the section call alongside the region calls makes it explicit whether a region
#' was sampled from an inflamed, excluded or desert tumour - which is the
#' variable the bulk predictor must actually recover.
#'
#' @param scores Programme scores (spots x programmes).
#' @param ... Passed to [derive_spatial_archetypes()].
#' @return The section-level archetype call as a length-one character vector.
#' @export
section_archetype <- function(scores, ...) {
  d <- derive_spatial_archetypes(scores, ...)
  d$archetype[1]
}
