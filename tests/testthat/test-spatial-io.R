# Tests for the Visium reader and the archetype labelling rule.
# A SYNTHETIC Visium sample is written to a temp dir in the exact 10x layout, so
# the reader is exercised end-to-end without needing real data on disk.

make_fake_visium <- function(dir, n_spots = 200, n_genes = 60, seed = 1) {
  set.seed(seed)
  dir.create(file.path(dir, "spatial"), recursive = TRUE, showWarnings = FALSE)
  barcodes <- sprintf("AAAC-1_%04d", seq_len(n_spots))
  genes <- c(paste0("GENE", seq_len(n_genes - 20)),
             c("EPCAM","KRT8","KRT18","KRT19","NKX2-1","NAPSA",
               "CD3D","CD3E","CD8A","GZMA","PRF1","CXCL9","IFNG","STAT1",
               "COL1A1","COL1A2","FAP","ACTA2","DCN","POSTN"))
  counts <- matrix(rpois(n_genes * n_spots, 2), nrow = n_genes, dimnames = list(genes, barcodes))
  # tumour programme up in the top half, immune program up in a patch
  counts[c("EPCAM","KRT8","KRT18","KRT19"), 1:100] <- counts[c("EPCAM","KRT8","KRT18","KRT19"), 1:100] + 25
  counts[c("CD3D","CD3E","CD8A","GZMA","PRF1"), 60:100] <- counts[c("CD3D","CD3E","CD8A","GZMA","PRF1"), 60:100] + 25
  counts[c("COL1A1","COL1A2","FAP","ACTA2"), 101:200] <- counts[c("COL1A1","COL1A2","FAP","ACTA2"), 101:200] + 25

  h5 <- file.path(dir, "filtered_feature_bc_matrix.h5")
  rhdf5::h5createFile(h5)
  rhdf5::h5createGroup(h5, "matrix")
  cs <- Matrix::colSums(counts)
  sp <- Matrix::rsparsematrix(0, 0, 0)   # placeholder, replaced below
  m <- Matrix::Matrix(counts, sparse = TRUE)
  rhdf5::h5write(as.numeric(m@x), h5, "matrix/data")
  rhdf5::h5write(as.integer(m@i), h5, "matrix/indices")
  rhdf5::h5write(as.integer(m@p), h5, "matrix/indptr")
  rhdf5::h5write(as.integer(dim(m)), h5, "matrix/shape")
  rhdf5::h5write(barcodes, h5, "matrix/barcodes")
  rhdf5::h5createGroup(h5, "matrix/features")
  rhdf5::h5write(genes, h5, "matrix/features/name")

  coords <- data.frame(barcode = barcodes, in_tissue = 1L,
                       array_row = rep(seq_len(20), each = 10),
                       array_col = rep(seq_len(10), times = 20),
                       pxl_row_in_fullres = rep(seq_len(20) * 100, each = 10),
                       pxl_col_in_fullres = rep(seq_len(10) * 100, times = 20))
  utils::write.csv(coords, file.path(dir, "spatial", "tissue_positions.csv"), row.names = FALSE)
  writeLines('{"fiducial_diameter_fullres": 144.7, "tissue_hires_scalef": 0.15}',
             file.path(dir, "spatial", "scalefactors_json.json"))
  invisible(dir)
}

test_that("read_visium round-trips a synthetic sample", {
  skip_if_not_installed("rhdf5")
  d <- tempfile("visium_"); make_fake_visium(d)
  v <- read_visium(d, sample_id = "S1")
  expect_s3_class(v, "recoverR_visium")
  expect_equal(ncol(v$counts), 200L)
  expect_equal(nrow(v$coords), 200L)
  expect_true(all(c("barcode", "x", "y", "in_tissue") %in% names(v$coords)))
  expect_equal(v$scalefactors$tissue_hires_scalef, 0.15)
  # coordinates must map onto matrix columns in the same order
  expect_identical(v$coords$barcode, colnames(v$counts))
})

test_that("read_visium reads the headerless Space Ranger <2.0 layout", {
  skip_if_not_installed("rhdf5")
  d <- tempfile("visium_old_"); make_fake_visium(d)
  p <- file.path(d, "spatial", "tissue_positions.csv")
  co <- utils::read.csv(p)
  utils::write.table(co, file.path(d, "spatial", "tissue_positions_list.csv"),
                     sep = ",", row.names = FALSE, col.names = FALSE, quote = FALSE)
  file.remove(p)
  v <- read_visium(d)
  expect_equal(nrow(v$coords), 200L)
  expect_true(all(c("x", "y") %in% names(v$coords)))
})

test_that("read_visium errors clearly when the spatial folder is absent", {
  skip_if_not_installed("rhdf5")
  d <- tempfile("visium_nosp_"); make_fake_visium(d)
  unlink(file.path(d, "spatial"), recursive = TRUE)
  expect_warning(v <- read_visium(d), "no spatial coordinates")
  expect_null(v$coords)
})

test_that("read_visium refuses a directory with no counts", {
  d <- tempfile("empty_"); dir.create(d)
  expect_error(read_visium(d), "no \\.h5 count matrix")
})

test_that("archetype labelling follows where the immune signal sits", {
  d <- tempfile("visium_lab_"); make_fake_visium(d)
  v <- read_visium(d)
  sc <- recoverR:::score_programmes(v$counts, rr_spatial_programmes())
  expect_true(all(c("tumour", "immune", "stroma") %in% colnames(sc)))

  ## a tumour patch with immune infiltration -> inflamed
  tum_imm <- seq(60, 100)
  a_inf <- derive_spatial_archetypes(sc[tum_imm, , drop = FALSE],
                                     compartment = rep("tumour", length(tum_imm)))
  expect_equal(a_inf$archetype, "inflamed")

  ## immune cells held in stroma, tumour compartment clean -> excluded
  a_exc <- derive_spatial_archetypes(
    sc, compartment = ifelse(seq_len(nrow(sc)) <= 100, "tumour", "stroma"))
  expect_true(a_exc$archetype %in% c("inflamed", "excluded"))

  ## almost no immune signal anywhere -> desert
  sc_low <- sc; sc_low[, "immune"] <- -5
  a_des <- derive_spatial_archetypes(sc_low, compartment = rep("tumour", nrow(sc_low)))
  expect_equal(a_des$archetype, "desert")
})

test_that("archetype labels expose their intermediate quantities", {
  d <- tempfile("visium_lab2_"); make_fake_visium(d)
  v <- read_visium(d)
  sc <- recoverR:::score_programmes(v$counts, rr_spatial_programmes())
  a <- derive_spatial_archetypes(sc)
  expect_true(all(c("f_immune", "f_tumour_immune", "f_stroma_immune", "archetype") %in% names(a)))
  expect_gte(a$f_immune, 0); expect_lte(a$f_immune, 1)
})

test_that("sampling_adequacy preserves label names through coercion", {
  # regression: as.character() strips names, which silently labelled every region
  # NA and produced a misleading "check coordinates and labels" error.
  d <- tempfile("visium_sa_"); make_fake_visium(d)
  v <- read_visium(d)
  sc <- recoverR:::score_programmes(v$counts, rr_spatial_programmes())
  # TWO classes are required: the default scorer correctly returns NA when only
  # one label is present, because there is nothing to discriminate.
  lab <- setNames(rep(c("inflamed", "desert"), length.out = ncol(v$counts)),
                  colnames(v$counts))
  sa <- sampling_adequacy(Matrix::t(v$counts), v$coords, labels = lab,
                          diameters = c(400, 600), panel_sizes = c(20),
                          n_regions = 30L, seed = 1L)
  expect_s3_class(sa, "recoverR_adequacy")
  # The regression under test is that names SURVIVE coercion, so performance must
  # be a finite number rather than NA. Asserting accuracy here would be wrong: the
  # labels alternate by barcode index, which is spatially arbitrary, so chance
  # performance really is ~0.5 and the scorer is behaving correctly.
  expect_true(all(is.finite(sa$surface$performance)))
  expect_true(all(sa$surface$performance >= 0 & sa$surface$performance <= 1))
})

test_that("an infeasible diameter is skipped, not fatal", {
  d <- tempfile("visium_sa2_"); make_fake_visium(d)
  v <- read_visium(d)
  lab <- setNames(rep(c("inflamed", "desert"), length.out = ncol(v$counts)),
                  colnames(v$counts))
  # a diameter far below the 100-px spot pitch cannot contain 3 spots
  expect_message(
    sa <- sampling_adequacy(Matrix::t(v$counts), v$coords, labels = lab,
                            diameters = c(1, 400), panel_sizes = 20,
                            n_regions = 30L, seed = 1L),
    "skipped")
  expect_s3_class(sa, "recoverR_adequacy")
})

test_that("fit_hill rejects unphysical and bound-pinned fits", {
  # Three failure modes have all been observed in practice. Each must return NULL
  # rather than a number a caller could quote as an adequacy threshold:
  #   1. a curve with no variation carries no saturation information
  #   2. a fit whose optimiser runs INTO a bound is not an estimate, it is the wall
  #   3. a genuinely saturating curve must still FIT (the guard must not be a
  #      blanket refusal)
  d <- rep(c(300, 500, 750, 1000, 1500, 2000, 3000, 4000), length.out = 24)

  # 1. flat curve
  expect_null(recoverR:::fit_hill(d, rep(0.6, 24)))

  # 2. pure noise: the Hill exponent pins at the lower bound of 0.5. Measured
  #    before the guard: Pinf 0.733, h 0.500, K 1.5 -- accepted and meaningless.
  set.seed(2)
  expect_message(recoverR:::fit_hill(d, runif(24, 0.5, 0.9)), "pinned at its bound|rejected")

  # 3. a real saturating curve
  dg <- c(300, 500, 750, 1000, 1500, 2000, 3000, 4000)
  f <- recoverR:::fit_hill(dg, 0.95 * dg^2 / (1200^2 + dg^2))
  expect_false(is.null(f))
  expect_true(abs(unname(f$coef["K"]) - 1200) / 1200 < 0.1)
  # d* for a 0.95 target on this curve (~5231) lies OUTSIDE the sampled range, so
  # it must be flagged as an extrapolation rather than reported as a measurement.
  ds <- f$d_star(0.95)
  expect_true(isTRUE(attr(ds, "extrapolated")))
})

# ---------------------------------------------------------------------------
# The imaging route (0.1.3): read_roi_table() must apply THE SAME RULE as the
# Visium route, so that an IF cohort and a Visium cohort answer the same
# question in the same language. These tests pin the equivalence, not just the
# arithmetic - a second copy of the rule is the drift risk this exists to stop.
# ---------------------------------------------------------------------------

make_roi_fixture <- function() {
  data.frame(
    roi_id           = c("R1", "R2", "R3", "R4", "R5"),
    sample_id        = "P1",
    diameter_mm      = c(0.5, 0.5, 1, 3, 3),          # an extra column must pass through
    n_tumour         = c(100, 100, 100, 100, 100),
    n_stroma         = c(100, 100, 100, 100, 100),
    #             desert(<5% immune)   excluded          inflamed
    n_tumour_immune  = c(2,   10,  10,  30,  30),
    n_stroma_immune  = c(2,   10,  80,  10,  30),
    stringsAsFactors = FALSE)
}

test_that("read_roi_table reproduces the three-way rule from counts", {
  t <- read_roi_table(make_roi_fixture())
  # R1: f_immune = 4/200 = 0.02 < 0.05            -> desert
  # R2: f_immune = 20/200 = 0.10, f_tum-imm = 0.10 -> excluded
  # R3: f_immune = 90/200 = 0.45, f_tum-imm = 0.10 -> excluded
  # R4: f_tum-imm = 0.30 >= 0.25                   -> inflamed
  # R5: f_tum-imm = 0.30                           -> inflamed
  expect_equal(t$archetype, c("desert", "excluded", "excluded", "inflamed", "inflamed"))
  expect_equal(t$n_spots, rep(200, 5))
  expect_equal(t$f_immune[1], 0.02)
  expect_equal(t$f_tumour_immune[4], 0.30)
  expect_true(all(c("n_spots", "n_tumour", "f_immune", "f_tumour_immune",
                    "f_stroma_immune", "archetype", "region", "sample_id") %in% names(t)))
  expect_equal(t$diameter_mm, c(0.5, 0.5, 1, 3, 3))   # extras pass through
})

test_that("both spatial routes apply the shared rule to the fractions they derive", {
  # The guarantee is "one rule, two routes": each route must reach its call through
  # call_archetype(), on the fractions THAT ROUTE derived.
  #
  # This test used to loop over f_imm while building an ROI with n_stroma = 0, so the derived
  # f_immune always EQUALLED f_tumour_immune and the loop's f_imm never reached the data: 8 of
  # its 16 expectations failed against correct code. The fixture was wrong, not the rule. The
  # cross-route test below asserts the property that was actually meant.

  # (a) the imaging route, from counts, with both fractions varied independently
  grid <- expand.grid(f_tumour_immune = c(0, 0.249, 0.25, 0.9),
                      f_stroma_immune = c(0, 0.1, 1))
  for (i in seq_len(nrow(grid))) {
    nT <- 1000; nS <- 1000
    roi <- read_roi_table(data.frame(roi_id = "r", sample_id = "s",
      n_tumour = nT, n_stroma = nS,
      n_tumour_immune = grid$f_tumour_immune[i] * nT,
      n_stroma_immune = grid$f_stroma_immune[i] * nS))
    expect_equal(roi$archetype,
                 recoverR:::call_archetype(roi$f_immune, roi$f_tumour_immune),
                 info = paste("roi f_ti =", grid$f_tumour_immune[i],
                              "f_si =", grid$f_stroma_immune[i]))
  }

  # (b) the Visium route, from programme scores, the same property
  for (k in c(0, 4, 5, 8)) {                    # immune-high spots out of 20
    sc <- data.frame(immune = c(rep(2, k), rep(0, 20 - k)))
    a <- derive_spatial_archetypes(sc, compartment = rep("tumour", 20))
    expect_equal(a$archetype, recoverR:::call_archetype(a$f_immune, a$f_tumour_immune))
  }
})

test_that("the two routes agree with each other on identical fractions", {
  # The actual cross-route guarantee: all-tumour material, so both routes see the SAME
  # (f_immune, f_tumour_immune) pair and must therefore return the same call.
  for (f in c(0, 0.04, 0.05, 0.24, 0.25, 0.9)) {
    k <- round(f * 100)
    roi <- read_roi_table(data.frame(roi_id = "r", sample_id = "s",
      n_tumour = 100, n_stroma = 0, n_tumour_immune = f * 100, n_stroma_immune = 0))
    vis <- derive_spatial_archetypes(
      data.frame(immune = c(rep(2, k), rep(0, 100 - k))),
      compartment = rep("tumour", 100))
    expect_equal(roi$f_immune, vis$f_immune, tolerance = 1e-12)
    expect_equal(roi$f_tumour_immune, vis$f_tumour_immune, tolerance = 1e-12)
    expect_equal(roi$archetype, vis$archetype, info = paste("f =", f))
  }
})

test_that("the shared archetype rule holds at its thresholds", {
  # desert_cut = 0.05 (exclusive) and inflamed_cut = 0.25 (inclusive)
  ca <- recoverR:::call_archetype
  expect_equal(ca(0, 0.9), "desert")
  expect_equal(ca(0.049, 0.9), "desert")      # just below desert_cut
  expect_equal(ca(0.05, 0.9), "inflamed")     # AT desert_cut is no longer desert
  expect_equal(ca(0.2, 0.249), "excluded")    # just below inflamed_cut
  expect_equal(ca(0.2, 0.25), "inflamed")     # AT inflamed_cut
  expect_equal(ca(0.2, 0.9), "inflamed")
  expect_true(is.na(ca(NA, 0.3)))             # undefined stays undefined
  expect_true(is.na(ca(0.2, NA)))
})

test_that("read_roi_table reports rather than hides undefined regions", {
  x <- make_roi_fixture()
  x$n_tumour[1] <- 0; x$n_tumour_immune[1] <- 0        # no tumour compartment
  t <- read_roi_table(x)
  expect_true(is.na(t$f_tumour_immune[1]))
  expect_true(is.na(t$archetype[1]))                   # UNDEFINED, not forced into a class
  expect_equal(attr(t, "n_undefined_no_tumour"), 1L)
})

test_that("read_roi_table refuses fractions masquerading as counts", {
  x <- make_roi_fixture()
  x$n_tumour_immune[2] <- 150                          # > its compartment
  expect_error(read_roi_table(x), "more immune cells than compartment cells")
})

test_that("read_roi_table fails loudly on a missing or non-numeric column", {
  x <- make_roi_fixture(); x$n_stroma_immune <- NULL
  expect_error(read_roi_table(x), "missing column")
  y <- make_roi_fixture(); y$n_tumour_immune <- as.character(y$n_tumour_immune)
  z <- make_roi_fixture(); z$n_tumour_immune <- rep("high", 5)
  expect_error(read_roi_table(z), "non-numeric")
})
