# =============================================================================
# Phase 7 — peakvalues C++ port: mathematical-equivalence tests.
#
# `.peakvalues_r` (kept as a private helper in R/UtilitiesFunc.R) and
# `peakvalues_cpp` (Rcpp-exported, src/peakvalues.cpp) must return
# byte-identical NumericVectors for every input.
# =============================================================================

ns <- asNamespace("sangeranalyseR")

# Reach the private R helper:
peakvalues_r <- ns$.peakvalues_r

# Convenience: build the matrix shape that `getpeaks` produces.
#   col 1 = peak indexes (integer-valued doubles)
#   col 2 = trace amplitudes at those indexes
mk_peaks <- function(idxs, vals) cbind(as.numeric(idxs), as.numeric(vals))

# -----------------------------------------------------------------------------
# Edge cases enumerated in the Phase 7 plan
# -----------------------------------------------------------------------------

test_that("empty region: c(0, NA)", {
    x <- mk_peaks(c(10, 20, 30), c(5, 7, 6))
    r <- peakvalues_r(x,   100, 200)
    c <- peakvalues_cpp(x, 100, 200)
    expect_identical(r, c)
    expect_equal(c[[1]], 0)
    expect_true(is.na(c[[2]]))
})

test_that("zero-row matrix: c(0, NA)", {
    x <- matrix(numeric(0), nrow = 0, ncol = 2)
    r <- peakvalues_r(x,   0, 100)
    c <- peakvalues_cpp(x, 0, 100)
    expect_identical(r, c)
})

test_that("single peak in region", {
    x <- mk_peaks(c(10, 50, 90), c(3, 9, 4))
    r <- peakvalues_r(x,   30, 60)   # only idx=50 (val=9) qualifies
    c <- peakvalues_cpp(x, 30, 60)
    expect_identical(r, c)
    expect_equal(c, c(9, 50))
})

test_that("multiple peaks in region, distinct max", {
    x <- mk_peaks(c(10, 20, 30, 40, 50), c(2, 8, 5, 9, 3))
    r <- peakvalues_r(x,   15, 45)   # idxs 20,30,40 -> vals 8,5,9
    c <- peakvalues_cpp(x, 15, 45)
    expect_identical(r, c)
    expect_equal(c, c(9, 40))
})

test_that("tied max takes the FIRST occurrence (matches which.max)", {
    x <- mk_peaks(c(10, 20, 30, 40), c(5, 9, 9, 7))
    r <- peakvalues_r(x,   5, 50)    # both idx=20 and idx=30 share max=9
    c <- peakvalues_cpp(x, 5, 50)
    expect_identical(r, c)
    expect_equal(c, c(9, 20))
})

test_that("strict boundary inequality: pstart and pstop are excluded", {
    x <- mk_peaks(c(10, 20, 30), c(5, 7, 6))
    r <- peakvalues_r(x,   10, 30)   # excludes 10 and 30 -> only 20
    c <- peakvalues_cpp(x, 10, 30)
    expect_identical(r, c)
    expect_equal(c, c(7, 20))
})

test_that("NA in column 2 is skipped", {
    x <- mk_peaks(c(10, 20, 30), c(NA, 7, 4))
    r <- peakvalues_r(x,   5, 35)
    c <- peakvalues_cpp(x, 5, 35)
    expect_identical(r, c)
    expect_equal(c, c(7, 20))
})

test_that("non-integer trace values are preserved exactly", {
    x <- mk_peaks(c(11, 22, 33), c(3.14159, 2.71828, 1.61803))
    r <- peakvalues_r(x,   0, 100)
    c <- peakvalues_cpp(x, 0, 100)
    expect_identical(r, c)
})

# -----------------------------------------------------------------------------
# Batch API: peakvalues_batch_cpp matches a vectorised .peakvalues_r loop
# -----------------------------------------------------------------------------
test_that("peakvalues_batch_cpp[, j] == peakvalues_cpp(x, starts[j], stops[j])", {
    set.seed(20260503L)
    idxs   <- sort(sample.int(2000L, 200L))
    vals   <- runif(200L, 0, 1000)
    x      <- mk_peaks(idxs, vals)
    starts <- sort(sample(1:1900L, 50L, replace = TRUE))
    stops  <- starts + sample(20:80L, 50L, replace = TRUE)

    batch <- peakvalues_batch_cpp(x, starts, stops)
    expect_equal(dim(batch), c(2L, length(starts)))

    for (j in seq_along(starts)) {
        single <- peakvalues_cpp(x, starts[j], stops[j])
        expect_identical(batch[, j], single,
                          info = sprintf("window j=%d", j))
    }
})

test_that("peakvalues_batch_cpp[, j] matches .peakvalues_r(x, starts[j], stops[j])", {
    set.seed(11L)
    idxs   <- sort(sample.int(500L, 80L))
    vals   <- runif(80L, 0, 100)
    x      <- mk_peaks(idxs, vals)
    starts <- sort(sample(1:450L, 30L, replace = TRUE))
    stops  <- starts + 25L

    batch <- peakvalues_batch_cpp(x, starts, stops)
    for (j in seq_along(starts)) {
        r <- peakvalues_r(x, starts[j], stops[j])
        expect_identical(batch[, j], r,
                          info = sprintf("window j=%d", j))
    }
})

test_that("peakvalues_batch_cpp validates pstarts / pstops length match", {
    x <- mk_peaks(c(10, 20, 30), c(5, 7, 6))
    expect_error(peakvalues_batch_cpp(x, c(1, 2), c(3)),
                  regexp = "equal length")
})

# -----------------------------------------------------------------------------
# Larger random-input fuzz: identical results across 200 random configs
# -----------------------------------------------------------------------------
test_that("fuzz: 200 random matrices give identical R / C++ results", {
    set.seed(20260503L)
    mismatches <- 0L
    for (trial in seq_len(200L)) {
        n <- sample(1:50, 1L)
        idxs <- sort(sample.int(2000L, n, replace = FALSE))
        vals <- runif(n, min = 0, max = 1000)
        x <- mk_peaks(idxs, vals)
        pstart <- runif(1L, 0, 2000)
        pstop  <- pstart + runif(1L, 1, 500)
        r <- peakvalues_r(x,   pstart, pstop)
        c <- peakvalues_cpp(x, pstart, pstop)
        if (!identical(r, c)) mismatches <- mismatches + 1L
    }
    expect_equal(mismatches, 0L)
})

# -----------------------------------------------------------------------------
# End-to-end determinism: SangerAlignment consensus on the ACHLO fixture
# must equal the Phase-6 result (locked in by recomputing twice).
# -----------------------------------------------------------------------------
test_that("SangerAlignment consensus is deterministic across two C++-driven runs", {
    ab1_dir <- system.file("extdata", "Allolobophora_chlorotica", "ACHLO",
                           package = "sangeranalyseR")
    args <- list(
        inputSource         = "ABIF",
        processMethod       = "REGEX",
        ABIF_Directory      = ab1_dir,
        REGEX_SuffixForward = "_[0-9]*_F.ab1$",
        REGEX_SuffixReverse = "_[0-9]*_R.ab1$",
        TrimmingMethod      = "M1",
        M1TrimmingCutoff    = 0.0001,
        BPPARAM             = BiocParallel::SerialParam()
    )
    sa1 <- do.call(SangerAlignment, args)
    sa2 <- do.call(SangerAlignment, args)
    expect_true(sa1@objectResults@creationResult)
    expect_true(sa2@objectResults@creationResult)
    expect_identical(as.character(sa1@contigsConsensus),
                     as.character(sa2@contigsConsensus))
    expect_identical(sort(names(sa1@contigList)),
                     sort(names(sa2@contigList)))
})

# -----------------------------------------------------------------------------
# Single-read constructed via R-only path vs C++-driven path: identical
# numeric output. (This is *the* invariant — base calling must not drift.)
# -----------------------------------------------------------------------------
test_that("MakeBaseCallsInside: R-loop and C++-loop produce identical primarySeq", {
    ab1 <- system.file("extdata", "Allolobophora_chlorotica", "ACHLO",
                       "Achl_ACHLO006-09_1_F.ab1", package = "sangeranalyseR")
    sr <- new("SangerRead",
              inputSource    = "ABIF",
              readFeature    = "Forward Read",
              readFileName   = ab1,
              TrimmingMethod = "M1")

    # Re-run MakeBaseCallsInside under a temporary monkey-patch that swaps
    # peakvalues_cpp for the R helper. Both invocations must produce the
    # same primarySeq / qualityPhredScores.
    abif <- sangerseqR::read.abif(ab1)
    ss   <- sangerseqR::sangerseq(abif)

    # The `peakvalues_cpp` symbol lives in the namespace; replace via a
    # local binding inside the package env. testthat::with_mocked_bindings
    # is the supported entry point.
    skip_if_not_installed("testthat")

    res_cpp <- ns$MakeBaseCallsInside(ss@traceMatrix,
                                       ss@peakPosMatrix,
                                       abif@data$PCON.2,
                                       0.33, "Forward Read", "")
    res_r <- testthat::with_mocked_bindings(
        peakvalues_cpp = ns$.peakvalues_r,
        .package       = "sangeranalyseR",
        code           = ns$MakeBaseCallsInside(ss@traceMatrix,
                                                ss@peakPosMatrix,
                                                abif@data$PCON.2,
                                                0.33, "Forward Read", "")
    )

    expect_identical(as.character(res_cpp$primarySeq),
                     as.character(res_r$primarySeq))
    expect_identical(as.character(res_cpp$secondarySeq),
                     as.character(res_r$secondarySeq))
    expect_identical(res_cpp$qualityPhredScores, res_r$qualityPhredScores)
})
