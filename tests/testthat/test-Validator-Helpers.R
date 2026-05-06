# =============================================================================
# Phase 4 — Direct unit tests for the new internal helpers in
# R/UtilitiesFuncInputChecker.R (.errAppend, .requireType, .requireEnum,
# .requireRange, .requireExt). Helpers are dot-prefixed and not exported, so
# we reach them via `:::`.
#
# These tests guard the consolidation work — if a future refactor breaks the
# helper contracts, every public `check*` wrapper that depends on them will
# also fail; these tests narrow the diagnosis.
# =============================================================================

ns <- asNamespace("sangeranalyseR")

test_that(".errAppend grows both vectors in lockstep", {
    out <- ns$.errAppend(character(0), character(0), "msg1", "T1")
    expect_equal(out[[1]], "msg1")
    expect_equal(out[[2]], "T1")

    out <- ns$.errAppend(out[[1]], out[[2]], "msg2", "T2")
    expect_equal(out[[1]], c("msg1", "msg2"))
    expect_equal(out[[2]], c("T1", "T2"))
})

test_that(".requireType returns inputs unchanged on success", {
    out <- ns$.requireType(123, "x", is.numeric, "numeric",
                            character(0), character(0))
    expect_length(out[[1]], 0L)
    expect_length(out[[2]], 0L)
})

test_that(".requireType appends an error on predicate failure", {
    out <- ns$.requireType("abc", "x", is.numeric, "numeric",
                            character(0), character(0))
    expect_length(out[[1]], 1L)
    expect_match(out[[1]], "must be numeric")
    expect_equal(out[[2]], "PARAMETER_TYPE_ERROR")
})

test_that(".requireEnum accepts allowed values", {
    out <- ns$.requireEnum("ABIF", "inputSource", c("ABIF", "FASTA"),
                            character(0), character(0))
    expect_length(out[[1]], 0L)
})

test_that(".requireEnum rejects values outside the allowed set", {
    out <- ns$.requireEnum("FAST5", "inputSource", c("ABIF", "FASTA"),
                            character(0), character(0))
    expect_length(out[[1]], 1L)
    expect_equal(out[[2]], "PARAMETER_VALUE_ERROR")
})

test_that(".requireEnum rejects non-character input", {
    out <- ns$.requireEnum(1L, "inputSource", c("ABIF", "FASTA"),
                            character(0), character(0))
    expect_length(out[[1]], 1L)
    expect_match(out[[1]], "must be character")
})

test_that(".requireRange accepts in-range numerics", {
    out <- ns$.requireRange(0.5, "x", 0, 1, character(0), character(0))
    expect_length(out[[1]], 0L)
})

test_that(".requireRange rejects out-of-range high", {
    out <- ns$.requireRange(1.5, "x", 0, 1, character(0), character(0))
    expect_equal(out[[2]], "PARAMETER_RANGE_ERROR")
})

test_that(".requireRange rejects out-of-range low", {
    out <- ns$.requireRange(-1, "x", 0, 1, character(0), character(0))
    expect_equal(out[[2]], "PARAMETER_RANGE_ERROR")
})

test_that(".requireRange rejects non-integer when integer=TRUE", {
    out <- ns$.requireRange(1.5, "x", 0, 10, character(0), character(0),
                             integer = TRUE)
    expect_equal(out[[2]], "PARAMETER_TYPE_ERROR")
})

test_that(".requireRange rejects non-numeric input", {
    out <- ns$.requireRange("abc", "x", 0, 1, character(0), character(0))
    expect_equal(out[[2]], "PARAMETER_TYPE_ERROR")
})

test_that(".requireExt accepts files matching the (escaped) regex", {
    out <- ns$.requireExt("/tmp/some_F.ab1", "\\.ab1$",
                           character(0), character(0),
                           message = "bad extension")
    expect_length(out[[1]], 0L)
})

test_that(".requireExt rejects .Xab1 (the original regex bug)", {
    out <- ns$.requireExt("/tmp/some_F.Xab1", "\\.ab1$",
                           character(0), character(0),
                           message = "bad extension")
    expect_length(out[[1]], 1L)
    expect_equal(out[[2]], "FILE_TYPE_ERROR")
})

test_that(".requireExt accepts both .fa and .fasta with the FASTA pattern", {
    out1 <- ns$.requireExt("/tmp/x.fa", "\\.(fa|fasta)$",
                            character(0), character(0),
                            message = "bad")
    out2 <- ns$.requireExt("/tmp/x.fasta", "\\.(fa|fasta)$",
                            character(0), character(0),
                            message = "bad")
    expect_length(out1[[1]], 0L)
    expect_length(out2[[1]], 0L)
})

test_that(".requireExt rejects .Xfa with the FASTA pattern", {
    out <- ns$.requireExt("/tmp/x.Xfa", "\\.(fa|fasta)$",
                           character(0), character(0),
                           message = "bad")
    expect_length(out[[1]], 1L)
})

# -----------------------------------------------------------------------------
# Public-API smoke: every consolidated `check*` wrapper still returns the
# documented (errors, errorTypes) shape and the right tag on failure.
# -----------------------------------------------------------------------------

test_that("checkInputSource accepts ABIF/FASTA, rejects others", {
    expect_length(checkInputSource("ABIF", character(0), character(0))[[1]], 0L)
    expect_length(checkInputSource("FASTA", character(0), character(0))[[1]], 0L)
    out <- checkInputSource("FAST5", character(0), character(0))
    expect_equal(out[[2]], "PARAMETER_VALUE_ERROR")
})

test_that("checkProcessMethod accepts REGEX/CSV, rejects others", {
    expect_length(checkProcessMethod("ABIF", "REGEX",
                                      character(0), character(0))[[1]], 0L)
    expect_length(checkProcessMethod("ABIF", "CSV",
                                      character(0), character(0))[[1]], 0L)
    out <- checkProcessMethod("ABIF", "GLOB", character(0), character(0))
    expect_equal(out[[2]], "PARAMETER_VALUE_ERROR")
})

test_that("checkSignalRatioCutoff range guard", {
    expect_length(checkSignalRatioCutoff(0.5, character(0), character(0))[[1]], 0L)
    out <- checkSignalRatioCutoff(1.5, character(0), character(0))
    expect_equal(out[[2]], "PARAMETER_RANGE_ERROR")
})

test_that("checkProcessorsNum accepts NULL", {
    out <- checkProcessorsNum(NULL, character(0), character(0))
    expect_length(out[[1]], 0L)
})

test_that("checkProcessorsNum rejects non-numeric", {
    out <- checkProcessorsNum("eight", character(0), character(0))
    expect_equal(out[[2]], "PARAMETER_TYPE_ERROR")
})

test_that("checkReadingFrame accepts 1/2/3, rejects 4", {
    for (rf in 1:3) {
        out <- checkReadingFrame(rf, character(0), character(0))
        expect_length(out[[1]], 0L)
    }
    out <- checkReadingFrame(4, character(0), character(0))
    expect_equal(out[[2]], "PARAMETER_VALUE_ERROR")
})

test_that("checkMinReadsNum rejects 0 with PARAMETER_VALUE_ERROR", {
    out <- checkMinReadsNum(0, character(0), character(0))
    expect_true("PARAMETER_VALUE_ERROR" %in% out[[2]])
})

test_that("checkMinFractionCall / checkMaxFractionLost in [0,1]", {
    expect_length(checkMinFractionCall(0.5, character(0), character(0))[[1]], 0L)
    expect_equal(checkMinFractionCall(1.5, character(0), character(0))[[2]],
                 "PARAMETER_RANGE_ERROR")
    expect_length(checkMaxFractionLost(0.5, character(0), character(0))[[1]], 0L)
    expect_equal(checkMaxFractionLost(-1, character(0), character(0))[[2]],
                 "PARAMETER_RANGE_ERROR")
})
