# =============================================================================
# Phase 4 — setValidity tests for ChromatogramParam, QualityReport,
# ObjectResults.
#
# These post-construction invariants only fire on direct `validObject()`
# calls or on `slot<-` assignments after construction; the package's
# Sanger* `initialize` methods produce only valid objects, so the suite is
# unaffected on the happy path.
# =============================================================================

# -----------------------------------------------------------------------------
# ChromatogramParam
# -----------------------------------------------------------------------------
test_that("ChromatogramParam: valid construction passes validObject", {
    cp <- new("ChromatogramParam",
              baseNumPerRow     = 100,
              heightPerRow      = 200,
              signalRatioCutoff = 0.33,
              showTrimmed       = TRUE)
    expect_true(validObject(cp, test = TRUE) == TRUE)
})

test_that("ChromatogramParam: out-of-range baseNumPerRow rejected by validity", {
    cp <- new("ChromatogramParam",
              baseNumPerRow     = 100,
              heightPerRow      = 200,
              signalRatioCutoff = 0.33,
              showTrimmed       = TRUE)
    cp@baseNumPerRow <- 999
    res <- validObject(cp, test = TRUE)
    expect_false(isTRUE(res))
    expect_match(paste(res, collapse = " "), "baseNumPerRow")
})

test_that("ChromatogramParam: out-of-range signalRatioCutoff rejected", {
    cp <- new("ChromatogramParam",
              baseNumPerRow     = 100,
              heightPerRow      = 200,
              signalRatioCutoff = 0.33,
              showTrimmed       = TRUE)
    cp@signalRatioCutoff <- 1.5
    res <- validObject(cp, test = TRUE)
    expect_false(isTRUE(res))
    expect_match(paste(res, collapse = " "), "signalRatioCutoff")
})

test_that("ChromatogramParam: non-logical showTrimmed rejected", {
    cp <- new("ChromatogramParam",
              baseNumPerRow     = 100,
              heightPerRow      = 200,
              signalRatioCutoff = 0.33,
              showTrimmed       = TRUE)
    # Direct slot assignment of an invalid type is itself blocked by S4 type
    # checking, so we test by pre-corrupting the prototype: build, mutate
    # to the wrong-length-logical case.
    cp@showTrimmed <- logical(0)
    res <- validObject(cp, test = TRUE)
    expect_false(isTRUE(res))
})

# -----------------------------------------------------------------------------
# QualityReport
# -----------------------------------------------------------------------------
test_that("QualityReport: empty default object passes validity", {
    # The empty/default state must be valid so vignettes can `new("QualityReport")`.
    qr <- new("QualityReport",
              qualityPhredScores = c(40L, 40L, 40L, 30L, 25L),
              TrimmingMethod     = "M1",
              M1TrimmingCutoff   = 0.0001)
    expect_true(validObject(qr, test = TRUE) == TRUE)
})

test_that("QualityReport: trimmedFinishPos < trimmedStartPos rejected", {
    qr <- new("QualityReport",
              qualityPhredScores = c(40L, 40L, 40L, 30L, 25L),
              TrimmingMethod     = "M1",
              M1TrimmingCutoff   = 0.0001)
    qr@trimmedStartPos  <- 4
    qr@trimmedFinishPos <- 1
    res <- validObject(qr, test = TRUE)
    expect_false(isTRUE(res))
    expect_match(paste(res, collapse = " "), "trimmedFinishPos")
})

test_that("QualityReport: remainingRatio outside [0,1] rejected", {
    qr <- new("QualityReport",
              qualityPhredScores = c(40L, 40L, 40L, 30L, 25L),
              TrimmingMethod     = "M1",
              M1TrimmingCutoff   = 0.0001)
    qr@remainingRatio <- 1.5
    res <- validObject(qr, test = TRUE)
    expect_false(isTRUE(res))
    expect_match(paste(res, collapse = " "), "remainingRatio")
})

# -----------------------------------------------------------------------------
# ObjectResults
# -----------------------------------------------------------------------------
test_that("ObjectResults: valid default passes validity", {
    or <- new("ObjectResults")
    expect_true(validObject(or, test = TRUE) == TRUE)
})

test_that("ObjectResults: mismatched errorMessages/errorTypes lengths rejected", {
    or <- new("ObjectResults")
    or@errorMessages <- c("msg1", "msg2", "msg3")
    or@errorTypes    <- c("T1")
    res <- validObject(or, test = TRUE)
    expect_false(isTRUE(res))
    expect_match(paste(res, collapse = " "), "errorMessages")
})

test_that("ObjectResults: mismatched warning lengths rejected", {
    or <- new("ObjectResults")
    or@warningMessages <- c("a", "b")
    or@warningTypes    <- character(0)
    res <- validObject(or, test = TRUE)
    expect_false(isTRUE(res))
    expect_match(paste(res, collapse = " "), "warning")
})
