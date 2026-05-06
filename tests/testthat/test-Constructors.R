# =============================================================================
# Phase 3 — Constructor smoke tests.
#
# Quick tests for the three exported wrapper functions in R/Constructors.R
# (SangerRead(), SangerContig(), SangerAlignment()), confirming each minimum-
# args invocation produces a valid S4 object whose ObjectResults reports
# successful creation.
# =============================================================================

inputFilesPath <- system.file("extdata/", package = "sangeranalyseR")
ab1_forward    <- file.path(inputFilesPath, "Allolobophora_chlorotica",
                            "ACHLO", "Achl_ACHLO006-09_1_F.ab1")
ab1_dir        <- file.path(inputFilesPath, "Allolobophora_chlorotica", "ACHLO")

test_that("SangerRead() wrapper builds a valid SangerRead", {
    sr <- SangerRead(inputSource    = "ABIF",
                     readFeature    = "Forward Read",
                     readFileName   = ab1_forward,
                     TrimmingMethod = "M1")
    expect_s4_class(sr, "SangerRead")
    expect_true(sr@objectResults@creationResult)
    expect_equal(sr@inputSource, "ABIF")
    expect_equal(sr@readFeature, "Forward Read")
})

test_that("SangerContig() wrapper builds a valid SangerContig", {
    sc <- SangerContig(inputSource         = "ABIF",
                       processMethod       = "REGEX",
                       ABIF_Directory      = ab1_dir,
                       contigName          = "Achl_ACHLO006-09",
                       REGEX_SuffixForward = "_[0-9]*_F.ab1$",
                       REGEX_SuffixReverse = "_[0-9]*_R.ab1$",
                       TrimmingMethod      = "M1",
                       processorsNum       = 1)
    expect_s4_class(sc, "SangerContig")
    expect_true(sc@objectResults@creationResult)
    expect_equal(sc@inputSource, "ABIF")
    expect_equal(sc@processMethod, "REGEX")
})

test_that("SangerAlignment() wrapper builds a valid SangerAlignment", {
    sa <- SangerAlignment(inputSource         = "ABIF",
                          processMethod       = "REGEX",
                          ABIF_Directory      = ab1_dir,
                          REGEX_SuffixForward = "_[0-9]*_F.ab1$",
                          REGEX_SuffixReverse = "_[0-9]*_R.ab1$",
                          TrimmingMethod      = "M1",
                          processorsNum       = 1)
    expect_s4_class(sa, "SangerAlignment")
    expect_true(sa@objectResults@creationResult)
    expect_gt(length(sa@contigList), 0L)
})
