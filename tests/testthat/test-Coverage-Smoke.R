# =============================================================================
# Phase 3 — Coverage breadth smoke tests.
#
# Tickle code paths that no other test file exercises so the coverage report
# isn't blind to large swaths of MethodShared, MethodSangerRead, etc.
#
# Skip Shiny launchApp* (heavy / interactive) and generateReport* (writes to
# disk; covered by examples).
# =============================================================================

test_that("data() loads the four bundled fixtures", {
    expect_silent(utils::data("qualityReportData", package = "sangeranalyseR"))
    expect_silent(utils::data("sangerReadFData",   package = "sangeranalyseR"))
    expect_silent(utils::data("sangerContigData",  package = "sangeranalyseR"))
    expect_silent(utils::data("sangerAlignmentData", package = "sangeranalyseR"))
})

test_that("qualityBasePlot() returns a plotly object for QualityReport and SangerRead", {
    utils::data("qualityReportData", package = "sangeranalyseR",
                envir = environment())
    utils::data("sangerReadFData",   package = "sangeranalyseR",
                envir = environment())

    p1 <- qualityBasePlot(qualityReportData)
    p2 <- qualityBasePlot(sangerReadFData)
    expect_true(inherits(p1, "plotly") || inherits(p1, "htmlwidget"))
    expect_true(inherits(p2, "plotly") || inherits(p2, "htmlwidget"))
})

test_that("MakeBaseCalls() returns a SangerRead", {
    utils::data("sangerReadFData", package = "sangeranalyseR",
                envir = environment())
    bc <- MakeBaseCalls(sangerReadFData, signalRatioCutoff = 0.22)
    expect_s4_class(bc, "SangerRead")
})

test_that("readTable() runs without error for SangerRead/Contig/Alignment", {
    utils::data("sangerReadFData",   package = "sangeranalyseR",
                envir = environment())
    utils::data("sangerContigData",  package = "sangeranalyseR",
                envir = environment())
    utils::data("sangerAlignmentData", package = "sangeranalyseR",
                envir = environment())

    expect_invisible(invisible(readTable(sangerReadFData)))
    # readTable for SC/SA prints to stdout via log_info; just ensure no crash.
    expect_error(readTable(sangerContigData),    NA)
    expect_error(readTable(sangerAlignmentData), NA)
})

test_that("writeFasta() dispatches for SangerRead, SangerContig, SangerAlignment", {
    utils::data("sangerReadFData",   package = "sangeranalyseR",
                envir = environment())
    utils::data("sangerContigData",  package = "sangeranalyseR",
                envir = environment())
    utils::data("sangerAlignmentData", package = "sangeranalyseR",
                envir = environment())

    out <- withr::local_tempdir()
    expect_error(writeFasta(sangerReadFData,    outputDir = out), NA)
    expect_error(writeFasta(sangerContigData,   outputDir = out), NA)
    expect_error(writeFasta(sangerAlignmentData, outputDir = out), NA)

    written <- list.files(out, recursive = TRUE)
    expect_gt(length(written), 0L)
})

test_that("updateQualityParam() re-trims a SangerRead", {
    utils::data("sangerReadFData", package = "sangeranalyseR",
                envir = environment())
    sr2 <- updateQualityParam(sangerReadFData,
                              TrimmingMethod        = "M2",
                              M1TrimmingCutoff      = NULL,
                              M2CutoffQualityScore  = 30,
                              M2SlidingWindowSize   = 12)
    expect_s4_class(sr2, "SangerRead")
    expect_equal(sr2@QualityReport@TrimmingMethod, "M2")
    expect_equal(sr2@QualityReport@M2CutoffQualityScore, 30)
    expect_equal(sr2@QualityReport@M2SlidingWindowSize,  12)
})
