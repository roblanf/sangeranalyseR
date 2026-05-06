# =============================================================================
# Phase 8 — lazy-AA report compatibility.
#
# After Phase 6 introduced lazyAA = TRUE as the default, the RMarkdown
# templates and Shiny server files were silently reading `@primaryAASeqS*`
# slots directly — which return AAString("") on lazy objects, producing
# empty AA tables in reports. Phase 8 routes every read through the
# accessor methods (primaryAASeqS1/S2/S3) so reports populate correctly
# whether the user opted in to eager translation or not.
#
# These tests pin the lazy-default contract end-to-end.
# =============================================================================

inputFilesPath <- system.file("extdata/", package = "sangeranalyseR")
ab1_fwd <- file.path(inputFilesPath, "Allolobophora_chlorotica",
                     "ACHLO", "Achl_ACHLO006-09_1_F.ab1")

# -----------------------------------------------------------------------------
# 1. Accessor smoke: lazy-built SangerRead returns AA via accessor with
#    same content as eager-built SangerRead's slot.
# -----------------------------------------------------------------------------
test_that("lazy SangerRead: accessor result == eager slot for all 3 frames", {
    sr_lazy  <- new("SangerRead",
                    inputSource    = "ABIF",
                    readFeature    = "Forward Read",
                    readFileName   = ab1_fwd,
                    TrimmingMethod = "M1",
                    lazyAA         = TRUE)
    sr_eager <- new("SangerRead",
                    inputSource    = "ABIF",
                    readFeature    = "Forward Read",
                    readFileName   = ab1_fwd,
                    TrimmingMethod = "M1",
                    lazyAA         = FALSE)
    expect_identical(as.character(primaryAASeqS1(sr_lazy)),
                     as.character(sr_eager@primaryAASeqS1))
    expect_identical(as.character(primaryAASeqS2(sr_lazy)),
                     as.character(sr_eager@primaryAASeqS2))
    expect_identical(as.character(primaryAASeqS3(sr_lazy)),
                     as.character(sr_eager@primaryAASeqS3))
})

# -----------------------------------------------------------------------------
# 2. Lazy SangerRead has empty *AA slots but accessors yield non-empty
# -----------------------------------------------------------------------------
test_that("default lazy SangerRead has empty @primaryAASeq slots", {
    sr <- new("SangerRead",
              inputSource    = "ABIF",
              readFeature    = "Forward Read",
              readFileName   = ab1_fwd,
              TrimmingMethod = "M1")
    expect_equal(length(sr@primaryAASeqS1), 0L)
    expect_equal(length(sr@primaryAASeqS2), 0L)
    expect_equal(length(sr@primaryAASeqS3), 0L)
    # And that accessors still produce a well-formed AAString
    expect_s4_class(primaryAASeqS1(sr), "AAString")
    expect_gt(length(primaryAASeqS1(sr)), 0L)
})

# -----------------------------------------------------------------------------
# 3. Critical: generateReport() succeeds on the default lazyAA = TRUE path
#    (the regression that motivated Phase 8).
# -----------------------------------------------------------------------------
test_that("generateReport on lazy SangerRead produces a non-empty HTML file", {
    skip_on_cran()
    skip_if_not_installed("rmarkdown")
    skip_if_not(rmarkdown::pandoc_available(),
                "pandoc is required for rmarkdown::render")

    sr <- new("SangerRead",
              inputSource    = "ABIF",
              readFeature    = "Forward Read",
              readFileName   = ab1_fwd,
              TrimmingMethod = "M1",
              lazyAA         = TRUE)
    expect_true(sr@objectResults@creationResult)

    out_dir <- withr::local_tempdir()
    out <- tryCatch(
        suppressWarnings(suppressMessages(
            generateReportSR(sr, outputDir = out_dir)
        )),
        error = function(e) e
    )
    expect_false(inherits(out, "error"),
                 info = if (inherits(out, "error"))
                            paste("render error:", conditionMessage(out))
                        else "")

    # The Rmd writes the report to outputDir/SangerRead_Report/.
    # Find any generated .html under out_dir.
    html_files <- list.files(out_dir, pattern = "\\.html$",
                              recursive = TRUE, full.names = TRUE)
    expect_gt(length(html_files), 0L)
    if (length(html_files) > 0L) {
        # Non-empty file
        expect_gt(file.info(html_files[[1]])$size, 1000L)
    }
})

# -----------------------------------------------------------------------------
# 4. The internal RMD code paths (which we patched) produce non-empty AA
#    sequences for lazy reads. We can't import the .Rmd directly, but we
#    can simulate the (now-fixed) chunk: data.frame(AAString(primaryAASeqS1(sr)))
# -----------------------------------------------------------------------------
test_that("Phase-8 RMD pattern: data.frame(AAString(primaryAASeqS1(sr))) is non-empty", {
    sr <- new("SangerRead",
              inputSource    = "ABIF",
              readFeature    = "Forward Read",
              readFileName   = ab1_fwd,
              TrimmingMethod = "M1")
    df1 <- data.frame(Biostrings::AAString(primaryAASeqS1(sr)))
    df2 <- data.frame(Biostrings::AAString(primaryAASeqS2(sr)))
    df3 <- data.frame(Biostrings::AAString(primaryAASeqS3(sr)))
    # Pre-Phase-8: these would have been empty (zero rows).
    expect_gt(nrow(df1), 0L)
    expect_gt(nrow(df2), 0L)
    expect_gt(nrow(df3), 0L)
})
