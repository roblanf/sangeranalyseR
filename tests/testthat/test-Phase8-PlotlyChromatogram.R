# =============================================================================
# Phase 8 — chromatogram_plotly correctness + downsampling.
# =============================================================================

inputFilesPath <- system.file("extdata/", package = "sangeranalyseR")
ab1_fwd <- file.path(inputFilesPath, "Allolobophora_chlorotica",
                     "ACHLO", "Achl_ACHLO006-09_1_F.ab1")

test_that("chromatogram_plotly returns a plotly htmlwidget", {
    sr <- new("SangerRead",
              inputSource    = "ABIF",
              readFeature    = "Forward Read",
              readFileName   = ab1_fwd,
              TrimmingMethod = "M1")
    p <- chromatogram_plotly(sr)
    expect_true(inherits(p, "plotly") || inherits(p, "htmlwidget"))
})

test_that("chromatogram_plotly downsamples when trace exceeds max_points", {
    sr <- new("SangerRead",
              inputSource    = "ABIF",
              readFeature    = "Forward Read",
              readFileName   = ab1_fwd,
              TrimmingMethod = "M1")
    n_total <- nrow(sr@traceMatrix)
    if (n_total > 1000L) {
        p <- chromatogram_plotly(sr, max_points = 500L)
        info <- attr(p, "downsample_info")
        expect_lte(info$rendered_points, 500L + 1L)
        expect_equal(info$original_points, n_total)
        expect_gte(info$downsample_stride, 2L)
    } else {
        skip("trace too short to test downsampling")
    }
})

test_that("chromatogram_plotly preserves all points below the threshold", {
    sr <- new("SangerRead",
              inputSource    = "ABIF",
              readFeature    = "Forward Read",
              readFileName   = ab1_fwd,
              TrimmingMethod = "M1")
    n_total <- nrow(sr@traceMatrix)
    p <- chromatogram_plotly(sr, max_points = n_total + 100L)
    info <- attr(p, "downsample_info")
    expect_equal(info$rendered_points, n_total)
    expect_equal(info$downsample_stride, 1L)
})

test_that("chromatogram_plotly accepts cb_friendly and custom palettes", {
    sr <- new("SangerRead",
              inputSource    = "ABIF",
              readFeature    = "Forward Read",
              readFileName   = ab1_fwd,
              TrimmingMethod = "M1")
    p1 <- chromatogram_plotly(sr, colors = "cb_friendly")
    expect_true(inherits(p1, "plotly") || inherits(p1, "htmlwidget"))

    custom <- c("#ff0000", "#00ff00", "#0000ff", "#ffff00", "#000000")
    p2 <- chromatogram_plotly(sr, colors = custom)
    expect_true(inherits(p2, "plotly") || inherits(p2, "htmlwidget"))
})

test_that("chromatogram_plotly rejects bad inputs", {
    expect_error(chromatogram_plotly(list()),
                  regexp = "sangerseq")
    sr <- new("SangerRead",
              inputSource    = "ABIF",
              readFeature    = "Forward Read",
              readFileName   = ab1_fwd,
              TrimmingMethod = "M1")
    expect_error(chromatogram_plotly(sr, colors = "no_such_palette"),
                  regexp = "default|cb_friendly|character")
})
