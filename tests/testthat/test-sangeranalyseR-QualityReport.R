test_that("QualityReport update quality trimming parameters 1 (Smaller M1TrimmingCutoff)", {
    qualityReport <- sangerRead@QualityReport
    newQualityReport <- updateQualityParam(qualityReport,
                                           TrimmingMethod         = "M1",
                                           M1TrimmingCutoff       = 0.000005,
                                           M2CutoffQualityScore   = NULL,
                                           M2SlidingWindowSize    = NULL)
    expect_type(newQualityReport, "S4")
    expect_s4_class(newQualityReport, "QualityReport")
    expect_equal(newQualityReport@TrimmingMethod, "M1")
    expect_equal(newQualityReport@M1TrimmingCutoff, 0.000005, tolerance=1e-10)
    expect_equal(newQualityReport@M2CutoffQualityScore, NULL)
    expect_equal(newQualityReport@M2SlidingWindowSize, NULL)
    expect_equal(length(newQualityReport@qualityPhredScoresRaw), 703)
    expect_equal(length(newQualityReport@qualityPhredScores), 702)
    expect_equal(length(newQualityReport@qualityBaseScores), 702)

    expect_equal(newQualityReport@rawSeqLength, 702)
    expect_equal(newQualityReport@trimmedSeqLength, 461)
    expect_equal(newQualityReport@trimmedStartPos, 16)
    expect_equal(newQualityReport@trimmedFinishPos, 477)
    expect_equal(newQualityReport@rawMeanQualityScore, 52.87607, tolerance=1e-6)
    expect_equal(newQualityReport@trimmedMeanQualityScore, 58.21041, tolerance=1e-6)
    expect_equal(newQualityReport@rawMinQualityScore, 1)
    expect_equal(newQualityReport@trimmedMinQualityScore, 13)
    expect_equal(newQualityReport@remainingRatio, 0.6566952, tolerance=1e-6)
})

test_that("QualityReport update quality trimming parameters 2 (bigger M1TrimmingCutoff)", {
    qualityReport <- sangerRead@QualityReport
    newQualityReport <- updateQualityParam(qualityReport,
                                           TrimmingMethod         = "M1",
                                           M1TrimmingCutoff       = 0.4,
                                           M2CutoffQualityScore   = NULL,
                                           M2SlidingWindowSize    = NULL)
    expect_type(newQualityReport, "S4")
    expect_s4_class(newQualityReport, "QualityReport")
    expect_equal(newQualityReport@TrimmingMethod, "M1")
    expect_equal(newQualityReport@M1TrimmingCutoff, 0.4, tolerance=1e-10)
    expect_equal(newQualityReport@M2CutoffQualityScore, NULL)
    expect_equal(newQualityReport@M2SlidingWindowSize, NULL)
    expect_equal(length(newQualityReport@qualityPhredScoresRaw), 703)
    expect_equal(length(newQualityReport@qualityPhredScores), 702)
    expect_equal(length(newQualityReport@qualityBaseScores), 702)

    expect_equal(newQualityReport@rawSeqLength, 702)
    expect_equal(newQualityReport@trimmedSeqLength, 700)
    expect_equal(newQualityReport@trimmedStartPos, 2)
    expect_equal(newQualityReport@trimmedFinishPos, 702)
    expect_equal(newQualityReport@rawMeanQualityScore, 52.87607, tolerance=1e-6)
    expect_equal(newQualityReport@trimmedMeanQualityScore, 53.01857, tolerance=1e-6)
    expect_equal(newQualityReport@rawMinQualityScore, 1)
    expect_equal(newQualityReport@trimmedMinQualityScore, 1)
    expect_equal(newQualityReport@remainingRatio, 0.997151, tolerance=1e-6)
})

test_that("QualityReport update quality trimming parameters 3 (to M2)", {
    qualityReport <- sangerRead@QualityReport
    newQualityReport <- updateQualityParam(qualityReport,
                                           TrimmingMethod         = "M2",
                                           M1TrimmingCutoff       = NULL,
                                           M2CutoffQualityScore   = 30,
                                           M2SlidingWindowSize    = 10)
    expect_type(newQualityReport, "S4")
    expect_s4_class(newQualityReport, "QualityReport")
    expect_equal(newQualityReport@TrimmingMethod, "M2")
    expect_equal(newQualityReport@M1TrimmingCutoff, NULL)
    expect_equal(newQualityReport@M2CutoffQualityScore, 30)
    expect_equal(newQualityReport@M2SlidingWindowSize, 10)
    expect_equal(length(newQualityReport@qualityPhredScoresRaw), 703)
    expect_equal(length(newQualityReport@qualityPhredScores), 702)
    expect_equal(length(newQualityReport@qualityBaseScores), 702)

    expect_equal(newQualityReport@rawSeqLength, 702)
    expect_equal(newQualityReport@trimmedSeqLength, 665)
    expect_equal(newQualityReport@trimmedStartPos, 19)
    expect_equal(newQualityReport@trimmedFinishPos, 684)
    expect_equal(newQualityReport@rawMeanQualityScore, 52.87607, tolerance=1e-6)
    expect_equal(newQualityReport@trimmedMeanQualityScore, 54.74887, tolerance=1e-6)
    expect_equal(newQualityReport@rawMinQualityScore, 1)
    expect_equal(newQualityReport@trimmedMinQualityScore, 12)
    expect_equal(newQualityReport@remainingRatio, 0.9472934, tolerance=1e-6)
})

