test_that("SangerRead creation test", {
    inputFilesPath <- system.file("extdata/", package = "sangeranalyseR")
    A_chloroticaFdReadFN <-
        file.path(inputFilesPath, "Allolobophora_chlorotica",
                  "RBNII396-13[C_LepFolF,C_LepFolR]_F_1.ab1")
    A_chloroticaRead <- SangerRead(
                            readFeature           = "Forward Read",
                            readFileName          = A_chloroticaFdReadFN,
                            geneticCode           = GENETIC_CODE,
                            TrimmingMethod        = "M1",
                            M1TrimmingCutoff      = 0.0001,
                            M2CutoffQualityScore  = NULL,
                            M2SlidingWindowSize   = NULL,
                            baseNumPerRow         = 80,
                            heightPerRow          = 200,
                            signalRatioCutoff     = 0.33,
                            showTrimmed           = TRUE)

    # expect_type(A_chloroticaRead, "S4")
    expect_equal(A_chloroticaRead@readFeature, "Forward Read")
    # expect_s4_class(A_chloroticaRead, "SangerRead")
})
