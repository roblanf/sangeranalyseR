inputFilesPath <- system.file("extdata/", package = "sangeranalyseR")
SRab1FFN <- file.path(inputFilesPath,
                      "Allolobophora_chlorotica",
                      "ACHLO",
                      "Achl_ACHLO006-09_1_F.ab1")
sangerReadF <- new("SangerRead",
                   readFeature           = "Forward Read",
                   readFileName          = SRab1FFN,
                   geneticCode           = GENETIC_CODE,
                   TrimmingMethod        = "M1",
                   M1TrimmingCutoff      = 0.0001,
                   M2CutoffQualityScore  = NULL,
                   M2SlidingWindowSize   = NULL,
                   baseNumPerRow         = 100,
                   heightPerRow          = 200,
                   signalRatioCutoff     = 0.33,
                   showTrimmed           = TRUE)

