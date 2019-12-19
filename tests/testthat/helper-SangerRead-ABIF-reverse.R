inputFilesPath <- system.file("extdata/", package = "sangeranalyseR")
SRab1RFN <- file.path(inputFilesPath,
                      "Allolobophora_chlorotica",
                      "ACHLO",
                      "ACHLO006-09[LCO1490_t1,HCO2198_t1]_R.ab1")
sangerReadR <- new("SangerRead",
                   readFeature           = "Reverse Read",
                   readFileName          = SRab1RFN,
                   geneticCode           = GENETIC_CODE,
                   TrimmingMethod        = "M1",
                   M1TrimmingCutoff      = 0.0001,
                   M2CutoffQualityScore  = NULL,
                   M2SlidingWindowSize   = NULL,
                   baseNumPerRow         = 100,
                   heightPerRow          = 200,
                   signalRatioCutoff     = 0.33,
                   showTrimmed           = TRUE)

