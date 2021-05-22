inputFilesPath <- system.file("extdata/", package = "sangeranalyseR")
# 1. FILE_NOT_EXISTS_ERROR
SRab1_File_Not_Exist_Error <- file.path(inputFilesPath,
                                        "Allolobophora_chlorotica",
                                        "ACHLO",
                                        "TEST_Achl_ACHLO006-09_2_R.ab1")
sangerReadR_File_Not_Exist_Error <- new("SangerRead",
                                        readFeature           = "Reverse Read",
                                        readFileName          = SRab1_File_Not_Exist_Error,
                                        geneticCode           = GENETIC_CODE,
                                        TrimmingMethod        = "M1",
                                        M1TrimmingCutoff      = 0.0001,
                                        M2CutoffQualityScore  = NULL,
                                        M2SlidingWindowSize   = NULL,
                                        baseNumPerRow         = 100,
                                        heightPerRow          = 200,
                                        signalRatioCutoff     = 0.33,
                                        showTrimmed           = TRUE)


# 2. FILE_TYPE_ERROR: ab1 file
SRab1_FILE_TYPE_ERROR <- file.path(inputFilesPath,
                                   "Allolobophora_chlorotica",
                                   "ACHLO",
                                   "Achl_ACHLO006-09_1_F.ab2")
sangerReadF_FILE_TYPE_ERROR <- new("SangerRead",
                                   readFeature           = "Forward Read",
                                   readFileName          = SRab1_FILE_TYPE_ERROR,
                                   geneticCode           = GENETIC_CODE,
                                   TrimmingMethod        = "M1",
                                   M1TrimmingCutoff      = 0.0001,
                                   M2CutoffQualityScore  = NULL,
                                   M2SlidingWindowSize   = NULL,
                                   baseNumPerRow         = 100,
                                   heightPerRow          = 200,
                                   signalRatioCutoff     = 0.33,
                                   showTrimmed           = TRUE)

# 3. FILE_TYPE_ERROR: fasta file
SRfa_FILE_TYPE_ERROR <- file.path(inputFilesPath,
                                  "fasta",
                                  "SangerRead",
                                  "Achl_ACHLO006-09_1_F.fast")
sangerReadF_FILE_TYPE_ERROR <- new("SangerRead",
                                   readFeature           = "Forward Read",
                                   inputSource           = "FASTA",
                                   readFileName          = SRfa_FILE_TYPE_ERROR,
                                   geneticCode           = GENETIC_CODE,
                                   TrimmingMethod        = "M1",
                                   M1TrimmingCutoff      = 0.0001,
                                   M2CutoffQualityScore  = NULL,
                                   M2SlidingWindowSize   = NULL,
                                   baseNumPerRow         = 100,
                                   heightPerRow          = 200,
                                   signalRatioCutoff     = 0.33,
                                   showTrimmed           = TRUE)

