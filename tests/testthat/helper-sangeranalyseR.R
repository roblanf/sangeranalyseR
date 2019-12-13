# inputFilesPath <- system.file("extdata/", package = "sangeranalyseR")
# A_chloroticaFdReadFN <- file.path(inputFilesPath,
#                                   "Allolobophora_chlorotica",
#                                   "RBNII384-13[C_LepFolF,C_LepFolR]_R_2.ab1")
# A_chloroticaRead <- new("SangerRead",
#                         readFeature           = "Forward Read",
#                         readFileName          = A_chloroticaFdReadFN)

rawDataDir <- system.file("extdata", package = "sangeranalyseR")
inputFilesParentDir <- file.path(rawDataDir, "Allolobophora_chlorotica")
contigName <- "ACHLO006-09[LCO1490_t1,HCO2198_t1]"
suffixForwardRegExp <- "_[F]_[0-9]*.ab1"
suffixReverseRegExp <- "_[R]_[0-9]*.ab1"
A_chloroticContig <- new("SangerContig",
                          parentDirectory       = inputFilesParentDir,
                          contigName            = contigName,
                          suffixForwardRegExp   = suffixForwardRegExp,
                          suffixReverseRegExp   = suffixReverseRegExp,
                          refAminoAcidSeq = "SRQWLFSTNHKDIGTLYFIFGAWAGMVGTSLSILIRAELGHPGALIGDDQIYNVIVTAHAFIMIFFMVMPIMIGGFGNWLVPLMLGAPDMAFPRMNNMSFWLLPPALSLLLVSSMVENGAGTGWTVYPPLSAGIAHGGASVDLAIFSLHLAGISSILGAVNFITTVINMRSTGISLDRMPLFVWSVVITALLLLLSLPVLAGAITMLLTDRNLNTSFFDPAGGGDPILYQHLFWFFGHPEVYILILPGFGMISHIISQESGKKETFGSLGMIYAMLAIGLLGFIVWAHHMFTVGMDVDTRAYFTSATMIIAVPTGIKIFSWLATLHGTQLSYSPAILWALGFVFLFTVGGLTGVVLANSSVDIILHDTYYVVAHFHYVLSMGAVFAIMAGFIHWYPLFTGLTLNNKWLKSHFIIMFIGVNLTFFPQHFLGLAGMPRRYSDYPDAYTTWNIVSTIGSTISLLGILFFFFIIWESLVSQRQVIYPIQLNSSIEWYQNTPPAEHSYSELPLLTN",
                          TrimmingMethod        = "M2",
                          M1TrimmingCutoff      = NULL,
                          M2CutoffQualityScore  = 40,
                          M2SlidingWindowSize   = 10,
                          baseNumPerRow         = 100,
                          heightPerRow          = 200,
                          signalRatioCutoff     = 0.33,
                          showTrimmed           = TRUE)
