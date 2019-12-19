inputFilesPath <- system.file("extdata/", package = "sangeranalyseR")
inputFilesParentDir <- file.path(inputFilesPath, "Allolobophora_chlorotica", "ACHLO")
contigName <- "ACHLO006-09[LCO1490_t1,HCO2198_t1]"
suffixForwardRegExp <- "_F.ab1"
suffixReverseRegExp <- "_R.ab1"
sangerContig <- new("SangerContig",
                    inputSource           = "ABIF",
                    parentDirectory       = inputFilesParentDir,
                    contigName            = contigName,
                    suffixForwardRegExp   = suffixForwardRegExp,
                    suffixReverseRegExp   = suffixReverseRegExp,
                    refAminoAcidSeq = "SRQWLFSTNHKDIGTLYFIFGAWAGMVGTSLSILIRAELGHPGALIGDDQIYNVIVTAHAFIMIFFMVMPIMIGGFGNWLVPLMLGAPDMAFPRMNNMSFWLLPPALSLLLVSSMVENGAGTGWTVYPPLSAGIAHGGASVDLAIFSLHLAGISSILGAVNFITTVINMRSTGISLDRMPLFVWSVVITALLLLLSLPVLAGAITMLLTDRNLNTSFFDPAGGGDPILYQHLFWFFGHPEVYILILPGFGMISHIISQESGKKETFGSLGMIYAMLAIGLLGFIVWAHHMFTVGMDVDTRAYFTSATMIIAVPTGIKIFSWLATLHGTQLSYSPAILWALGFVFLFTVGGLTGVVLANSSVDIILHDTYYVVAHFHYVLSMGAVFAIMAGFIHWYPLFTGLTLNNKWLKSHFIIMFIGVNLTFFPQHFLGLAGMPRRYSDYPDAYTTWNIVSTIGSTISLLGILFFFFIIWESLVSQRQVIYPIQLNSSIEWYQNTPPAEHSYSELPLLTN",
                    TrimmingMethod        = "M1",
                    M1TrimmingCutoff      = 0.0001,
                    M2CutoffQualityScore  = NULL,
                    M2SlidingWindowSize   = NULL,
                    baseNumPerRow         = 100,
                    heightPerRow          = 200,
                    signalRatioCutoff     = 0.33,
                    showTrimmed           = TRUE)

inputFilesParentDir2 <- file.path(inputFilesPath, "Allolobophora_chlorotica", "RBNII")
contigName2 <- "RBNII395-13[C_LepFolF,C_LepFolR]"

sangerContig2 <- new("SangerContig",
                    inputSource           = "ABIF",
                    parentDirectory       = inputFilesParentDir2,
                    contigName            = contigName2,
                    suffixForwardRegExp   = suffixForwardRegExp,
                    suffixReverseRegExp   = suffixReverseRegExp,
                    refAminoAcidSeq = "SRQWLFSTNHKDIGTLYFIFGAWAGMVGTSLSILIRAELGHPGALIGDDQIYNVIVTAHAFIMIFFMVMPIMIGGFGNWLVPLMLGAPDMAFPRMNNMSFWLLPPALSLLLVSSMVENGAGTGWTVYPPLSAGIAHGGASVDLAIFSLHLAGISSILGAVNFITTVINMRSTGISLDRMPLFVWSVVITALLLLLSLPVLAGAITMLLTDRNLNTSFFDPAGGGDPILYQHLFWFFGHPEVYILILPGFGMISHIISQESGKKETFGSLGMIYAMLAIGLLGFIVWAHHMFTVGMDVDTRAYFTSATMIIAVPTGIKIFSWLATLHGTQLSYSPAILWALGFVFLFTVGGLTGVVLANSSVDIILHDTYYVVAHFHYVLSMGAVFAIMAGFIHWYPLFTGLTLNNKWLKSHFIIMFIGVNLTFFPQHFLGLAGMPRRYSDYPDAYTTWNIVSTIGSTISLLGILFFFFIIWESLVSQRQVIYPIQLNSSIEWYQNTPPAEHSYSELPLLTN",
                    TrimmingMethod        = "M1",
                    M1TrimmingCutoff      = 0.0001,
                    M2CutoffQualityScore  = NULL,
                    M2SlidingWindowSize   = NULL,
                    baseNumPerRow         = 100,
                    heightPerRow          = 200,
                    signalRatioCutoff     = 0.33,
                    showTrimmed           = TRUE)


