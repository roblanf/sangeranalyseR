rawDataDir <- system.file("extdata", package = "sangeranalyseR")
inputFilesParentDir <- file.path(rawDataDir, "Allolobophora_chlorotica", "ACHLO")
contigName <- "ACHLO006-09[LCO1490_t1,HCO2198_t1]"
suffixForwardRegExp <- "_[F]_[0-9]*.ab1"
suffixReverseRegExp <- "_[R]_[0-9]*.ab1"
sangerContig <- new("SangerContig",
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
forwardReadList <- sangerContig@forwardReadList
reverseReadList <- sangerContig@reverseReadList
refAminoAcidSeq <- sangerContig@refAminoAcidSeq
minFractionCall <- sangerContig@minFractionCall
maxFractionLost <- sangerContig@maxFractionLost
geneticCode <- sangerContig@geneticCode
acceptStopCodons <- sangerContig@acceptStopCodons
readingFrame <- sangerContig@readingFrame
