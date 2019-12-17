rawDataDir <- system.file("extdata", package = "sangeranalyseR")
parentDir <- file.path(rawDataDir, "Allolobophora_chlorotica", "ACHLO")
suffixForwardRegExp <- "_[F]_[0-9]*.ab1"
suffixReverseRegExp <- "_[R]_[0-9]*.ab1"
sangerAlignment <- new("SangerAlignment",
                       inputSource           = "ABIF",
                       parentDirectory       = parentDir,
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
