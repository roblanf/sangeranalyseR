rawDataDir <- system.file("extdata", package = "sangeranalyseR")
inputDir <- file.path(rawDataDir, "Allolobophora_chlorotica", "ACHLO")
suffixForwardRegExp <- "_[F]_[0-9]*.ab1"
suffixReverseRegExp <- "_[R]_[0-9]*.ab1"
sangerAlignment <- new("SangerAlignment",
                       parentDirectory       = inputDir,
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



