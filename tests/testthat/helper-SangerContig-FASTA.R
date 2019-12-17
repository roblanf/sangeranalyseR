inputFilesPath <- system.file("extdata/", package = "sangeranalyseR")
SCfastaFN <- file.path(inputFilesPath, "fasta",
                     "SangerContig", "ACHLO006-09[LCO1490_t1,HCO2198_t1].fa")
contigName <- "ACHLO006-09[LCO1490_t1,HCO2198_t1]"
suffixForwardRegExpFa <- "_[F]_[0-9]*"
suffixReverseRegExpFa <- "_[R]_[0-9]*"
sangerContigFa <- new("SangerContig",
                      inputSource           = "FASTA",
                      fastaFileName         = SCfastaFN,
                      contigName            = contigName,
                      suffixForwardRegExp   = suffixForwardRegExpFa,
                      suffixReverseRegExp   = suffixReverseRegExpFa,
                      refAminoAcidSeq = "SRQWLFSTNHKDIGTLYFIFGAWAGMVGTSLSILIRAELGHPGALIGDDQIYNVIVTAHAFIMIFFMVMPIMIGGFGNWLVPLMLGAPDMAFPRMNNMSFWLLPPALSLLLVSSMVENGAGTGWTVYPPLSAGIAHGGASVDLAIFSLHLAGISSILGAVNFITTVINMRSTGISLDRMPLFVWSVVITALLLLLSLPVLAGAITMLLTDRNLNTSFFDPAGGGDPILYQHLFWFFGHPEVYILILPGFGMISHIISQESGKKETFGSLGMIYAMLAIGLLGFIVWAHHMFTVGMDVDTRAYFTSATMIIAVPTGIKIFSWLATLHGTQLSYSPAILWALGFVFLFTVGGLTGVVLANSSVDIILHDTYYVVAHFHYVLSMGAVFAIMAGFIHWYPLFTGLTLNNKWLKSHFIIMFIGVNLTFFPQHFLGLAGMPRRYSDYPDAYTTWNIVSTIGSTISLLGILFFFFIIWESLVSQRQVIYPIQLNSSIEWYQNTPPAEHSYSELPLLTN"
                      )

