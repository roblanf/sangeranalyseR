inputFilesPath <- system.file("extdata/", package = "sangeranalyseR")
SCfastaFN <- file.path(inputFilesPath, "fasta",
                     "SangerContig", "Achl_ACHLO006-09.fa")
contigName <- "Achl_ACHLO006-09"
suffixForwardRegExpFa <- "_[0-9]*_F"
suffixReverseRegExpFa <- "_[0-9]*_R"
namesConversionCSV <- file.path(rawDataDir, "fasta", "SangerContig", "names_conversion_1.csv")
sangerContigFa <- new("SangerContig",
                      inputSource           = "FASTA",
                      fastaFileName         = SCfastaFN,
                      contigName            = contigName,
                      suffixForwardRegExp   = suffixForwardRegExpFa,
                      suffixReverseRegExp   = suffixReverseRegExpFa,
                      namesConversionCSV    = namesConversionCSV,
                      refAminoAcidSeq = "SRQWLFSTNHKDIGTLYFIFGAWAGMVGTSLSILIRAELGHPGALIGDDQIYNVIVTAHAFIMIFFMVMPIMIGGFGNWLVPLMLGAPDMAFPRMNNMSFWLLPPALSLLLVSSMVENGAGTGWTVYPPLSAGIAHGGASVDLAIFSLHLAGISSILGAVNFITTVINMRSTGISLDRMPLFVWSVVITALLLLLSLPVLAGAITMLLTDRNLNTSFFDPAGGGDPILYQHLFWFFGHPEVYILILPGFGMISHIISQESGKKETFGSLGMIYAMLAIGLLGFIVWAHHMFTVGMDVDTRAYFTSATMIIAVPTGIKIFSWLATLHGTQLSYSPAILWALGFVFLFTVGGLTGVVLANSSVDIILHDTYYVVAHFHYVLSMGAVFAIMAGFIHWYPLFTGLTLNNKWLKSHFIIMFIGVNLTFFPQHFLGLAGMPRRYSDYPDAYTTWNIVSTIGSTISLLGILFFFFIIWESLVSQRQVIYPIQLNSSIEWYQNTPPAEHSYSELPLLTN",
                      processorsNum         = 2)


SCfastaFN2 <- file.path(inputFilesPath, "fasta",
                       "SangerContig", "Achl_RBNII384-13.fa")
contigName2 <- "Achl_RBNII384-13"
namesConversionCSV <- file.path(rawDataDir, "fasta", "SangerContig", "names_conversion_2.csv")
sangerContigFa2 <- new("SangerContig",
                      inputSource           = "FASTA",
                      fastaFileName         = SCfastaFN2,
                      namesConversionCSV    = namesConversionCSV,
                      contigName            = contigName2,
                      suffixForwardRegExp   = suffixForwardRegExpFa,
                      suffixReverseRegExp   = suffixReverseRegExpFa,
                      refAminoAcidSeq = "SRQWLFSTNHKDIGTLYFIFGAWAGMVGTSLSILIRAELGHPGALIGDDQIYNVIVTAHAFIMIFFMVMPIMIGGFGNWLVPLMLGAPDMAFPRMNNMSFWLLPPALSLLLVSSMVENGAGTGWTVYPPLSAGIAHGGASVDLAIFSLHLAGISSILGAVNFITTVINMRSTGISLDRMPLFVWSVVITALLLLLSLPVLAGAITMLLTDRNLNTSFFDPAGGGDPILYQHLFWFFGHPEVYILILPGFGMISHIISQESGKKETFGSLGMIYAMLAIGLLGFIVWAHHMFTVGMDVDTRAYFTSATMIIAVPTGIKIFSWLATLHGTQLSYSPAILWALGFVFLFTVGGLTGVVLANSSVDIILHDTYYVVAHFHYVLSMGAVFAIMAGFIHWYPLFTGLTLNNKWLKSHFIIMFIGVNLTFFPQHFLGLAGMPRRYSDYPDAYTTWNIVSTIGSTISLLGILFFFFIIWESLVSQRQVIYPIQLNSSIEWYQNTPPAEHSYSELPLLTN",
                      processorsNum         = 2)

