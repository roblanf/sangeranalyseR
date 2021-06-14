inputFilesPath <- system.file("extdata/", package = "sangeranalyseR")
SCfastaFN <- file.path(inputFilesPath, "fasta",
                     "SangerContig", "Achl_ACHLO006-09.fa")
contigName <- "Achl_ACHLO006-09"
suffixForwardRegExpFa <- "_[0-9]*_F"
suffixReverseRegExpFa <- "_[0-9]*_R"
namesConversionCSV <- file.path(inputFilesPath, "fasta", "SangerContig", "names_conversion_1.csv")
sangerContigFa <- new("SangerContig",
                      inputSource           = "FASTA",
                      processMethod         = "REGEX",
                      FASTA_File            = SCfastaFN,
                      contigName            = contigName,
                      REGEX_SuffixForward   = suffixForwardRegExpFa,
                      REGEX_SuffixReverse   = suffixReverseRegExpFa,
                      CSV_NamesConversion    = namesConversionCSV,
                      refAminoAcidSeq = "SRQWLFSTNHKDIGTLYFIFGAWAGMVGTSLSILIRAELGHPGALIGDDQIYNVIVTAHAFIMIFFMVMPIMIGGFGNWLVPLMLGAPDMAFPRMNNMSFWLLPPALSLLLVSSMVENGAGTGWTVYPPLSAGIAHGGASVDLAIFSLHLAGISSILGAVNFITTVINMRSTGISLDRMPLFVWSVVITALLLLLSLPVLAGAITMLLTDRNLNTSFFDPAGGGDPILYQHLFWFFGHPEVYILILPGFGMISHIISQESGKKETFGSLGMIYAMLAIGLLGFIVWAHHMFTVGMDVDTRAYFTSATMIIAVPTGIKIFSWLATLHGTQLSYSPAILWALGFVFLFTVGGLTGVVLANSSVDIILHDTYYVVAHFHYVLSMGAVFAIMAGFIHWYPLFTGLTLNNKWLKSHFIIMFIGVNLTFFPQHFLGLAGMPRRYSDYPDAYTTWNIVSTIGSTISLLGILFFFFIIWESLVSQRQVIYPIQLNSSIEWYQNTPPAEHSYSELPLLTN",
                      processorsNum         = 2)


SCfastaFN2 <- file.path(inputFilesPath, "fasta",
                       "SangerContig", "Achl_RBNII384-13.fa")
contigName2 <- "Achl_RBNII384-13"
namesConversionCSV <- file.path(inputFilesPath, "fasta", "SangerContig", "names_conversion_2.csv")
sangerContigFa2 <- new("SangerContig",
                      inputSource           = "FASTA",
                      processMethod         = "CSV",
                      FASTA_File            = SCfastaFN2,
                      CSV_NamesConversion   = namesConversionCSV,
                      contigName            = contigName2,
                      REGEX_SuffixForward   = suffixForwardRegExpFa,
                      REGEX_SuffixReverse   = suffixReverseRegExpFa,
                      refAminoAcidSeq = "SRQWLFSTNHKDIGTLYFIFGAWAGMVGTSLSILIRAELGHPGALIGDDQIYNVIVTAHAFIMIFFMVMPIMIGGFGNWLVPLMLGAPDMAFPRMNNMSFWLLPPALSLLLVSSMVENGAGTGWTVYPPLSAGIAHGGASVDLAIFSLHLAGISSILGAVNFITTVINMRSTGISLDRMPLFVWSVVITALLLLLSLPVLAGAITMLLTDRNLNTSFFDPAGGGDPILYQHLFWFFGHPEVYILILPGFGMISHIISQESGKKETFGSLGMIYAMLAIGLLGFIVWAHHMFTVGMDVDTRAYFTSATMIIAVPTGIKIFSWLATLHGTQLSYSPAILWALGFVFLFTVGGLTGVVLANSSVDIILHDTYYVVAHFHYVLSMGAVFAIMAGFIHWYPLFTGLTLNNKWLKSHFIIMFIGVNLTFFPQHFLGLAGMPRRYSDYPDAYTTWNIVSTIGSTISLLGILFFFFIIWESLVSQRQVIYPIQLNSSIEWYQNTPPAEHSYSELPLLTN",
                      processorsNum         = 2)

