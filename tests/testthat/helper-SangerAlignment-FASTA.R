rawDataDir <- system.file("extdata", package = "sangeranalyseR")
fastaFN <- file.path(rawDataDir, "fasta",
                     "SangerAlignment", "Sanger_all_reads.fa")
suffixForwardRegExpFa <- "_[0-9]*_F"
suffixReverseRegExpFa <- "_[0-9]*_R"
namesConversionCSV <- file.path(rawDataDir, "fasta", "SangerAlignment", "names_conversion.csv")
sangerAlignmentFa <- new("SangerAlignment",
                         inputSource           = "FASTA",
                         processMethod         = "REGEX",
                         FASTA_File            = fastaFN,
                         CSV_NamesConversion    = namesConversionCSV,
                         REGEX_SuffixForward   = suffixForwardRegExpFa,
                         REGEX_SuffixReverse   = suffixReverseRegExpFa,
                         refAminoAcidSeq = "SRQWLFSTNHKDIGTLYFIFGAWAGMVGTSLSILIRAELGHPGALIGDDQIYNVIVTAHAFIMIFFMVMPIMIGGFGNWLVPLMLGAPDMAFPRMNNMSFWLLPPALSLLLVSSMVENGAGTGWTVYPPLSAGIAHGGASVDLAIFSLHLAGISSILGAVNFITTVINMRSTGISLDRMPLFVWSVVITALLLLLSLPVLAGAITMLLTDRNLNTSFFDPAGGGDPILYQHLFWFFGHPEVYILILPGFGMISHIISQESGKKETFGSLGMIYAMLAIGLLGFIVWAHHMFTVGMDVDTRAYFTSATMIIAVPTGIKIFSWLATLHGTQLSYSPAILWALGFVFLFTVGGLTGVVLANSSVDIILHDTYYVVAHFHYVLSMGAVFAIMAGFIHWYPLFTGLTLNNKWLKSHFIIMFIGVNLTFFPQHFLGLAGMPRRYSDYPDAYTTWNIVSTIGSTISLLGILFFFFIIWESLVSQRQVIYPIQLNSSIEWYQNTPPAEHSYSELPLLTN",
                         processorsNum         = 2)
