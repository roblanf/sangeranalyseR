inputFilesPath <- system.file("extdata/", package = "sangeranalyseR")
SRfastaFFN <- file.path(inputFilesPath,
                        "fasta",
                        "SangerRead",
                        "Achl_ACHLO006-09_1_F.fa")
namesConversionCSV <- file.path(inputFilesPath, "fasta", "SangerRead", "names_conversion_1.csv")
fastaReadFName <- "Achl_ACHLO006-09_1_F"
sangerReadFFa <- new("SangerRead",
                     inputSource        = "FASTA",
                     readFeature        = "Forward Read",
                     readFileName       = SRfastaFFN,
                     fastaReadName      = fastaReadFName,
                     namesConversionCSV = namesConversionCSV,
                     geneticCode        = GENETIC_CODE)
