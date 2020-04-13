inputFilesPath <- system.file("extdata/", package = "sangeranalyseR")
SRfastaRFN <- file.path(inputFilesPath,
                        "fasta",
                        "SangerRead",
                        "Achl_ACHLO006-09_2_R.fa")
fastaReadRName <- "Achl_ACHLO006-09_2_R"
sangerReadRFa <- new("SangerRead",
                     inputSource        = "FASTA",
                     readFeature        = "Reverse Read",
                     readFileName       = SRfastaRFN,
                     fastaReadName      = fastaReadRName,
                     geneticCode        = GENETIC_CODE)
