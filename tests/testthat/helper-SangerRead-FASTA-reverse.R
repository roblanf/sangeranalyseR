inputFilesPath <- system.file("extdata/", package = "sangeranalyseR")
SRfastaRFN <- file.path(inputFilesPath,
                        "fasta",
                        "SangerRead",
                        "ACHLO006-09[LCO1490_t1,HCO2198_t1]_2_R.fa")
fastaReadRName <- "ACHLO006-09[LCO1490_t1,HCO2198_t1]_2_R"
sangerReadRFa <- new("SangerRead",
                     inputSource   = "FASTA",
                     readFeature   = "Reverse Read",
                     readFileName  = SRfastaRFN,
                     fastaReadName = fastaReadRName,
                     geneticCode   = GENETIC_CODE)
