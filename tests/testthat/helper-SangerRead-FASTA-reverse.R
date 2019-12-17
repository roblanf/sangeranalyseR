inputFilesPath <- system.file("extdata/", package = "sangeranalyseR")
SRfastaRFN <- file.path(inputFilesPath,
                        "fasta",
                        "SangerRead",
                        "ACHLO006-09[LCO1490_t1,HCO2198_t1]_R_2.fa")
fastaReadRName <- "ACHLO006-09[LCO1490_t1,HCO2198_t1]_R_2"
sangerReadRFa <- new("SangerRead",
                     inputSource   = "FASTA",
                     readFeature   = "Reverse Read",
                     readFileName  = SRfastaRFN,
                     fastaReadName = fastaReadRName,
                     geneticCode   = GENETIC_CODE)
