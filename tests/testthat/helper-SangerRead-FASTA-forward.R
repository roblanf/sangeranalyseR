inputFilesPath <- system.file("extdata/", package = "sangeranalyseR")
SRfastaFFN <- file.path(inputFilesPath,
                        "fasta",
                        "SangerRead",
                        "ACHLO006-09[LCO1490_t1,HCO2198_t1]_F.fa")
fastaReadFName <- "ACHLO006-09[LCO1490_t1,HCO2198_t1]_F"
sangerReadFFa <- new("SangerRead",
                     inputSource   = "FASTA",
                     readFeature   = "Forward Read",
                     readFileName  = SRfastaFFN,
                     fastaReadName = fastaReadFName,
                     geneticCode   = GENETIC_CODE)
