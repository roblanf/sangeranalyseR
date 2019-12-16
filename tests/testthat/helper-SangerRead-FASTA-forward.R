inputFilesPath <- system.file("extdata/", package = "sangeranalyseR")
SRfastaFFN <- file.path(inputFilesPath,
                        "fasta",
                        "SangerRead",
                        "ACHLO006-09[LCO1490_t1,HCO2198_t1]_F_1.fa")
fastaReadName <- "ACHLO006-09[LCO1490_t1,HCO2198_t1]_F_1.fa"
sangerReadFa <- new("SangerRead",
                    inputSource   = "FASTA",
                    readFeature   = "Forward Read",
                    readFileName  = SRfastaFFN,
                    fastaReadName = fastaReadName,
                    geneticCode   = GENETIC_CODE)
