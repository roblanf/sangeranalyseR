inputFilesPath <- system.file("extdata/", package = "sangeranalyseR")
A_chloroticaFNfa <- file.path(inputFilesPath,
                              "fasta",
                              "SangerRead",
                              "ACHLO006-09[LCO1490_t1,HCO2198_t1]_F_1.fa")
fastaReadName <- "ACHLO006-09[LCO1490_t1,HCO2198_t1]_F_1.fa"
sangerReadFa <- new("SangerRead",
                     inputSource   = "FASTA",
                     readFeature   = "Forward Read",
                     readFileName  = A_chloroticaFNfa,
                     fastaReadName = fastaReadName,
                     geneticCode   = GENETIC_CODE)
