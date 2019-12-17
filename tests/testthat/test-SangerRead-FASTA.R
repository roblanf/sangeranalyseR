### ============================================================================
### SangerRead FASTA Initial test
### ============================================================================
test_that("SangerRead - FASTA Initial test", {
    expect_type(sangerReadFFa, "S4")
    expect_s4_class(sangerReadFFa, "SangerRead")
    expect_equal(sangerReadFFa@inputSource, "FASTA")
    expect_equal(sangerReadFFa@readFeature, "Forward Read")
    expect_equal(basename(sangerReadFFa@readFileName), "ACHLO006-09[LCO1490_t1,HCO2198_t1]_F_1.fa")
    expect_equal(as.character(sangerReadFFa@primarySeq), "CTGGGCGTCTGAGCAGGAATGGTTGGAGCCGGTATAAGACTTCTAATTCGAATCGAGCTAAGACAACCAGGAGCGTTCCTGGGCAGAGACCAACTATACAATACTATCGTTACTGCACACGCATTTGTAATAATCTTCTTTCTAGTAATGCCTGTATTCATCGGGGGATTCGGAAACTGGCTTTTACCTTTAATACTTGGAGCCCCCGATATAGCATTCCCTCGACTCAACAACATGAGATTCTGACTACTTCCCCCATCACTGATCCTTTTAGTGTCCTCTGCGGCGGTAGAAAAAGGCGCTGGTACGGGGTGAACTGTTTATCCGCCTCTAGCAAGAAATCTTGCCCACGCAGGCCCGTCTGTAGATTTAGCCATCTTTTCCCTTCATTTAGCGGGTGCGTCTTCTATTCTAGGGGCTATTAATTTTATCACCACAGTTATTAATATGCGTTGAAGAGG")
    expect_equal(as.character(sangerReadFFa@secondarySeq), "")
    expect_equal(as.character(sangerReadFFa@primaryAASeqS1), "LGV*AGMVGAGIRLLIRIELRQPGAFLGRDQLYNTIVTAHAFVIIFFLVMPVFIGGFGNWLLPLILGAPDIAFPRLNNMRF*LLPPSLILLVSSAAVEKGAGTG*TVYPPLARNLAHAGPSVDLAIFSLHLAGASSILGAINFITTVINMR*R")
    expect_equal(as.character(sangerReadFFa@primaryAASeqS2), "WASEQEWLEPV*DF*FESS*DNQERSWAETNYTILSLLHTHL**SSF**CLYSSGDSETGFYL*YLEPPI*HSLDSTT*DSDYFPHH*SF*CPLRR*KKALVRGELFIRL*QEILPTQARL*I*PSFPFI*RVRLLF*GLLILSPQLLICVEE")
    expect_equal(as.character(sangerReadFFa@primaryAASeqS3), "GRLSRNGWSRYKTSNSNRAKTTRSVPGQRPTIQYYRYCTRICNNLLSSNACIHRGIRKLAFTFNTWSPRYSIPSTQQHEILTTSPITDPFSVLCGGRKRRWYGVNCLSASSKKSCPRRPVCRFSHLFPSFSGCVFYSRGY*FYHHSY*YALKR")
    expect_equal(is.null(sangerReadFFa@ChromatogramParam), TRUE)
    expect_equal(is.null(sangerReadFFa@QualityReport), TRUE)
    expect_equal(is.null(sangerReadFFa@abifRawData), TRUE)
    expect_equal(as.character(sangerReadFFa@primarySeqRaw), "")
    expect_equal(as.character(sangerReadFFa@secondarySeqRaw), "")
    expect_equal(sangerReadFFa@primarySeqID, "From fasta file")
    expect_equal(sangerReadFFa@secondarySeqID, "")
    expect_equal(sangerReadFFa@traceMatrix, matrix())
    expect_equal(sangerReadFFa@peakPosMatrix, matrix())
    expect_equal(sangerReadFFa@peakAmpMatrix, matrix())
})

### ============================================================================
### SangerRead FASTA Update Quality Trimming test
### ============================================================================
test_that("SangerRead update quality trimming parameters 10 (M2CutoffQualityScore smaller than threashold)", {
    expect_error(new("SangerRead",
                     inputSource   = "FASTA",
                     readFeature   = "Forward Read",
                     readFileName  = SRfastaFFN,
                     fastaReadName = "Random_read_name",
                     geneticCode   = GENETIC_CODE),
                 "The name 'Random_read_name' is not in the 'ACHLO006-09[LCO1490_t1,HCO2198_t1]_F_1.fa' FASTA file", fixed = TRUE)
})

### ============================================================================
### SangerRead Functions test
### ============================================================================
test_that("SangerRead - FASTA functions test", {
    expect_message(qualityBasePlot(sangerReadFFa),
                 paste("SangerRead with 'FASTA' inputSource",
                       "cannot plot quality plots"))
    expect_message(updateQualityParam(sangerReadFFa),
                 paste("SangerRead with 'FASTA' inputSource",
                 "cannot update quality parameters"))
    expect_message(MakeBaseCalls(sangerReadFFa),
                 paste("SangerRead with 'FASTA'",
                       "inputSource cannot do base calling"))

    outputFasta <- writeFastaSR(sangerReadFFa)
    expect_true(file.exists(outputFasta))
    con = file(outputFasta, "r")
    line = readLines(con, n = 1)
    expect_equal(line, ">ACHLO006-09[LCO1490_t1,HCO2198_t1]_F_1.fa")
    line = readLines(con, n = 1)
    expect_equal(line, "CTGGGCGTCTGAGCAGGAATGGTTGGAGCCGGTATAAGACTTCTAATTCGAATCGAGCTAAGACAACCAGGAGCGTTCCT")
    line = readLines(con, n = 1)
    expect_equal(line, "GGGCAGAGACCAACTATACAATACTATCGTTACTGCACACGCATTTGTAATAATCTTCTTTCTAGTAATGCCTGTATTCA")
    line = readLines(con, n = 1)
    expect_equal(line, "TCGGGGGATTCGGAAACTGGCTTTTACCTTTAATACTTGGAGCCCCCGATATAGCATTCCCTCGACTCAACAACATGAGA")
    line = readLines(con, n = 1)
    expect_equal(line, "TTCTGACTACTTCCCCCATCACTGATCCTTTTAGTGTCCTCTGCGGCGGTAGAAAAAGGCGCTGGTACGGGGTGAACTGT")
    line = readLines(con, n = 1)
    expect_equal(line, "TTATCCGCCTCTAGCAAGAAATCTTGCCCACGCAGGCCCGTCTGTAGATTTAGCCATCTTTTCCCTTCATTTAGCGGGTG")
    line = readLines(con, n = 1)
    expect_equal(line, "CGTCTTCTATTCTAGGGGCTATTAATTTTATCACCACAGTTATTAATATGCGTTGAAGAGG")
    close(con)
    unlink(outputFasta)
})
