### ============================================================================
### SangerRead FASTA Initial test
### ============================================================================
test_that("SangerRead - FASTA Initial test - forward", {
    expect_type(sangerReadFFa, "S4")
    expect_s4_class(sangerReadFFa, "SangerRead")
    expect_equal(sangerReadFFa@inputSource, "FASTA")
    expect_equal(sangerReadFFa@readFeature, "Forward Read")
    expect_equal(basename(sangerReadFFa@readFileName), "Achl_ACHLO006-09_1_F.fa")
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
test_that("SangerRead - FASTA Initial test - reverse", {
    expect_type(sangerReadRFa, "S4")
    expect_s4_class(sangerReadRFa, "SangerRead")
    expect_equal(sangerReadRFa@inputSource, "FASTA")
    expect_equal(sangerReadRFa@readFeature, "Reverse Read")
    expect_equal(basename(sangerReadRFa@readFileName), "Achl_ACHLO006-09_2_R.fa")
    expect_equal(as.character(sangerReadRFa@primarySeq), "GTATATCGACGGCCAGTGGTCAACAAATCATAAAGATATTGGAACTTTATATTTTATTCTGGGCGTCTGAGCAGGAATGGTTGGAGCCGGTATAAGACTTCTAATTCGAATCGAGCTAAGACAACCAGGAGCGTTCCTGGGCAGAGACCAACTATACAATACTATCGTTACTGCACACGCATTTGTAATAATCTTCTTTCTAGTAATGCCTGTATTCATCGGGGGATTCGGAAACTGGCTTTTACCTTTAATACTTGGAGCCCCCGATATAGCATTCCCTCGACTCAACAACATGAGATTCTGACTACTTCCCCCATCACTGATCCTTTTAGTGTCCTCTGCGGCGGTAGAAAAAGGCGCTGGTACGGGGTGAACTGTTTATCCGCCTCTAGCAAGAAATCTTGCCCACGCAGGCCCGTCTGTAGATTTAGCCATCTTTTCCCTTCATTTAGCGGGTGCGTCTTCTATTCTAGGGGCTATTAATTTTATCACCACAGTTATTAATATGCGTTGAAGAGGATTACGTCTTGAACGAATTCCCCTGTTTGTCTGAGCTGTGCTAATTACAGTTGTTCTTCTACTTCTATCTTTACCAGTGCTAGCAGGTGCCATTACCATACTTCTTACCGACCG")
    expect_equal(as.character(sangerReadRFa@secondarySeq), "")
    expect_equal(as.character(sangerReadRFa@primaryAASeqS1), "VYRRPVVNKS*RYWNFIFYSGRLSRNGWSRYKTSNSNRAKTTRSVPGQRPTIQYYRYCTRICNNLLSSNACIHRGIRKLAFTFNTWSPRYSIPSTQQHEILTTSPITDPFSVLCGGRKRRWYGVNCLSASSKKSCPRRPVCRFSHLFPSFSGCVFYSRGY*FYHHSY*YALKRITS*TNSPVCLSCANYSCSSTSIFTSASRCHYHTSYRP")
    expect_equal(as.character(sangerReadRFa@primaryAASeqS2), "YIDGQWSTNHKDIGTLYFILGV*AGMVGAGIRLLIRIELRQPGAFLGRDQLYNTIVTAHAFVIIFFLVMPVFIGGFGNWLLPLILGAPDIAFPRLNNMRF*LLPPSLILLVSSAAVEKGAGTG*TVYPPLARNLAHAGPSVDLAIFSLHLAGASSILGAINFITTVINMR*RGLRLERIPLFV*AVLITVVLLLLSLPVLAGAITILLTD")
    expect_equal(as.character(sangerReadRFa@primaryAASeqS3), "ISTASGQQIIKILELYILFWASEQEWLEPV*DF*FESS*DNQERSWAETNYTILSLLHTHL**SSF**CLYSSGDSETGFYL*YLEPPI*HSLDSTT*DSDYFPHH*SF*CPLRR*KKALVRGELFIRL*QEILPTQARL*I*PSFPFI*RVRLLF*GLLILSPQLLICVEEDYVLNEFPCLSELC*LQLFFYFYLYQC*QVPLPYFLPT")
    expect_equal(is.null(sangerReadRFa@ChromatogramParam), TRUE)
    expect_equal(is.null(sangerReadRFa@QualityReport), TRUE)
    expect_equal(is.null(sangerReadRFa@abifRawData), TRUE)
    expect_equal(as.character(sangerReadRFa@primarySeqRaw), "")
    expect_equal(as.character(sangerReadRFa@secondarySeqRaw), "")
    expect_equal(sangerReadRFa@primarySeqID, "From fasta file")
    expect_equal(sangerReadRFa@secondarySeqID, "")
    expect_equal(sangerReadRFa@traceMatrix, matrix())
    expect_equal(sangerReadRFa@peakPosMatrix, matrix())
    expect_equal(sangerReadRFa@peakAmpMatrix, matrix())
})


### ============================================================================
### SangerRead FASTA Update Quality Trimming test
### ============================================================================
# test_that("SangerRead update quality trimming parameters 10 (M2CutoffQualityScore smaller than threashold)", {
#     expect_error(new("SangerRead",
#                      inputSource   = "FASTA",
#                      readFeature   = "Forward Read",
#                      readFileName  = SRfastaFFN,
#                      fastaReadName = "Random_read_name",
#                      geneticCode   = GENETIC_CODE),
#                  "The name 'Random_read_name' is not in the 'Achl_ACHLO006-09_1_F.fa' FASTA file", fixed = TRUE)
# })

### ============================================================================
### SangerRead Functions test
### ============================================================================
test_that("SangerRead - FASTA functions test -forward", {
    # expect_message(qualityBasePlot(sangerReadFFa),
    #              paste("SangerRead with 'FASTA' inputSource",
    #                    "cannot plot quality plots"))
    # expect_message(updateQualityParam(sangerReadFFa),
    #              paste("SangerRead with 'FASTA' inputSource",
    #              "cannot update quality parameters"))
    # expect_message(MakeBaseCalls(sangerReadFFa),
    #              paste("SangerRead with 'FASTA'",
    #                    "inputSource cannot do base calling"))

    outputFasta <- writeFastaSR(sangerReadFFa)
    expect_true(file.exists(outputFasta))
    con = file(outputFasta, "r")
    line = readLines(con, n = 1)
    expect_equal(line, ">Achl_ACHLO006-09_1_F.fa")
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
test_that("SangerRead - FASTA functions test -reverse", {
    # expect_message(qualityBasePlot(sangerReadRFa),
    #                paste("SangerRead with 'FASTA' inputSource",
    #                      "cannot plot quality plots"))
    # expect_message(updateQualityParam(sangerReadRFa),
    #                paste("SangerRead with 'FASTA' inputSource",
    #                      "cannot update quality parameters"))
    # expect_message(MakeBaseCalls(sangerReadRFa),
    #                paste("SangerRead with 'FASTA'",
    #                      "inputSource cannot do base calling"))

    outputFasta <- writeFastaSR(sangerReadRFa)
    expect_true(file.exists(outputFasta))
    con = file(outputFasta, "r")
    line = readLines(con, n = 1)
    expect_equal(line, ">Achl_ACHLO006-09_2_R.fa")
    line = readLines(con, n = 1)
    expect_equal(line, "CGGTCGGTAAGAAGTATGGTAATGGCACCTGCTAGCACTGGTAAAGATAGAAGTAGAAGAACAACTGTAATTAGCACAGC")
    line = readLines(con, n = 1)
    expect_equal(line, "TCAGACAAACAGGGGAATTCGTTCAAGACGTAATCCTCTTCAACGCATATTAATAACTGTGGTGATAAAATTAATAGCCC")
    line = readLines(con, n = 1)
    expect_equal(line, "CTAGAATAGAAGACGCACCCGCTAAATGAAGGGAAAAGATGGCTAAATCTACAGACGGGCCTGCGTGGGCAAGATTTCTT")
    line = readLines(con, n = 1)
    expect_equal(line, "GCTAGAGGCGGATAAACAGTTCACCCCGTACCAGCGCCTTTTTCTACCGCCGCAGAGGACACTAAAAGGATCAGTGATGG")
    line = readLines(con, n = 1)
    expect_equal(line, "GGGAAGTAGTCAGAATCTCATGTTGTTGAGTCGAGGGAATGCTATATCGGGGGCTCCAAGTATTAAAGGTAAAAGCCAGT")
    line = readLines(con, n = 1)
    expect_equal(line, "TTCCGAATCCCCCGATGAATACAGGCATTACTAGAAAGAAGATTATTACAAATGCGTGTGCAGTAACGATAGTATTGTAT")
    line = readLines(con, n = 1)
    expect_equal(line, "AGTTGGTCTCTGCCCAGGAACGCTCCTGGTTGTCTTAGCTCGATTCGAATTAGAAGTCTTATACCGGCTCCAACCATTCC")
    line = readLines(con, n = 1)
    expect_equal(line, "TGCTCAGACGCCCAGAATAAAATATAAAGTTCCAATATCTTTATGATTTGTTGACCACTGGCCGTCGATATAC")
    close(con)
    unlink(outputFasta)
})
