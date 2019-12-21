test_that("SangerContig Initial test", {
    expect_type(sangerContig, "S4")
    expect_s4_class(sangerContig, "SangerContig")
    expect_equal(sangerContig@inputSource, "ABIF")
    expect_equal(sangerContig@contigName, "ACHLO006-09[LCO1490_t1,HCO2198_t1]")
    expect_equal(sangerContig@suffixForwardRegExp, "_[0-9]*_F.ab1")
    expect_equal(sangerContig@suffixReverseRegExp, "_[0-9]*_R.ab1")
    expect_equal(sangerContig@trimmingMethodSC, "M1")
    expect_equal(sangerContig@minReadsNum, 2)
    expect_equal(sangerContig@minReadLength, 20)
    expect_equal(sangerContig@refAminoAcidSeq, "SRQWLFSTNHKDIGTLYFIFGAWAGMVGTSLSILIRAELGHPGALIGDDQIYNVIVTAHAFIMIFFMVMPIMIGGFGNWLVPLMLGAPDMAFPRMNNMSFWLLPPALSLLLVSSMVENGAGTGWTVYPPLSAGIAHGGASVDLAIFSLHLAGISSILGAVNFITTVINMRSTGISLDRMPLFVWSVVITALLLLLSLPVLAGAITMLLTDRNLNTSFFDPAGGGDPILYQHLFWFFGHPEVYILILPGFGMISHIISQESGKKETFGSLGMIYAMLAIGLLGFIVWAHHMFTVGMDVDTRAYFTSATMIIAVPTGIKIFSWLATLHGTQLSYSPAILWALGFVFLFTVGGLTGVVLANSSVDIILHDTYYVVAHFHYVLSMGAVFAIMAGFIHWYPLFTGLTLNNKWLKSHFIIMFIGVNLTFFPQHFLGLAGMPRRYSDYPDAYTTWNIVSTIGSTISLLGILFFFFIIWESLVSQRQVIYPIQLNSSIEWYQNTPPAEHSYSELPLLTN")
    expect_equal(sangerContig@minFractionCall, 0.5)
    expect_equal(sangerContig@maxFractionLost, 0.5)
    expect_equal(sangerContig@acceptStopCodons, TRUE)
    expect_equal(sangerContig@readingFrame, 1)
    expect_equal(as.character(sangerContig@contigSeq), "TATATCGACGGCCAGTGGTCAACAAATCATAAAGATATTGGAACTTTATATTTTATTCTGGGCGTCTGAGCAGGAATGGTTGGAGCCGGTATAAGACTTCTAATTCGAATCGAGCTAAGACAACCAGGAGCGTTCCTGGGCAGAGACCAACTATACAATACTATCGTTACTGCACACGCATTTGTAATAATCTTCTTTCTAGTAATGCCTGTATTCATCGGGGGATTCGGAAACTGGCTTTTACCTTTAATACTTGGAGCCCCCGATATAGCATTCCCTCGACTCAACAACATGAGATTCTGACTACTTCCCCCATCACTGATCCTTTTAGTGTCCTCTGCGGCGGTAGAAAAAGGCGCTGGTACGGGGTGAACTGTTTATCCGCCTCTAGCAAGAAATCTTGCCCACGCAGGCCCGTCTGTAGATTTAGCCATCTTTTCCCTTCATTTAGCGGGTGCGTCTTCTATTCTAGGGGCTATTAATTTTATCACCACAGTTATTAATATGCGTTGAAGAGGATTACGTCTTGAACGAATTCCCCTGTTTGTCTGAGCTGTGCTAATTACAGTTGTTCTTCTACTTCTATCTTTACCAGTGCTAGCAGGTGCCATTACCATACTTCTTACCGAC")
})

test_that("SangerContig update quality trimming parameters test 1", {
    newTrimmedSangerContig <- updateQualityParam(sangerContig,
                                                 TrimmingMethod       = "M1",
                                                 M1TrimmingCutoff     = 0.1,
                                                 M2CutoffQualityScore = NULL,
                                                 M2SlidingWindowSize  = NULL)
    sapply(newTrimmedSangerContig@forwardReadList, function(forwardRead) {
        expect_equal(forwardRead@QualityReport@TrimmingMethod, "M1")
        expect_equal(forwardRead@QualityReport@M1TrimmingCutoff, 0.1, tolerance=1e-10)
        expect_equal(forwardRead@QualityReport@M2CutoffQualityScore, NULL)
        expect_equal(forwardRead@QualityReport@M2SlidingWindowSize, NULL)
    })
    sapply(newTrimmedSangerContig@reverseReadList, function(reverseRead) {
        expect_equal(reverseRead@QualityReport@TrimmingMethod, "M1")
        expect_equal(reverseRead@QualityReport@M1TrimmingCutoff, 0.1, tolerance=1e-10)
        expect_equal(reverseRead@QualityReport@M2CutoffQualityScore, NULL)
        expect_equal(reverseRead@QualityReport@M2SlidingWindowSize, NULL)
    })
})
