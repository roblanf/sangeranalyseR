test_that("SangerContig Initial test", {
    expect_type(sangerContig, "S4")
    expect_s4_class(sangerContig, "SangerContig")
    expect_equal(sangerContig@inputSource, "ABIF")
    expect_equal(sangerContig@contigName, "Achl_ACHLO006-09")
    expect_equal(sangerContig@REGEX_SuffixForward, "_[0-9]*_F.ab1")
    expect_equal(sangerContig@REGEX_SuffixReverse, "_[0-9]*_R.ab1")
    expect_equal(sangerContig@minReadsNum, 2)
    expect_equal(sangerContig@minReadLength, 20)
    expect_equal(sangerContig@refAminoAcidSeq, "SRQWLFSTNHKDIGTLYFIFGAWAGMVGTSLSILIRAELGHPGALIGDDQIYNVIVTAHAFIMIFFMVMPIMIGGFGNWLVPLMLGAPDMAFPRMNNMSFWLLPPALSLLLVSSMVENGAGTGWTVYPPLSAGIAHGGASVDLAIFSLHLAGISSILGAVNFITTVINMRSTGISLDRMPLFVWSVVITALLLLLSLPVLAGAITMLLTDRNLNTSFFDPAGGGDPILYQHLFWFFGHPEVYILILPGFGMISHIISQESGKKETFGSLGMIYAMLAIGLLGFIVWAHHMFTVGMDVDTRAYFTSATMIIAVPTGIKIFSWLATLHGTQLSYSPAILWALGFVFLFTVGGLTGVVLANSSVDIILHDTYYVVAHFHYVLSMGAVFAIMAGFIHWYPLFTGLTLNNKWLKSHFIIMFIGVNLTFFPQHFLGLAGMPRRYSDYPDAYTTWNIVSTIGSTISLLGILFFFFIIWESLVSQRQVIYPIQLNSSIEWYQNTPPAEHSYSELPLLTN")
    expect_equal(sangerContig@minFractionCall, 0.5)
    expect_equal(sangerContig@maxFractionLost, 0.5)
    expect_equal(sangerContig@acceptStopCodons, TRUE)
    expect_equal(sangerContig@readingFrame, 1)
    expect_equal(as.character(sangerContig@contigSeq), "TTATATTTTATTCTGGGCGTCTGAGCAGGAATGGTTGGAGCCGGTATAAGACTTCTAATTCGAATCGAGCTAAGACAACCAGGAGCGTTCCTGGGCAGAGACCAACTATACAATACTATCGTTACTGCACACGCATTTGTAATAATCTTCTTTCTAGTAATGCCTGTATTCATCGGGGGATTCGGAAACTGGCTTTTACCTTTAATACTTGGAGCCCCCGATATAGCATTCCCTCGACTCAACAACATGAGATTCTGACTACTTCCCCCATCACTGATCCTTTTAGTGTCCTCTGCGGCGGTAGAAAAAGGCGCTGGTACGGGGTGAACTGTTTATCCGCCTCTAGCAAGAAATCTTGCCCACGCAGGCCCGTCTGTAGATTTAGCCATCTTTTCCCTTCATTTAGCGGGTGCGTCTTCTATTCTAGGGGCTATTAATTTTATCACCACAGTTATTAATATGCGTTGAAGAGGATTACGTCTTGAACGAATTCCCCTGTTTGTCTGAGCTGTGCTAATTACAGTTGTTCTTCTACTTCTATCTTTACCAGTGCTAGCAGGTGCCATTACCATACTTCTTACCGACCGAAACCTCAATACTTCATTCTTTGATCCTGCCGGTGGTGGAGACCCCATCCTC")
})

test_that("SangerContig update quality trimming parameters test 1 - M1", {
    newTrimmedSangerContig <- updateQualityParam(sangerContig,
                                                 TrimmingMethod       = "M1",
                                                 M1TrimmingCutoff     = 0.1,
                                                 M2CutoffQualityScore = NULL,
                                                 M2SlidingWindowSize  = NULL,
                                                 processorsNum        = 1)
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
test_that("SangerContig update quality trimming parameters test 2 - M2", {
    newTrimmedSangerContig <- updateQualityParam(sangerContig,
                                                 TrimmingMethod       = "M2",
                                                 M1TrimmingCutoff     = NULL,
                                                 M2CutoffQualityScore = 30,
                                                 M2SlidingWindowSize  = 15,
                                                 processorsNum        = 1)
    sapply(newTrimmedSangerContig@forwardReadList, function(forwardRead) {
        expect_equal(forwardRead@QualityReport@TrimmingMethod, "M2")
        expect_equal(forwardRead@QualityReport@M1TrimmingCutoff, NULL)
        expect_equal(forwardRead@QualityReport@M2CutoffQualityScore, 30, tolerance=1e-10)
        expect_equal(forwardRead@QualityReport@M2SlidingWindowSize, 15, tolerance=1e-10)
    })
    sapply(newTrimmedSangerContig@reverseReadList, function(reverseRead) {
        expect_equal(reverseRead@QualityReport@TrimmingMethod, "M2")
        expect_equal(reverseRead@QualityReport@M1TrimmingCutoff, NULL)
        expect_equal(reverseRead@QualityReport@M2CutoffQualityScore, 30, tolerance=1e-10)
        expect_equal(reverseRead@QualityReport@M2SlidingWindowSize, 15, tolerance=1e-10)
    })
})

# test_that("SangerContig update quality trimming parameters 3 (M2CutoffQualityScore bigger than threashold)", {
#     expect_error(updateQualityParam(sangerContig,
#                                     TrimmingMethod       = "M2",
#                                     M1TrimmingCutoff     = NULL,
#                                     M2CutoffQualityScore = 61,
#                                     M2SlidingWindowSize  = 41,
#                                     processorsNum        = 1),
#                  "\nYour input M2CutoffQualityScore is: '61' is invalid.'M2CutoffQualityScore' shouldbe between 0 and 60.\n\nYour input M2SlidingWindowSize is: '41' is invalid.'M2SlidingWindowSize' shouldbe between 0 and 40.\n", fixed = TRUE)
# })
# test_that("SangerContig update quality trimming parameters 4 (M2CutoffQualityScore smaller than threashold)", {
#     expect_error(updateQualityParam(sangerContig,
#                                     TrimmingMethod       = "M2",
#                                     M1TrimmingCutoff     = NULL,
#                                     M2CutoffQualityScore = -1,
#                                     M2SlidingWindowSize  = -1,
#                                     processorsNum        = 1),
#                  "\nYour input M2CutoffQualityScore is: '-1' is invalid.'M2CutoffQualityScore' shouldbe between 0 and 60.\n\nYour input M2SlidingWindowSize is: '-1' is invalid.'M2SlidingWindowSize' shouldbe between 0 and 40.\n", fixed = TRUE)
# })

