test_that("SangerRead update quality trimming parameters 1 (smaller M1)", {
    trimmedSangerRead <- updateQualityParam(sangerRead,
                                            TrimmingMethod       = "M1",
                                            M1TrimmingCutoff     = 0.000001,
                                            M2CutoffQualityScore = NULL,
                                            M2SlidingWindowSize  = NULL)
    expect_type(trimmedSangerRead, "S4")
    expect_s4_class(trimmedSangerRead, "SangerRead")
    expect_equal(trimmedSangerRead@inputSource, "ABIF")
    expect_equal(trimmedSangerRead@readFeature, "Forward Read")
    expect_equal(basename(trimmedSangerRead@readFileName), "ACHLO006-09[LCO1490_t1,HCO2198_t1]_F_1.ab1")
    expect_equal(as.character(trimmedSangerRead@primarySeq), "CACTTTATATTTTATTCTGGGCGTCTGAGCAGGAATGGTTGGAGCCGGTATAAGACTTCTAATTCGAATCGAGCTAAGACAACCAGGAGCGTTCCTGGGCAGAGACCAACTATACAATACTATCGTTACTGCACACGCATTTGTAATAATCTTCTTTCTAGTAATGCCTGTATTCATCGGGGGATTCGGAAACTGGCTTTTACCTTTAATACTTGGAGCCCCCGATATAGCATTCCCTCGACTCAACAACATGAGATTCTGACTACTTCCCCCATCACTGATCCTTTTAGTGTCCTCTGCGGCGGTAGAAAAAGGCGCTGGTACGGGGTGAACTGTTTATCCGCCTCTAGCAAGAAATCTTGCCCACGCAGGCCCGTCTGTAGATTTAGCCATCTTTTCCCTTCATTTAGCGGGTGCGTCTTCTATTCTAGGGGCTATTAATTTTATCACCACAGTTATTAATATGCGTTGAAGAGGATTACGTCTTGAACGAATTCCCCTGTTTGTCTGAGCTGTGCTAATTACAGTTGTTCTTCTACTTCTATCTTTACCAGTGCTAGCAGGTGCCATTACCATACTTCTTACCGACCGAAACCTCAATACTTCATTCTTTGATCCTGCCGGTGGTGGAGACCCCATCCTCTACTAGCACTTATTCTGATTTTTAGATCACCCTGATGTTGAGTCATACTGAATTCCTGA")
    expect_equal(as.character(trimmedSangerRead@secondarySeq), "TACTTTATATTTTATTCTGGGCGTCTGAGCAGGAATGGTTGGAGCCGGTATAAGACTTCTAATTCGAATCGAGCTAAGACAACCAGGAGCGTTCCTGGGCAGAGACCAACTATACAATACTATCGTTACTGCACACGCATTTGTAATAATCTTCTTTCTAGTAATGCCTGTATTCATCGGGGGATTCGGAAACTGGCTTTTACCTTTAATACTTGGAGCCCCCGATATAGCATTCCCTCGACTCAACAACATGAGATTCTGACTACTTCCCCCATCACTGATCCTTTTAGTGTCCTCTGCGGCGGTAGAAAAAGGCGCTGGTACGGGGTGAACTGTTTATCCGCCTCTAGCAAGAAATCTTGCCCACGCAGGCCCGTCTGTAGATTTAGCCATCTTTTCCCTTCATTTAGCGGGTGCGTCTTCTATTCTAGGGGCTATTAATTTTATCACCACAGTTATTAATATGCGTTGAAGAGGATTACGTCTTGAACGAATTCCCCTGTTTGTCTGAGCTGTGCTAATTACAGTTGTTCTTCTACTTCTATCTTTACCAGTGCTAGCAGGTGCCATTACCATACTTCTTACCGACCGAAACCTCAATACTTCATTCTTTGATCCTGCCGGTGGTGGAGACCCCATCCTCTACTAGCACTTATTCTGATTTTTCGATCACCCTGATGTTGAGTCATAGTGAATTCCTGA")
    expect_equal(as.character(trimmedSangerRead@primaryAASeqS1), "HFIFYSGRLSRNGWSRYKTSNSNRAKTTRSVPGQRPTIQYYRYCTRICNNLLSSNACIHRGIRKLAFTFNTWSPRYSIPSTQQHEILTTSPITDPFSVLCGGRKRRWYGVNCLSASSKKSCPRRPVCRFSHLFPSFSGCVFYSRGY*FYHHSY*YALKRITS*TNSPVCLSCANYSCSSTSIFTSASRCHYHTSYRPKPQYFIL*SCRWWRPHPLLALILIFRSP*C*VILNS*")
    expect_equal(as.character(trimmedSangerRead@primaryAASeqS2), "TLYFILGV*AGMVGAGIRLLIRIELRQPGAFLGRDQLYNTIVTAHAFVIIFFLVMPVFIGGFGNWLLPLILGAPDIAFPRLNNMRF*LLPPSLILLVSSAAVEKGAGTG*TVYPPLARNLAHAGPSVDLAIFSLHLAGASSILGAINFITTVINMR*RGLRLERIPLFV*AVLITVVLLLLSLPVLAGAITILLTDRNLNTSFFDPAGGGDPILY*HLF*FLDHPDVESY*IP")
    expect_equal(as.character(trimmedSangerRead@primaryAASeqS3), "LYILFWASEQEWLEPV*DF*FESS*DNQERSWAETNYTILSLLHTHL**SSF**CLYSSGDSETGFYL*YLEPPI*HSLDSTT*DSDYFPHH*SF*CPLRR*KKALVRGELFIRL*QEILPTQARL*I*PSFPFI*RVRLLF*GLLILSPQLLICVEEDYVLNEFPCLSELC*LQLFFYFYLYQC*QVPLPYFLPTETSILHSLILPVVETPSSTSTYSDF*ITLMLSHTEFL")
    expect_equal(trimmedSangerRead@ChromatogramParam@baseNumPerRow, 100)
    expect_equal(trimmedSangerRead@ChromatogramParam@heightPerRow, 200)
    expect_equal(trimmedSangerRead@ChromatogramParam@signalRatioCutoff, 0.33)
    expect_equal(trimmedSangerRead@ChromatogramParam@showTrimmed, TRUE)

    expect_equal(trimmedSangerRead@QualityReport@TrimmingMethod, "M1")
    expect_equal(trimmedSangerRead@QualityReport@M1TrimmingCutoff, 0.000001, tolerance=1e-10)
    expect_equal(trimmedSangerRead@QualityReport@M2CutoffQualityScore, NULL)
    expect_equal(trimmedSangerRead@QualityReport@M2SlidingWindowSize, NULL)
    expect_equal(length(trimmedSangerRead@QualityReport@qualityPhredScoresRaw), 703)
    expect_equal(length(trimmedSangerRead@QualityReport@qualityPhredScores), 702)
    expect_equal(length(trimmedSangerRead@QualityReport@qualityBaseScores), 702)

    expect_equal(trimmedSangerRead@QualityReport@rawSeqLength, 702)
    expect_equal(trimmedSangerRead@QualityReport@trimmedSeqLength, 364)
    expect_equal(trimmedSangerRead@QualityReport@trimmedStartPos, 34)
    expect_equal(trimmedSangerRead@QualityReport@trimmedFinishPos, 398)
    expect_equal(trimmedSangerRead@QualityReport@rawMeanQualityScore, 52.87607, tolerance=1e-6)
    expect_equal(trimmedSangerRead@QualityReport@trimmedMeanQualityScore, 58.99451, tolerance=1e-6)
    expect_equal(trimmedSangerRead@QualityReport@rawMinQualityScore, 1)
    expect_equal(trimmedSangerRead@QualityReport@trimmedMinQualityScore, 15)
    expect_equal(trimmedSangerRead@QualityReport@remainingRatio, 0.5185185, tolerance=1e-6)
})

test_that("SangerRead update quality trimming parameters 2 (bigger M1)", {
    trimmedSangerRead2 <- updateQualityParam(sangerRead,
                                            TrimmingMethod       = "M1",
                                            M1TrimmingCutoff     = 0.1,
                                            M2CutoffQualityScore = NULL,
                                            M2SlidingWindowSize  = NULL)
    expect_type(trimmedSangerRead2, "S4")
    expect_s4_class(trimmedSangerRead2, "SangerRead")
    expect_equal(trimmedSangerRead2@readFeature, "Forward Read")
    expect_equal(basename(trimmedSangerRead2@readFileName), "ACHLO006-09[LCO1490_t1,HCO2198_t1]_F_1.ab1")
    expect_equal(as.character(trimmedSangerRead2@primarySeq), "CACTTTATATTTTATTCTGGGCGTCTGAGCAGGAATGGTTGGAGCCGGTATAAGACTTCTAATTCGAATCGAGCTAAGACAACCAGGAGCGTTCCTGGGCAGAGACCAACTATACAATACTATCGTTACTGCACACGCATTTGTAATAATCTTCTTTCTAGTAATGCCTGTATTCATCGGGGGATTCGGAAACTGGCTTTTACCTTTAATACTTGGAGCCCCCGATATAGCATTCCCTCGACTCAACAACATGAGATTCTGACTACTTCCCCCATCACTGATCCTTTTAGTGTCCTCTGCGGCGGTAGAAAAAGGCGCTGGTACGGGGTGAACTGTTTATCCGCCTCTAGCAAGAAATCTTGCCCACGCAGGCCCGTCTGTAGATTTAGCCATCTTTTCCCTTCATTTAGCGGGTGCGTCTTCTATTCTAGGGGCTATTAATTTTATCACCACAGTTATTAATATGCGTTGAAGAGGATTACGTCTTGAACGAATTCCCCTGTTTGTCTGAGCTGTGCTAATTACAGTTGTTCTTCTACTTCTATCTTTACCAGTGCTAGCAGGTGCCATTACCATACTTCTTACCGACCGAAACCTCAATACTTCATTCTTTGATCCTGCCGGTGGTGGAGACCCCATCCTCTACTAGCACTTATTCTGATTTTTAGATCACCCTGATGTTGAGTCATACTGAATTCCTGA")
    expect_equal(as.character(trimmedSangerRead2@secondarySeq), "TACTTTATATTTTATTCTGGGCGTCTGAGCAGGAATGGTTGGAGCCGGTATAAGACTTCTAATTCGAATCGAGCTAAGACAACCAGGAGCGTTCCTGGGCAGAGACCAACTATACAATACTATCGTTACTGCACACGCATTTGTAATAATCTTCTTTCTAGTAATGCCTGTATTCATCGGGGGATTCGGAAACTGGCTTTTACCTTTAATACTTGGAGCCCCCGATATAGCATTCCCTCGACTCAACAACATGAGATTCTGACTACTTCCCCCATCACTGATCCTTTTAGTGTCCTCTGCGGCGGTAGAAAAAGGCGCTGGTACGGGGTGAACTGTTTATCCGCCTCTAGCAAGAAATCTTGCCCACGCAGGCCCGTCTGTAGATTTAGCCATCTTTTCCCTTCATTTAGCGGGTGCGTCTTCTATTCTAGGGGCTATTAATTTTATCACCACAGTTATTAATATGCGTTGAAGAGGATTACGTCTTGAACGAATTCCCCTGTTTGTCTGAGCTGTGCTAATTACAGTTGTTCTTCTACTTCTATCTTTACCAGTGCTAGCAGGTGCCATTACCATACTTCTTACCGACCGAAACCTCAATACTTCATTCTTTGATCCTGCCGGTGGTGGAGACCCCATCCTCTACTAGCACTTATTCTGATTTTTCGATCACCCTGATGTTGAGTCATAGTGAATTCCTGA")
    expect_equal(as.character(trimmedSangerRead2@primaryAASeqS1), "HFIFYSGRLSRNGWSRYKTSNSNRAKTTRSVPGQRPTIQYYRYCTRICNNLLSSNACIHRGIRKLAFTFNTWSPRYSIPSTQQHEILTTSPITDPFSVLCGGRKRRWYGVNCLSASSKKSCPRRPVCRFSHLFPSFSGCVFYSRGY*FYHHSY*YALKRITS*TNSPVCLSCANYSCSSTSIFTSASRCHYHTSYRPKPQYFIL*SCRWWRPHPLLALILIFRSP*C*VILNS*")
    expect_equal(as.character(trimmedSangerRead2@primaryAASeqS2), "TLYFILGV*AGMVGAGIRLLIRIELRQPGAFLGRDQLYNTIVTAHAFVIIFFLVMPVFIGGFGNWLLPLILGAPDIAFPRLNNMRF*LLPPSLILLVSSAAVEKGAGTG*TVYPPLARNLAHAGPSVDLAIFSLHLAGASSILGAINFITTVINMR*RGLRLERIPLFV*AVLITVVLLLLSLPVLAGAITILLTDRNLNTSFFDPAGGGDPILY*HLF*FLDHPDVESY*IP")
    expect_equal(as.character(trimmedSangerRead2@primaryAASeqS3), "LYILFWASEQEWLEPV*DF*FESS*DNQERSWAETNYTILSLLHTHL**SSF**CLYSSGDSETGFYL*YLEPPI*HSLDSTT*DSDYFPHH*SF*CPLRR*KKALVRGELFIRL*QEILPTQARL*I*PSFPFI*RVRLLF*GLLILSPQLLICVEEDYVLNEFPCLSELC*LQLFFYFYLYQC*QVPLPYFLPTETSILHSLILPVVETPSSTSTYSDF*ITLMLSHTEFL")
    expect_equal(trimmedSangerRead2@ChromatogramParam@baseNumPerRow, 100)
    expect_equal(trimmedSangerRead2@ChromatogramParam@heightPerRow, 200)
    expect_equal(trimmedSangerRead2@ChromatogramParam@signalRatioCutoff, 0.33)
    expect_equal(trimmedSangerRead2@ChromatogramParam@showTrimmed, TRUE)

    expect_equal(trimmedSangerRead2@QualityReport@TrimmingMethod, "M1")
    expect_equal(trimmedSangerRead2@QualityReport@M1TrimmingCutoff, 0.1, tolerance=1e-10)
    expect_equal(trimmedSangerRead2@QualityReport@M2CutoffQualityScore, NULL)
    expect_equal(trimmedSangerRead2@QualityReport@M2SlidingWindowSize, NULL)
    expect_equal(length(trimmedSangerRead2@QualityReport@qualityPhredScoresRaw), 703)
    expect_equal(length(trimmedSangerRead2@QualityReport@qualityPhredScores), 702)
    expect_equal(length(trimmedSangerRead2@QualityReport@qualityBaseScores), 702)

    expect_equal(trimmedSangerRead2@QualityReport@rawSeqLength, 702)
    expect_equal(trimmedSangerRead2@QualityReport@trimmedSeqLength, 694)
    expect_equal(trimmedSangerRead2@QualityReport@trimmedStartPos, 5)
    expect_equal(trimmedSangerRead2@QualityReport@trimmedFinishPos, 699)
    expect_equal(trimmedSangerRead2@QualityReport@rawMeanQualityScore, 52.87607, tolerance=1e-6)
    expect_equal(trimmedSangerRead2@QualityReport@trimmedMeanQualityScore, 53.38329, tolerance=1e-6)
    expect_equal(trimmedSangerRead2@QualityReport@rawMinQualityScore, 1)
    expect_equal(trimmedSangerRead2@QualityReport@trimmedMinQualityScore, 5)
    expect_equal(trimmedSangerRead2@QualityReport@remainingRatio, 0.988604, tolerance=1e-6)
})

test_that("SangerRead update quality trimming parameters 3 (equal to 1). Only trim first base", {
    # First base will always be trimmed ~
    trimmedSangerRead3 <- updateQualityParam(sangerRead,
                                             TrimmingMethod       = "M1",
                                             M1TrimmingCutoff     = 1,
                                             M2CutoffQualityScore = NULL,
                                             M2SlidingWindowSize  = NULL)
    expect_type(trimmedSangerRead3, "S4")
    expect_s4_class(trimmedSangerRead3, "SangerRead")
    expect_equal(trimmedSangerRead3@readFeature, "Forward Read")
    expect_equal(basename(trimmedSangerRead3@readFileName), "ACHLO006-09[LCO1490_t1,HCO2198_t1]_F_1.ab1")
    expect_equal(as.character(trimmedSangerRead3@primarySeq), "CACTTTATATTTTATTCTGGGCGTCTGAGCAGGAATGGTTGGAGCCGGTATAAGACTTCTAATTCGAATCGAGCTAAGACAACCAGGAGCGTTCCTGGGCAGAGACCAACTATACAATACTATCGTTACTGCACACGCATTTGTAATAATCTTCTTTCTAGTAATGCCTGTATTCATCGGGGGATTCGGAAACTGGCTTTTACCTTTAATACTTGGAGCCCCCGATATAGCATTCCCTCGACTCAACAACATGAGATTCTGACTACTTCCCCCATCACTGATCCTTTTAGTGTCCTCTGCGGCGGTAGAAAAAGGCGCTGGTACGGGGTGAACTGTTTATCCGCCTCTAGCAAGAAATCTTGCCCACGCAGGCCCGTCTGTAGATTTAGCCATCTTTTCCCTTCATTTAGCGGGTGCGTCTTCTATTCTAGGGGCTATTAATTTTATCACCACAGTTATTAATATGCGTTGAAGAGGATTACGTCTTGAACGAATTCCCCTGTTTGTCTGAGCTGTGCTAATTACAGTTGTTCTTCTACTTCTATCTTTACCAGTGCTAGCAGGTGCCATTACCATACTTCTTACCGACCGAAACCTCAATACTTCATTCTTTGATCCTGCCGGTGGTGGAGACCCCATCCTCTACTAGCACTTATTCTGATTTTTAGATCACCCTGATGTTGAGTCATACTGAATTCCTGA")
    expect_equal(as.character(trimmedSangerRead3@secondarySeq), "TACTTTATATTTTATTCTGGGCGTCTGAGCAGGAATGGTTGGAGCCGGTATAAGACTTCTAATTCGAATCGAGCTAAGACAACCAGGAGCGTTCCTGGGCAGAGACCAACTATACAATACTATCGTTACTGCACACGCATTTGTAATAATCTTCTTTCTAGTAATGCCTGTATTCATCGGGGGATTCGGAAACTGGCTTTTACCTTTAATACTTGGAGCCCCCGATATAGCATTCCCTCGACTCAACAACATGAGATTCTGACTACTTCCCCCATCACTGATCCTTTTAGTGTCCTCTGCGGCGGTAGAAAAAGGCGCTGGTACGGGGTGAACTGTTTATCCGCCTCTAGCAAGAAATCTTGCCCACGCAGGCCCGTCTGTAGATTTAGCCATCTTTTCCCTTCATTTAGCGGGTGCGTCTTCTATTCTAGGGGCTATTAATTTTATCACCACAGTTATTAATATGCGTTGAAGAGGATTACGTCTTGAACGAATTCCCCTGTTTGTCTGAGCTGTGCTAATTACAGTTGTTCTTCTACTTCTATCTTTACCAGTGCTAGCAGGTGCCATTACCATACTTCTTACCGACCGAAACCTCAATACTTCATTCTTTGATCCTGCCGGTGGTGGAGACCCCATCCTCTACTAGCACTTATTCTGATTTTTCGATCACCCTGATGTTGAGTCATAGTGAATTCCTGA")
    expect_equal(as.character(trimmedSangerRead3@primaryAASeqS1), "HFIFYSGRLSRNGWSRYKTSNSNRAKTTRSVPGQRPTIQYYRYCTRICNNLLSSNACIHRGIRKLAFTFNTWSPRYSIPSTQQHEILTTSPITDPFSVLCGGRKRRWYGVNCLSASSKKSCPRRPVCRFSHLFPSFSGCVFYSRGY*FYHHSY*YALKRITS*TNSPVCLSCANYSCSSTSIFTSASRCHYHTSYRPKPQYFIL*SCRWWRPHPLLALILIFRSP*C*VILNS*")
    expect_equal(as.character(trimmedSangerRead3@primaryAASeqS2), "TLYFILGV*AGMVGAGIRLLIRIELRQPGAFLGRDQLYNTIVTAHAFVIIFFLVMPVFIGGFGNWLLPLILGAPDIAFPRLNNMRF*LLPPSLILLVSSAAVEKGAGTG*TVYPPLARNLAHAGPSVDLAIFSLHLAGASSILGAINFITTVINMR*RGLRLERIPLFV*AVLITVVLLLLSLPVLAGAITILLTDRNLNTSFFDPAGGGDPILY*HLF*FLDHPDVESY*IP")
    expect_equal(as.character(trimmedSangerRead3@primaryAASeqS3), "LYILFWASEQEWLEPV*DF*FESS*DNQERSWAETNYTILSLLHTHL**SSF**CLYSSGDSETGFYL*YLEPPI*HSLDSTT*DSDYFPHH*SF*CPLRR*KKALVRGELFIRL*QEILPTQARL*I*PSFPFI*RVRLLF*GLLILSPQLLICVEEDYVLNEFPCLSELC*LQLFFYFYLYQC*QVPLPYFLPTETSILHSLILPVVETPSSTSTYSDF*ITLMLSHTEFL")
    expect_equal(trimmedSangerRead3@ChromatogramParam@baseNumPerRow, 100)
    expect_equal(trimmedSangerRead3@ChromatogramParam@heightPerRow, 200)
    expect_equal(trimmedSangerRead3@ChromatogramParam@signalRatioCutoff, 0.33)
    expect_equal(trimmedSangerRead3@ChromatogramParam@showTrimmed, TRUE)

    expect_equal(trimmedSangerRead3@QualityReport@TrimmingMethod, "M1")
    expect_equal(trimmedSangerRead3@QualityReport@M1TrimmingCutoff, 1, tolerance=1e-10)
    expect_equal(trimmedSangerRead3@QualityReport@M2CutoffQualityScore, NULL)
    expect_equal(trimmedSangerRead3@QualityReport@M2SlidingWindowSize, NULL)
    expect_equal(length(trimmedSangerRead3@QualityReport@qualityPhredScoresRaw), 703)
    expect_equal(length(trimmedSangerRead3@QualityReport@qualityPhredScores), 702)
    expect_equal(length(trimmedSangerRead3@QualityReport@qualityBaseScores), 702)

    expect_equal(trimmedSangerRead3@QualityReport@rawSeqLength, 702)
    expect_equal(trimmedSangerRead3@QualityReport@trimmedSeqLength, 701)
    expect_equal(trimmedSangerRead3@QualityReport@trimmedStartPos, 1)
    expect_equal(trimmedSangerRead3@QualityReport@trimmedFinishPos, 702)
    expect_equal(trimmedSangerRead3@QualityReport@rawMeanQualityScore, 52.87607, tolerance=1e-6)
    expect_equal(trimmedSangerRead3@QualityReport@trimmedMeanQualityScore, 52.94864, tolerance=1e-6)
    expect_equal(trimmedSangerRead3@QualityReport@rawMinQualityScore, 1)
    expect_equal(trimmedSangerRead3@QualityReport@trimmedMinQualityScore, 1)
    expect_equal(trimmedSangerRead3@QualityReport@remainingRatio, 0.9985755, tolerance=1e-6)
})

test_that("SangerRead update quality trimming parameters 4 (very small). All trimmed", {
    # First base will always be trimmed ~
    trimmedSangerRead4 <- updateQualityParam(sangerRead,
                                             TrimmingMethod       = "M1",
                                             M1TrimmingCutoff     =  0.000000001,
                                             M2CutoffQualityScore = NULL,
                                             M2SlidingWindowSize  = NULL)
    expect_type(trimmedSangerRead4, "S4")
    expect_s4_class(trimmedSangerRead4, "SangerRead")
    expect_equal(trimmedSangerRead4@readFeature, "Forward Read")
    expect_equal(basename(trimmedSangerRead4@readFileName), "ACHLO006-09[LCO1490_t1,HCO2198_t1]_F_1.ab1")
    expect_equal(as.character(trimmedSangerRead4@primarySeq), "CACTTTATATTTTATTCTGGGCGTCTGAGCAGGAATGGTTGGAGCCGGTATAAGACTTCTAATTCGAATCGAGCTAAGACAACCAGGAGCGTTCCTGGGCAGAGACCAACTATACAATACTATCGTTACTGCACACGCATTTGTAATAATCTTCTTTCTAGTAATGCCTGTATTCATCGGGGGATTCGGAAACTGGCTTTTACCTTTAATACTTGGAGCCCCCGATATAGCATTCCCTCGACTCAACAACATGAGATTCTGACTACTTCCCCCATCACTGATCCTTTTAGTGTCCTCTGCGGCGGTAGAAAAAGGCGCTGGTACGGGGTGAACTGTTTATCCGCCTCTAGCAAGAAATCTTGCCCACGCAGGCCCGTCTGTAGATTTAGCCATCTTTTCCCTTCATTTAGCGGGTGCGTCTTCTATTCTAGGGGCTATTAATTTTATCACCACAGTTATTAATATGCGTTGAAGAGGATTACGTCTTGAACGAATTCCCCTGTTTGTCTGAGCTGTGCTAATTACAGTTGTTCTTCTACTTCTATCTTTACCAGTGCTAGCAGGTGCCATTACCATACTTCTTACCGACCGAAACCTCAATACTTCATTCTTTGATCCTGCCGGTGGTGGAGACCCCATCCTCTACTAGCACTTATTCTGATTTTTAGATCACCCTGATGTTGAGTCATACTGAATTCCTGA")
    expect_equal(as.character(trimmedSangerRead4@secondarySeq), "TACTTTATATTTTATTCTGGGCGTCTGAGCAGGAATGGTTGGAGCCGGTATAAGACTTCTAATTCGAATCGAGCTAAGACAACCAGGAGCGTTCCTGGGCAGAGACCAACTATACAATACTATCGTTACTGCACACGCATTTGTAATAATCTTCTTTCTAGTAATGCCTGTATTCATCGGGGGATTCGGAAACTGGCTTTTACCTTTAATACTTGGAGCCCCCGATATAGCATTCCCTCGACTCAACAACATGAGATTCTGACTACTTCCCCCATCACTGATCCTTTTAGTGTCCTCTGCGGCGGTAGAAAAAGGCGCTGGTACGGGGTGAACTGTTTATCCGCCTCTAGCAAGAAATCTTGCCCACGCAGGCCCGTCTGTAGATTTAGCCATCTTTTCCCTTCATTTAGCGGGTGCGTCTTCTATTCTAGGGGCTATTAATTTTATCACCACAGTTATTAATATGCGTTGAAGAGGATTACGTCTTGAACGAATTCCCCTGTTTGTCTGAGCTGTGCTAATTACAGTTGTTCTTCTACTTCTATCTTTACCAGTGCTAGCAGGTGCCATTACCATACTTCTTACCGACCGAAACCTCAATACTTCATTCTTTGATCCTGCCGGTGGTGGAGACCCCATCCTCTACTAGCACTTATTCTGATTTTTCGATCACCCTGATGTTGAGTCATAGTGAATTCCTGA")
    expect_equal(as.character(trimmedSangerRead4@primaryAASeqS1), "HFIFYSGRLSRNGWSRYKTSNSNRAKTTRSVPGQRPTIQYYRYCTRICNNLLSSNACIHRGIRKLAFTFNTWSPRYSIPSTQQHEILTTSPITDPFSVLCGGRKRRWYGVNCLSASSKKSCPRRPVCRFSHLFPSFSGCVFYSRGY*FYHHSY*YALKRITS*TNSPVCLSCANYSCSSTSIFTSASRCHYHTSYRPKPQYFIL*SCRWWRPHPLLALILIFRSP*C*VILNS*")
    expect_equal(as.character(trimmedSangerRead4@primaryAASeqS2), "TLYFILGV*AGMVGAGIRLLIRIELRQPGAFLGRDQLYNTIVTAHAFVIIFFLVMPVFIGGFGNWLLPLILGAPDIAFPRLNNMRF*LLPPSLILLVSSAAVEKGAGTG*TVYPPLARNLAHAGPSVDLAIFSLHLAGASSILGAINFITTVINMR*RGLRLERIPLFV*AVLITVVLLLLSLPVLAGAITILLTDRNLNTSFFDPAGGGDPILY*HLF*FLDHPDVESY*IP")
    expect_equal(as.character(trimmedSangerRead4@primaryAASeqS3), "LYILFWASEQEWLEPV*DF*FESS*DNQERSWAETNYTILSLLHTHL**SSF**CLYSSGDSETGFYL*YLEPPI*HSLDSTT*DSDYFPHH*SF*CPLRR*KKALVRGELFIRL*QEILPTQARL*I*PSFPFI*RVRLLF*GLLILSPQLLICVEEDYVLNEFPCLSELC*LQLFFYFYLYQC*QVPLPYFLPTETSILHSLILPVVETPSSTSTYSDF*ITLMLSHTEFL")
    expect_equal(trimmedSangerRead4@ChromatogramParam@baseNumPerRow, 100)
    expect_equal(trimmedSangerRead4@ChromatogramParam@heightPerRow, 200)
    expect_equal(trimmedSangerRead4@ChromatogramParam@signalRatioCutoff, 0.33)
    expect_equal(trimmedSangerRead4@ChromatogramParam@showTrimmed, TRUE)

    expect_equal(trimmedSangerRead4@QualityReport@TrimmingMethod, "M1")
    expect_equal(trimmedSangerRead4@QualityReport@M1TrimmingCutoff, 0.000000001, tolerance=1e-10)
    expect_equal(trimmedSangerRead4@QualityReport@M2CutoffQualityScore, NULL)
    expect_equal(trimmedSangerRead4@QualityReport@M2SlidingWindowSize, NULL)
    expect_equal(length(trimmedSangerRead4@QualityReport@qualityPhredScoresRaw), 703)
    expect_equal(length(trimmedSangerRead4@QualityReport@qualityPhredScores), 702)
    expect_equal(length(trimmedSangerRead4@QualityReport@qualityBaseScores), 702)

    expect_equal(trimmedSangerRead4@QualityReport@rawSeqLength, 702)
    expect_equal(trimmedSangerRead4@QualityReport@trimmedSeqLength, 1)
    expect_equal(trimmedSangerRead4@QualityReport@trimmedStartPos, 1)
    expect_equal(trimmedSangerRead4@QualityReport@trimmedFinishPos, 2)
    expect_equal(trimmedSangerRead4@QualityReport@rawMeanQualityScore, 52.87607, tolerance=1e-6)
    expect_equal(trimmedSangerRead4@QualityReport@trimmedMeanQualityScore, 4, tolerance=1e-6)
    expect_equal(trimmedSangerRead4@QualityReport@rawMinQualityScore, 1)
    expect_equal(trimmedSangerRead4@QualityReport@trimmedMinQualityScore, 4)
    expect_equal(trimmedSangerRead4@QualityReport@remainingRatio, 0.001424501, tolerance=1e-6)
})


test_that("SangerRead update quality trimming parameters 5 (smaller than 0).", {
    # First base will always be trimmed ~
    expect_error(updateQualityParam(sangerRead,
                                    TrimmingMethod       = "M1",
                                    M1TrimmingCutoff     =  -0.00001,
                                    M2CutoffQualityScore = NULL,
                                    M2SlidingWindowSize  = NULL),
                 "Your input M1TrimmingCutoff is: '-1e-05' is invalid.'M1TrimmingCutoff' shouldbe between 0 and 1.", fixed = TRUE)
})

test_that("SangerRead update quality trimming parameters 6 (bigger than 1).", {
    # First base will always be trimmed ~
    expect_error(updateQualityParam(sangerRead,
                                    TrimmingMethod       = "M1",
                                    M1TrimmingCutoff     =  1.0001,
                                    M2CutoffQualityScore = NULL,
                                    M2SlidingWindowSize  = NULL),
                 "Your input M1TrimmingCutoff is: '1.0001' is invalid.'M1TrimmingCutoff' shouldbe between 0 and 1.", fixed = TRUE)
})

test_that("SangerRead update quality trimming parameters 7 (Change from M1 to M2)", {
    trimmedSangerRead7 <- updateQualityParam(sangerRead,
                                             TrimmingMethod       = "M2",
                                             M1TrimmingCutoff     = NULL,
                                             M2CutoffQualityScore = 25,
                                             M2SlidingWindowSize  = 10)
    expect_type(trimmedSangerRead7, "S4")
    expect_s4_class(trimmedSangerRead7, "SangerRead")
    expect_equal(trimmedSangerRead7@readFeature, "Forward Read")
    expect_equal(basename(trimmedSangerRead7@readFileName), "ACHLO006-09[LCO1490_t1,HCO2198_t1]_F_1.ab1")
    expect_equal(as.character(trimmedSangerRead7@primarySeq), "CACTTTATATTTTATTCTGGGCGTCTGAGCAGGAATGGTTGGAGCCGGTATAAGACTTCTAATTCGAATCGAGCTAAGACAACCAGGAGCGTTCCTGGGCAGAGACCAACTATACAATACTATCGTTACTGCACACGCATTTGTAATAATCTTCTTTCTAGTAATGCCTGTATTCATCGGGGGATTCGGAAACTGGCTTTTACCTTTAATACTTGGAGCCCCCGATATAGCATTCCCTCGACTCAACAACATGAGATTCTGACTACTTCCCCCATCACTGATCCTTTTAGTGTCCTCTGCGGCGGTAGAAAAAGGCGCTGGTACGGGGTGAACTGTTTATCCGCCTCTAGCAAGAAATCTTGCCCACGCAGGCCCGTCTGTAGATTTAGCCATCTTTTCCCTTCATTTAGCGGGTGCGTCTTCTATTCTAGGGGCTATTAATTTTATCACCACAGTTATTAATATGCGTTGAAGAGGATTACGTCTTGAACGAATTCCCCTGTTTGTCTGAGCTGTGCTAATTACAGTTGTTCTTCTACTTCTATCTTTACCAGTGCTAGCAGGTGCCATTACCATACTTCTTACCGACCGAAACCTCAATACTTCATTCTTTGATCCTGCCGGTGGTGGAGACCCCATCCTCTACTAGCACTTATTCTGATTTTTAGATCACCCTGATGTTGAGTCATACTGAATTCCTGA")
    expect_equal(as.character(trimmedSangerRead7@secondarySeq), "TACTTTATATTTTATTCTGGGCGTCTGAGCAGGAATGGTTGGAGCCGGTATAAGACTTCTAATTCGAATCGAGCTAAGACAACCAGGAGCGTTCCTGGGCAGAGACCAACTATACAATACTATCGTTACTGCACACGCATTTGTAATAATCTTCTTTCTAGTAATGCCTGTATTCATCGGGGGATTCGGAAACTGGCTTTTACCTTTAATACTTGGAGCCCCCGATATAGCATTCCCTCGACTCAACAACATGAGATTCTGACTACTTCCCCCATCACTGATCCTTTTAGTGTCCTCTGCGGCGGTAGAAAAAGGCGCTGGTACGGGGTGAACTGTTTATCCGCCTCTAGCAAGAAATCTTGCCCACGCAGGCCCGTCTGTAGATTTAGCCATCTTTTCCCTTCATTTAGCGGGTGCGTCTTCTATTCTAGGGGCTATTAATTTTATCACCACAGTTATTAATATGCGTTGAAGAGGATTACGTCTTGAACGAATTCCCCTGTTTGTCTGAGCTGTGCTAATTACAGTTGTTCTTCTACTTCTATCTTTACCAGTGCTAGCAGGTGCCATTACCATACTTCTTACCGACCGAAACCTCAATACTTCATTCTTTGATCCTGCCGGTGGTGGAGACCCCATCCTCTACTAGCACTTATTCTGATTTTTCGATCACCCTGATGTTGAGTCATAGTGAATTCCTGA")
    expect_equal(as.character(trimmedSangerRead7@primaryAASeqS1), "HFIFYSGRLSRNGWSRYKTSNSNRAKTTRSVPGQRPTIQYYRYCTRICNNLLSSNACIHRGIRKLAFTFNTWSPRYSIPSTQQHEILTTSPITDPFSVLCGGRKRRWYGVNCLSASSKKSCPRRPVCRFSHLFPSFSGCVFYSRGY*FYHHSY*YALKRITS*TNSPVCLSCANYSCSSTSIFTSASRCHYHTSYRPKPQYFIL*SCRWWRPHPLLALILIFRSP*C*VILNS*")
    expect_equal(as.character(trimmedSangerRead7@primaryAASeqS2), "TLYFILGV*AGMVGAGIRLLIRIELRQPGAFLGRDQLYNTIVTAHAFVIIFFLVMPVFIGGFGNWLLPLILGAPDIAFPRLNNMRF*LLPPSLILLVSSAAVEKGAGTG*TVYPPLARNLAHAGPSVDLAIFSLHLAGASSILGAINFITTVINMR*RGLRLERIPLFV*AVLITVVLLLLSLPVLAGAITILLTDRNLNTSFFDPAGGGDPILY*HLF*FLDHPDVESY*IP")
    expect_equal(as.character(trimmedSangerRead7@primaryAASeqS3), "LYILFWASEQEWLEPV*DF*FESS*DNQERSWAETNYTILSLLHTHL**SSF**CLYSSGDSETGFYL*YLEPPI*HSLDSTT*DSDYFPHH*SF*CPLRR*KKALVRGELFIRL*QEILPTQARL*I*PSFPFI*RVRLLF*GLLILSPQLLICVEEDYVLNEFPCLSELC*LQLFFYFYLYQC*QVPLPYFLPTETSILHSLILPVVETPSSTSTYSDF*ITLMLSHTEFL")
    expect_equal(trimmedSangerRead7@ChromatogramParam@baseNumPerRow, 100)
    expect_equal(trimmedSangerRead7@ChromatogramParam@heightPerRow, 200)
    expect_equal(trimmedSangerRead7@ChromatogramParam@signalRatioCutoff, 0.33)
    expect_equal(trimmedSangerRead7@ChromatogramParam@showTrimmed, TRUE)

    expect_equal(trimmedSangerRead7@QualityReport@TrimmingMethod, "M2")
    expect_equal(trimmedSangerRead7@QualityReport@M1TrimmingCutoff, NULL, tolerance=1e-10)
    expect_equal(trimmedSangerRead7@QualityReport@M2CutoffQualityScore, 25, tolerance=1e-10)
    expect_equal(trimmedSangerRead7@QualityReport@M2SlidingWindowSize, 10, tolerance=1e-10)
    expect_equal(length(trimmedSangerRead7@QualityReport@qualityPhredScoresRaw), 703)
    expect_equal(length(trimmedSangerRead7@QualityReport@qualityPhredScores), 702)
    expect_equal(length(trimmedSangerRead7@QualityReport@qualityBaseScores), 702)

    expect_equal(trimmedSangerRead7@QualityReport@rawSeqLength, 702)
    expect_equal(trimmedSangerRead7@QualityReport@trimmedSeqLength, 681)
    expect_equal(trimmedSangerRead7@QualityReport@trimmedStartPos, 7)
    expect_equal(trimmedSangerRead7@QualityReport@trimmedFinishPos, 688)
    expect_equal(trimmedSangerRead7@QualityReport@rawMeanQualityScore, 52.87607, tolerance=1e-6)
    expect_equal(trimmedSangerRead7@QualityReport@trimmedMeanQualityScore, 54.05286, tolerance=1e-6)
    expect_equal(trimmedSangerRead7@QualityReport@rawMinQualityScore, 1)
    expect_equal(trimmedSangerRead7@QualityReport@trimmedMinQualityScore, 8)
    expect_equal(trimmedSangerRead7@QualityReport@remainingRatio, 0.9700855, tolerance=1e-6)
})

test_that("SangerRead update quality trimming parameters 8 (Update from M1 to M2)", {
    trimmedSangerRead <- updateQualityParam(sangerRead,
                                            TrimmingMethod       = "M1",
                                            M1TrimmingCutoff     = 0.000001,
                                            M2CutoffQualityScore = NULL,
                                            M2SlidingWindowSize  = NULL)
    trimmedSangerRead8 <- updateQualityParam(trimmedSangerRead,
                                             TrimmingMethod       = "M2",
                                             M1TrimmingCutoff     = NULL,
                                             M2CutoffQualityScore = 30,
                                             M2SlidingWindowSize  = 12)
    expect_type(trimmedSangerRead8, "S4")
    expect_s4_class(trimmedSangerRead8, "SangerRead")
    expect_equal(trimmedSangerRead8@readFeature, "Forward Read")
    expect_equal(basename(trimmedSangerRead8@readFileName), "ACHLO006-09[LCO1490_t1,HCO2198_t1]_F_1.ab1")
    expect_equal(as.character(trimmedSangerRead8@primarySeq), "CACTTTATATTTTATTCTGGGCGTCTGAGCAGGAATGGTTGGAGCCGGTATAAGACTTCTAATTCGAATCGAGCTAAGACAACCAGGAGCGTTCCTGGGCAGAGACCAACTATACAATACTATCGTTACTGCACACGCATTTGTAATAATCTTCTTTCTAGTAATGCCTGTATTCATCGGGGGATTCGGAAACTGGCTTTTACCTTTAATACTTGGAGCCCCCGATATAGCATTCCCTCGACTCAACAACATGAGATTCTGACTACTTCCCCCATCACTGATCCTTTTAGTGTCCTCTGCGGCGGTAGAAAAAGGCGCTGGTACGGGGTGAACTGTTTATCCGCCTCTAGCAAGAAATCTTGCCCACGCAGGCCCGTCTGTAGATTTAGCCATCTTTTCCCTTCATTTAGCGGGTGCGTCTTCTATTCTAGGGGCTATTAATTTTATCACCACAGTTATTAATATGCGTTGAAGAGGATTACGTCTTGAACGAATTCCCCTGTTTGTCTGAGCTGTGCTAATTACAGTTGTTCTTCTACTTCTATCTTTACCAGTGCTAGCAGGTGCCATTACCATACTTCTTACCGACCGAAACCTCAATACTTCATTCTTTGATCCTGCCGGTGGTGGAGACCCCATCCTCTACTAGCACTTATTCTGATTTTTAGATCACCCTGATGTTGAGTCATACTGAATTCCTGA")
    expect_equal(as.character(trimmedSangerRead8@secondarySeq), "TACTTTATATTTTATTCTGGGCGTCTGAGCAGGAATGGTTGGAGCCGGTATAAGACTTCTAATTCGAATCGAGCTAAGACAACCAGGAGCGTTCCTGGGCAGAGACCAACTATACAATACTATCGTTACTGCACACGCATTTGTAATAATCTTCTTTCTAGTAATGCCTGTATTCATCGGGGGATTCGGAAACTGGCTTTTACCTTTAATACTTGGAGCCCCCGATATAGCATTCCCTCGACTCAACAACATGAGATTCTGACTACTTCCCCCATCACTGATCCTTTTAGTGTCCTCTGCGGCGGTAGAAAAAGGCGCTGGTACGGGGTGAACTGTTTATCCGCCTCTAGCAAGAAATCTTGCCCACGCAGGCCCGTCTGTAGATTTAGCCATCTTTTCCCTTCATTTAGCGGGTGCGTCTTCTATTCTAGGGGCTATTAATTTTATCACCACAGTTATTAATATGCGTTGAAGAGGATTACGTCTTGAACGAATTCCCCTGTTTGTCTGAGCTGTGCTAATTACAGTTGTTCTTCTACTTCTATCTTTACCAGTGCTAGCAGGTGCCATTACCATACTTCTTACCGACCGAAACCTCAATACTTCATTCTTTGATCCTGCCGGTGGTGGAGACCCCATCCTCTACTAGCACTTATTCTGATTTTTCGATCACCCTGATGTTGAGTCATAGTGAATTCCTGA")
    expect_equal(as.character(trimmedSangerRead8@primaryAASeqS1), "HFIFYSGRLSRNGWSRYKTSNSNRAKTTRSVPGQRPTIQYYRYCTRICNNLLSSNACIHRGIRKLAFTFNTWSPRYSIPSTQQHEILTTSPITDPFSVLCGGRKRRWYGVNCLSASSKKSCPRRPVCRFSHLFPSFSGCVFYSRGY*FYHHSY*YALKRITS*TNSPVCLSCANYSCSSTSIFTSASRCHYHTSYRPKPQYFIL*SCRWWRPHPLLALILIFRSP*C*VILNS*")
    expect_equal(as.character(trimmedSangerRead8@primaryAASeqS2), "TLYFILGV*AGMVGAGIRLLIRIELRQPGAFLGRDQLYNTIVTAHAFVIIFFLVMPVFIGGFGNWLLPLILGAPDIAFPRLNNMRF*LLPPSLILLVSSAAVEKGAGTG*TVYPPLARNLAHAGPSVDLAIFSLHLAGASSILGAINFITTVINMR*RGLRLERIPLFV*AVLITVVLLLLSLPVLAGAITILLTDRNLNTSFFDPAGGGDPILY*HLF*FLDHPDVESY*IP")
    expect_equal(as.character(trimmedSangerRead8@primaryAASeqS3), "LYILFWASEQEWLEPV*DF*FESS*DNQERSWAETNYTILSLLHTHL**SSF**CLYSSGDSETGFYL*YLEPPI*HSLDSTT*DSDYFPHH*SF*CPLRR*KKALVRGELFIRL*QEILPTQARL*I*PSFPFI*RVRLLF*GLLILSPQLLICVEEDYVLNEFPCLSELC*LQLFFYFYLYQC*QVPLPYFLPTETSILHSLILPVVETPSSTSTYSDF*ITLMLSHTEFL")
    expect_equal(trimmedSangerRead8@ChromatogramParam@baseNumPerRow, 100)
    expect_equal(trimmedSangerRead8@ChromatogramParam@heightPerRow, 200)
    expect_equal(trimmedSangerRead8@ChromatogramParam@signalRatioCutoff, 0.33)
    expect_equal(trimmedSangerRead8@ChromatogramParam@showTrimmed, TRUE)

    expect_equal(trimmedSangerRead8@QualityReport@TrimmingMethod, "M2")
    expect_equal(trimmedSangerRead8@QualityReport@M1TrimmingCutoff, NULL)
    expect_equal(trimmedSangerRead8@QualityReport@M2CutoffQualityScore, 30, tolerance=1e-10)
    expect_equal(trimmedSangerRead8@QualityReport@M2SlidingWindowSize, 12, tolerance=1e-10)
    expect_equal(length(trimmedSangerRead8@QualityReport@qualityPhredScoresRaw), 703)
    expect_equal(length(trimmedSangerRead8@QualityReport@qualityPhredScores), 702)
    expect_equal(length(trimmedSangerRead8@QualityReport@qualityBaseScores), 702)

    expect_equal(trimmedSangerRead8@QualityReport@rawSeqLength, 702)
    expect_equal(trimmedSangerRead8@QualityReport@trimmedSeqLength, 667)
    expect_equal(trimmedSangerRead8@QualityReport@trimmedStartPos, 18)
    expect_equal(trimmedSangerRead8@QualityReport@trimmedFinishPos, 685)
    expect_equal(trimmedSangerRead8@QualityReport@rawMeanQualityScore, 52.87607, tolerance=1e-6)
    expect_equal(trimmedSangerRead8@QualityReport@trimmedMeanQualityScore, 54.63568, tolerance=1e-6)
    expect_equal(trimmedSangerRead8@QualityReport@rawMinQualityScore, 1)
    expect_equal(trimmedSangerRead8@QualityReport@trimmedMinQualityScore, 12)
    expect_equal(trimmedSangerRead8@QualityReport@remainingRatio, 0.9501425, tolerance=1e-6)
})

test_that("SangerRead update quality trimming parameters 9 (M2CutoffQualityScore bigger than threashold)", {
    expect_error(updateQualityParam(sangerRead,
                                    TrimmingMethod       = "M2",
                                    M1TrimmingCutoff     = NULL,
                                    M2CutoffQualityScore = 61,
                                    M2SlidingWindowSize  = 41),
                 "\nYour input M2CutoffQualityScore is: '61' is invalid.'M2CutoffQualityScore' shouldbe between 0 and 60.\n\nYour input M2SlidingWindowSize is: '41' is invalid.'M2SlidingWindowSize' shouldbe between 0 and 40.\n", fixed = TRUE)
})

test_that("SangerRead update quality trimming parameters 10 (M2CutoffQualityScore smaller than threashold)", {
    expect_error(updateQualityParam(sangerRead,
                                    TrimmingMethod       = "M2",
                                    M1TrimmingCutoff     = NULL,
                                    M2CutoffQualityScore = -1,
                                    M2SlidingWindowSize  = -1),
                 "\nYour input M2CutoffQualityScore is: '-1' is invalid.'M2CutoffQualityScore' shouldbe between 0 and 60.\n\nYour input M2SlidingWindowSize is: '-1' is invalid.'M2SlidingWindowSize' shouldbe between 0 and 40.\n", fixed = TRUE)
})
