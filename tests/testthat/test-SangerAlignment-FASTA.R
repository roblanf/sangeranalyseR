### ============================================================================
### SangerAlignment Initialization test
### ============================================================================
test_that("sangerAlignmentFa Initial test", {
    expect_type(sangerAlignmentFa, "S4")
    expect_s4_class(sangerAlignmentFa, "SangerAlignment")
    expect_equal(sangerAlignmentFa@inputSource, "FASTA")
    expect_equal(basename(sangerAlignmentFa@FASTA_File), "Sanger_all_reads.fa")
    expect_equal(sangerAlignmentFa@REGEX_SuffixForward, "_[0-9]*_F")
    expect_equal(sangerAlignmentFa@REGEX_SuffixReverse, "_[0-9]*_R")
    expect_equal(sangerAlignmentFa@trimmingMethodSA, "")
    expect_equal(sangerAlignmentFa@minFractionCallSA, 0.5)
    expect_equal(sangerAlignmentFa@maxFractionLostSA, 0.5)
    expect_equal(sangerAlignmentFa@refAminoAcidSeq, "SRQWLFSTNHKDIGTLYFIFGAWAGMVGTSLSILIRAELGHPGALIGDDQIYNVIVTAHAFIMIFFMVMPIMIGGFGNWLVPLMLGAPDMAFPRMNNMSFWLLPPALSLLLVSSMVENGAGTGWTVYPPLSAGIAHGGASVDLAIFSLHLAGISSILGAVNFITTVINMRSTGISLDRMPLFVWSVVITALLLLLSLPVLAGAITMLLTDRNLNTSFFDPAGGGDPILYQHLFWFFGHPEVYILILPGFGMISHIISQESGKKETFGSLGMIYAMLAIGLLGFIVWAHHMFTVGMDVDTRAYFTSATMIIAVPTGIKIFSWLATLHGTQLSYSPAILWALGFVFLFTVGGLTGVVLANSSVDIILHDTYYVVAHFHYVLSMGAVFAIMAGFIHWYPLFTGLTLNNKWLKSHFIIMFIGVNLTFFPQHFLGLAGMPRRYSDYPDAYTTWNIVSTIGSTISLLGILFFFFIIWESLVSQRQVIYPIQLNSSIEWYQNTPPAEHSYSELPLLTN")
    expect_equal(as.character(sangerAlignmentFa@contigsConsensus), "TTATAYTTTATTYTRGGCGTCTGAGCAGGAATGGTTGGAGCYGGTATAAGACTYCTAATTCGAATYGAGCTAAGACARCCRGGAGCRTTCCTRGGMAGRGAYCAACTMTAYAATACTATYGTWACTGCWCACGCATTTGTAATAATYTTCTTTCTAGTAATRCCTGTATTYATYGGGGGRTTCGGWAAYTGRCTTYTACCTTTAATACTTGGAGCCCCYGAYATRGCATTCCCWCGACTYAACAACATRAGATTCTGACTMCTTCCCCCATCACTRATCCTTYTAGTGTCCTCTGCKGCRGTAGAAAAAGGCGCTGGWACKGGRTGAACTGTTTATCCGCCYCTAGCAAGAAATMTTGCYCAYGCMGGCCCRTCTGTAGAYTTAGCYATYTTTTCYCTTCATTTAGCGGGTGCKTCWTCWATYYTAGGGGCYATTAATTTTATYACYACWGTTATTAAYATGCGWTGAAGAGGMTTACGWCTTGAACGAATYCCMYTRTTYGTYTGAGCYGTRCTAATTACAGTKGTTCTTCTACTYCTATCYTTACCAGTGYTAGCMGGTGCMATTACYATACTWCTTACCGAYCGAAAYCTCAATACYTCMTTCTTTGATCCTGCYGGTGGTGGAGAYCCCATCCTCTACCAACACTTATTCTGATTTTTTGGTCACCCTGAG")

    expect_equal(sangerAlignmentFa@contigsTree$tip.label[1], "Achl_ACHLO006-09")
    expect_equal(sangerAlignmentFa@contigsTree$tip.label[2], "Achl_ACHLO007-09")
    expect_equal(sangerAlignmentFa@contigsTree$tip.label[3], "Achl_ACHLO040-09")
    expect_equal(sangerAlignmentFa@contigsTree$tip.label[4], "Achl_ACHLO041-09")
})

### ============================================================================
### SangerAlignment Functions test
### ============================================================================
# test_that("sangerAlignmentFa update quality trimming parameters test 1 - M1", {
#     expect_message(updateQualityParam(sangerAlignmentFa),
#                    paste("SangerAlignment with 'FASTA' inputSource",
#                          "cannot update quality parameters"))
# })
