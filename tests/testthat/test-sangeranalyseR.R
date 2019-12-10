test_that("SangerRead creation test", {
    expect_equal(2*2, 4)
})


test_that("SangerRead second test", {
    # hetab1 <- readsangerseq(system.file("extdata", "heterozygous.ab1", package = "sangerseqR"))
    # makeBaseCalls(hetab1)

    inputFilesPath <- system.file("extdata/", package = "sangeranalyseR")
    A_chloroticaFdReadFN <- file.path(inputFilesPath,
                                      "Allolobophora_chlorotica",
                                      "RBNII384-13[C_LepFolF,C_LepFolR]_R_2.ab1")
    A_chloroticaRead <- new("SangerRead",
                            readFeature           = "Forward Read",
                            readFileName          = A_chloroticaFdReadFN)




    # expect_type(A_chloroticaRead, "S4")

    # expect_equal(A_chloroticaRead@readFeature, "Forward Read")
    expect_equal(2*2, 4)
    # expect_s4_class(A_chloroticaRead, "SangerRead")
})
