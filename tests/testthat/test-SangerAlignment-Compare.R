test_that("sangerAlignment vs sangerAlignmentFa Initial test", {
    expect_equal(as.character(sangerAlignment@contigsConsensus), as.character(sangerAlignmentFa@contigsConsensus))
})


