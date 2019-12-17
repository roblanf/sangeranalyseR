test_that("sangerContig vs sangerContigFa Initial test", {
    expect_equal(sangerContig@contigName, sangerContigFa@contigName)
    expect_equal(as.character(sangerContig@contigSeq), as.character(sangerContigFa@contigSeq))
})


