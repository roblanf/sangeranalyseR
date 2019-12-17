test_that("sangerContig vs sangerContigFa Initial test", {
    expect_equal(sangerContig@contigName, sangerContigFa@contigName)
    expect_equal(as.character(sangerContig@contigSeq), as.character(sangerContigFa@contigSeq))

    expect_equal(sangerContig2@contigName, sangerContigFa2@contigName)
    expect_equal(as.character(sangerContig2@contigSeq), as.character(sangerContigFa2@contigSeq))
})


