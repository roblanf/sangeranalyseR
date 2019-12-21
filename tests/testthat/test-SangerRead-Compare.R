test_that("sangerReadF vs sangerReadFFa Initial test", {
    expect_equal(sangerReadF@readFeature, sangerReadFFa@readFeature)
    expect_equal(substr(as.character(sangerReadF@primarySeq), sangerReadF@QualityReport@trimmedStartPos+1, sangerReadF@QualityReport@trimmedFinishPos), as.character(sangerReadFFa@primarySeq))
    expect_equal(as.character(sangerReadF@primaryAASeqS1), as.character(sangerReadFFa@primaryAASeqS1))
    expect_equal(as.character(sangerReadF@primaryAASeqS2), as.character(sangerReadFFa@primaryAASeqS2))
    expect_equal(as.character(sangerReadF@primaryAASeqS3), as.character(sangerReadFFa@primaryAASeqS3))
})

test_that("sangerReadR vs sangerReadRFa Initial test", {
    expect_equal(sangerReadR@readFeature, sangerReadRFa@readFeature)
    expect_equal(substr(as.character(sangerReadR@primarySeq), sangerReadR@QualityReport@trimmedStartPos+1, sangerReadR@QualityReport@trimmedFinishPos), as.character(sangerReadRFa@primarySeq))
    expect_equal(as.character(sangerReadR@primaryAASeqS1), as.character(sangerReadRFa@primaryAASeqS1))
    expect_equal(as.character(sangerReadR@primaryAASeqS2), as.character(sangerReadRFa@primaryAASeqS2))
    expect_equal(as.character(sangerReadR@primaryAASeqS3), as.character(sangerReadRFa@primaryAASeqS3))
})
