test_that("sangerReadF vs sangerReadFFa Initial test", {
    expect_equal(sangerReadF@readFeature, sangerReadFFa@readFeature)
    expect_equal(substr(as.character(sangerReadF@primarySeq), sangerReadF@QualityReport@trimmedStartPos+1, sangerReadF@QualityReport@trimmedFinishPos), as.character(sangerReadFFa@primarySeq))
})

test_that("sangerReadR vs sangerReadRFa Initial test", {
    expect_equal(sangerReadR@readFeature, sangerReadRFa@readFeature)
    expect_equal(substr(as.character(sangerReadR@primarySeq), sangerReadR@QualityReport@trimmedStartPos+1, sangerReadR@QualityReport@trimmedFinishPos), as.character(sangerReadRFa@primarySeq))
})
