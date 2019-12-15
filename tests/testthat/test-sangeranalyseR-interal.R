test_that("Get processor number test", {
    expect_equal(getProcessors(10), 1)
    expect_equal(getProcessors(15), 1)
    expect_equal(getProcessors(2), 1)
    expect_equal(getProcessors(NULL), 1)
})


