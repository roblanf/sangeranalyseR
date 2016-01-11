context("trim.mott function")
data(drosophila)


test_that("input variables", {
  expect_error(trim.mott("this is not abif"), "abif.seq must be an")
  expect_error(trim.mott("", "egg"), "cutoff must be a")
  expect_error(trim.mott("", "-1"), "cutoff must be a")
})

test_that("trimming works", {
    r1 = trim.mott(drosophila[[1]], 1000) # shouldn't trim anything
    expect_equal(r1$start, 1)
    expect_equal(r1$finish, length(drosophila[[1]]@data$PCON.2))

    r2 = trim.mott(drosophila[[1]]) # default trims a bit
    expect_equal(r2$start, 27)
    expect_equal(r2$finish, 477)
})