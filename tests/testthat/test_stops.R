context("count stop codons function")
data(drosophila)

s = DNAString(drosophila[[1]]@data$PBAS.2)
test_that("input variables", {
  expect_warning(count.stop.codons(s[1:2]), "Cannot calculate")
  expect_error(count.stop.codons("Not a DNAString"), "sequence must")
  expect_error(count.stop.codons(s, reading.frame = 10))
})

s2 = DNAString("TAATTAATAATTAATAATAA")
test_that("reading.frames work", {
    expect_equal(count.stop.codons(s2, reading.frame = 1), 1)
    expect_equal(count.stop.codons(s2, reading.frame = 2), 2)
    expect_equal(count.stop.codons(s2, reading.frame = 3), 3)
})

test_that("we spot bogus genetic codes", {
    expect_error(count.stop.codons(s, reading.frame = 1, genetic.code = GENETIC_CODE[1]), "Your genetic")
})