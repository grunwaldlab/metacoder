context("taxa")

x <- taxon(
  name = taxon_name("Poa annua"),
  rank = taxon_rank("species"),
  id = taxon_id(93036)
)

test_that("taxa works", {
  aa <- taxa(x, x, x)

  expect_is(aa, "taxa")
  expect_is(unlist(aa), "list")
  expect_is(aa[[1]], "Taxon")
  expect_is(aa[[2]], "Taxon")
  expect_is(aa[[3]], "Taxon")
})

test_that("taxa - empty", {
  aa <- taxa()

  expect_is(aa, "taxa")
  expect_is(unclass(aa), "list")
  expect_equal(length(aa), 0)
})

test_that("taxa - print method", {
  # no inputs
  expect_output(print(taxa()), "<taxa>")
  expect_output(print(taxa()), "no\\. taxa:  0")
  expect_output(print(taxa(x)), "no\\. taxa:  1")
  expect_output(print(taxa(x)), "Poa annua / species / 93036")
})

test_that("taxa fails well", {
  expect_error(taxa(5), "all inputs to 'taxa' must be of class 'Taxon'")
  expect_error(taxa(mtcars), "all inputs to 'taxa' must be of class 'Taxon'")
  expect_error(taxa(4, x, "adff"),
               "all inputs to 'taxa' must be of class 'Taxon'")
})

test_that("dots and .list return the same output", {
  expect_equal(taxa(x, x, x), taxa(.list = list(x, x, x)))
  expect_error(taxa(taxon("a"), .list = list(taxon("a"))),
               'Both `...` and `.list` were supplied')
})

