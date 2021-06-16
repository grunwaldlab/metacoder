context("taxon")

name <- taxon_name("Poa annua")
rank <- taxon_rank("species")
id <- taxon_id(93036)

test_that("taxa works", {
  aa <- taxon(name, rank, id)

  expect_is(aa, "Taxon")
  expect_null(aa$authority)
  expect_is(aa$id, "TaxonId")
  expect_is(aa$name, "TaxonName")
  expect_is(aa$rank, "TaxonRank")
  expect_output(aa$print())
})

test_that("taxa fails well", {
  expect_error(taxon(rank = "adfd"),
               "argument \"name\" is missing")
  expect_error(taxon(5),
               "name must be of class TaxonName, character")
  expect_error(taxon("adf", 5),
               "rank must be of class TaxonRank, character")
  expect_error(taxon("adfadsf", authority = 23),
               "authority must be of class character")
})

test_that("taxon can do null data", {
  x <- taxon(NULL)
  expect_is(x, "Taxon")
  expect_null(x$name)
  expect_null(x$id)
})
