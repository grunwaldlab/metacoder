context("taxon_rank")

test_that("taxa_rank - just rank", {
  aa <- taxon_rank("species")

  expect_is(aa, "TaxonRank")
  expect_null(aa$database)
  expect_type(aa$name, "character")
  expect_output(print(aa$print()), "<TaxonRank>")
})

test_that("taxa_rank - name and database (TaxonRank)", {
  aa <- taxon_rank("species", database_list$ncbi)

  expect_is(aa, "TaxonRank")
  expect_is(aa$database, "TaxonDatabase")
  expect_equal(aa$database$name, "ncbi")
  expect_equal(aa$database$url, "http://www.ncbi.nlm.nih.gov/taxonomy")
  expect_type(aa$name, "character")
})

test_that("taxa_rank - name and database (character)", {
  aa <- taxon_rank("genus", "ncbi")

  expect_is(aa, "TaxonRank")
  expect_is(aa$database, "TaxonDatabase")
  expect_equal(aa$database$name, "ncbi")
  expect_null(aa$database$url)
  expect_type(aa$name, "character")
})

test_that("taxon_rank fails well", {
  expect_error(taxon_rank(),
               "argument \"name\" is missing")
  expect_error(taxon_rank(mtcars),
               "name must be of class character, TaxonName")
  expect_error(taxon_rank("adf", 5),
               "database must be of class character, TaxonDatabase")
})

test_that("taxon_rank can do null data", {
  # empty taxon_rank() tested in above block
  x <- taxon_rank(NULL)
  expect_is(x, "TaxonRank")
  expect_null(x$name)
  expect_null(x$database)
})
