context("taxon_id")

test_that("taxa_id - just ID", {
  aa <- taxon_id(93036)

  expect_is(aa, "TaxonId")
  expect_null(aa$database)
  expect_type(aa$id, "double")
  expect_output(aa$print())
})

test_that("taxa_id - ID and database (TaxonDatabase)", {
  aa <- taxon_id(93036, database_list$ncbi)

  expect_is(aa, "TaxonId")
  expect_is(aa$database, "TaxonDatabase")
  expect_equal(aa$database$name, "ncbi")
  expect_equal(aa$database$url, "http://www.ncbi.nlm.nih.gov/taxonomy")
  expect_type(aa$id, "double")
  expect_output(aa$print())
})

test_that("taxa_id - ID and database (character)", {
  aa <- taxon_id(93036, "ncbi")

  expect_is(aa, "TaxonId")
  expect_is(aa$database, "TaxonDatabase")
  expect_equal(aa$database$name, "ncbi")
  expect_null(aa$database$url)
  expect_type(aa$id, "double")
  expect_output(aa$print())
})

test_that("taxon_id fails well", {
  expect_error(taxon_id(),
               "argument \"id\" is missing")
  expect_error(taxon_id(mtcars),
               "id must be of class character, integer, numeric")
  expect_error(taxon_id("adf", 5),
               "database must be of class character")
})

test_that("taxon_id can do null data", {
  # empty taxon_id() tested in above block
  x <- taxon_id(NULL)
  expect_is(x, "TaxonId")
  expect_null(x$name)
  expect_null(x$database)
})
