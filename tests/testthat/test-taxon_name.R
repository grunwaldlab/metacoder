context("taxon_nane")

test_that("taxa_name - just", {
  aa <- taxon_name("Poa")

  expect_is(aa, "TaxonName")
  expect_null(aa$database)
  expect_type(aa$name, "character")
  expect_output(aa$print())
})

test_that("taxa_name - name and database (TaxonDatabase)", {
  aa <- taxon_name("Poa", database_list$ncbi)

  expect_is(aa, "TaxonName")
  expect_is(aa$database, "TaxonDatabase")
  expect_equal(aa$database$name, "ncbi")
  expect_equal(aa$database$url, "http://www.ncbi.nlm.nih.gov/taxonomy")
  expect_type(aa$name, "character")
  expect_output(aa$print())
})

test_that("taxa_name - ID and database (character)", {
  aa <- taxon_name("Poa", "ncbi")

  expect_is(aa, "TaxonName")
  expect_is(aa$database, "TaxonDatabase")
  expect_equal(aa$database$name, "ncbi")
  expect_null(aa$database$url)
  expect_type(aa$name, "character")
  expect_output(aa$print())
})

test_that("taxon_name fails well", {
  expect_error(taxon_name(),
               "argument \"name\" is missing")
  expect_error(taxon_name(mtcars),
               "name must be of class character")
  expect_error(taxon_name("adf", 5),
               "database must be of class character, TaxonDatabase")
})

test_that("taxon_name can do null data", {
  # empty taxon_name() tested in above block
  x <- taxon_name(NULL)
  expect_is(x, "TaxonName")
  expect_null(x$name)
  expect_null(x$database)
})
