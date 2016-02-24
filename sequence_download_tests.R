context("Sequence download")

test_that("query_taxon returns the right type", {
  seqinr::choosebank("genbank")
  expect_true(class(query_taxon("Phytophthora ramorum")) == "qaw")
  expect_true(class(query_taxon("Phytophthora ramorum", execute = F)) == "character")
  expect_true(length(query_taxon("Phytophthora ramorum", execute = F)) == 1)
  expect_true(is.null(query_taxon(NULL)))
  expect_true(is.null(query_taxon(character(0))))
  expect_true(is.na(query_taxon(NA)))
  seqinr::closebank()
})