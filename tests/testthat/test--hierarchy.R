library(testthat)
context("hierarchy")

## Creating test data
plantae <- taxon(
  name = taxon_name("Plantae")
)
solanaceae <- taxon(
  name = taxon_name("Solanaceae")
)
solanum <- taxon(
  name = taxon_name("Solanum")
)
sl <- taxon(
  name = taxon_name("Solanum lycopersicum")
)


test_that("Characters as inputs", {
  aa <- hierarchy(plantae, solanaceae, solanum, sl)

  expect_is(aa, "Hierarchy")
  expect_is(aa$taxa, "list")
  expect_is(aa$taxa[[1]], "Taxon")
  expect_is(aa$print, "function")
  expect_equal(
    aa,
    hierarchy("Plantae", "Solanaceae", "Solanum", "Solanum lycopersicum")
  )
})


test_that("hierarchy - empty", {
  aa <- hierarchy()

  expect_is(aa, "Hierarchy")
  expect_null(aa$taxa)
  expect_null(aa$ranklist)

  # prints 'Empty hierarchy'
  expect_output(
    print(hierarchy()),
    "Empty hierarchy"
  )
  expect_output(
    print(hierarchy()),
    "<Hierarchy>"
  )
})

plantae <- taxon(
  name = taxon_name("Plantae"),
  rank = "kingdom"
)
solanaceae <- taxon(
  name = taxon_name("Solanaceae"),
  rank = "family"
)
solanum <- taxon(
  name = taxon_name("Solanum"),
  rank = "genus"
)
sl <- taxon(
  name = taxon_name("Solanum lycopersicum"),
  rank = "species"
)

test_that("hierarchy - print when not empty", {
  expect_output(
    print(hierarchy(plantae, solanaceae, solanum, sl)),
    "Plantae / kingdom /"
  )
  expect_output(
    print(hierarchy(plantae, solanaceae, solanum, sl)),
    "Solanaceae / family /"
  )
  expect_output(
    print(hierarchy(plantae, solanaceae, solanum, sl)),
    "Solanum / genus /"
  )
  expect_output(
    print(hierarchy(plantae, solanaceae, solanum, sl)),
    "Solanum lycopersicum / species /"
  )
})


test_that("hierarchy fails well", {
  expect_error(
    hierarchy(4),
    "all inputs to 'hierarchy' must be of class 'Taxon' or 'character'")
  expect_error(
    hierarchy(solanum, 5),
    "all inputs to 'hierarchy' must be of class 'Taxon' or 'character'")
})


test_that("dots and .list return the same output", {
  expect_equal(hierarchy(plantae, solanaceae, solanum, sl),
               hierarchy(.list = list(plantae, solanaceae, solanum, sl)))
})

test_that("hierarchy can do null data", {
  # empty hierarchy()
  x <- hierarchy()
  expect_is(x, "Hierarchy")
  expect_null(x$taxa)
  expect_null(x$ranklist)

  # specifying NULL
  x <- hierarchy(NULL)
  expect_is(x, "Hierarchy")
  expect_null(x$taxa)
  expect_null(x$ranklist)
})



test_that("hierarchy pop", {
  x <- taxon(
    name = taxon_name("Poaceae"),
    rank = taxon_rank("family"),
    id = taxon_id(4479)
  )

  y <- taxon(
    name = taxon_name("Poa"),
    rank = taxon_rank("genus"),
    id = taxon_id(4544)
  )

  z <- taxon(
    name = taxon_name("Poa annua"),
    rank = taxon_rank("species"),
    id = taxon_id(93036)
  )

  obj <- hierarchy(z, y, x)

  res <- pop(obj, ranks("species"))
  expect_equal(length(res$taxa), 2)
  expect_equal(res, pop(obj, nms("Poa annua")))
  expect_equal(res, pop(obj, ids("93036")))
  expect_error(pop(hierarchy(), nms("Poa annua")),
               "no taxa found")
})


test_that("hierarchy pick", {
  x <- taxon(
    name = taxon_name("Poaceae"),
    rank = taxon_rank("family"),
    id = taxon_id(4479)
  )

  y <- taxon(
    name = taxon_name("Poa"),
    rank = taxon_rank("genus"),
    id = taxon_id(4544)
  )

  z <- taxon(
    name = taxon_name("Poa annua"),
    rank = taxon_rank("species"),
    id = taxon_id(93036)
  )

  obj <- hierarchy(z, y, x)

  res <- pick(obj, ranks("species"))
  expect_equal(length(res$taxa), 1)
  expect_equal(res, pick(obj, nms("Poa annua")))
  expect_equal(res, pick(obj, ids("93036")))
  expect_error(pick(hierarchy(), nms("Poa annua")),
               "no taxa found")
})


test_that("hierarchy span", {
  x <- taxon(
    name = taxon_name("Poaceae"),
    rank = taxon_rank("family"),
    id = taxon_id(4479)
  )

  y <- taxon(
    name = taxon_name("Poa"),
    rank = taxon_rank("genus"),
    id = taxon_id(4544)
  )

  z <- taxon(
    name = taxon_name("Poa annua"),
    rank = taxon_rank("species"),
    id = taxon_id(93036)
  )

  obj <- hierarchy(z, y, x)

  res <- span(obj, ranks("family", "genus"))
  expect_equal(length(res$taxa), 2)
  expect_equal(res, span(obj, nms("Poaceae", "Poa")))
  expect_equal(res, span(obj, ids("4544", "4479")))
  expect_error(span(hierarchy(), nms("Poa annua")),
               "no taxa found")
})


