context("hierarchies")

## Creating test data
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
hier1 <- hierarchy(z, y, x)

a <- taxon(
  name = taxon_name("Felidae"),
  rank = taxon_rank("family"),
  id = taxon_id(9681)
)
b <- taxon(
  name = taxon_name("Puma"),
  rank = taxon_rank("genus"),
  id = taxon_id(146712)
)
c <- taxon(
  name = taxon_name("Puma concolor"),
  rank = taxon_rank("species"),
  id = taxon_id(9696)
)
hier2 <- hierarchy(c, b, a)


test_that("hierarchies works", {
  aa <- hierarchies(hier1, hier2)

  expect_is(aa, "hierarchies")
  expect_is(unclass(aa), "list")
  expect_is(aa[[1]], "Hierarchy")
  expect_is(aa[[2]], "Hierarchy")

  expect_equal(length(aa), 2)
})


test_that("hierarchies - empty", {
  aa <- hierarchies(hierarchy(), hierarchy())

  expect_is(aa, "hierarchies")
  expect_is(unclass(aa), "list")
  expect_is(aa[[1]], "Hierarchy")
  expect_is(aa[[2]], "Hierarchy")
  expect_null(aa[[1]]$taxa)
  expect_null(aa[[2]]$taxa)

  # prints 'Empty hierarchy'
  expect_output(
    print(hierarchies(hierarchy())),
    "Empty hierarchy"
  )
  expect_output(
    print(hierarchies(hierarchy(), hierarchy())),
    "Empty hierarchy\n  Empty hierarchy"
  )

  aa <- hierarchies()

  expect_match(paste0(capture.output(hierarchies()), collapse = ""),
               "<Hierarchies>")
  expect_is(aa, "hierarchies")
  expect_is(unclass(aa), "list")
  expect_equal(length(aa), 0)
})


test_that("hierarchies - print when not empty", {
  expect_output(
    print(hierarchies(hier1)),
    "Poaceae / Poa / Poa annua"
  )

  expect_output(
    print(hierarchies(hier1, hier2)),
    "Poaceae / Poa / Poa annua \n  Felidae / Puma / Puma concolor"
  )
})


test_that("hierarchies fails well", {
  expect_error(hierarchies(4),
               "all inputs to 'hierarchies' must be of class 'Hierarchy'")
  expect_error(hierarchies("a", "b", "c"),
               "all inputs to 'hierarchies' must be of class 'Hierarchy'")
  expect_error(hierarchies(hier1, "c"),
               "all inputs to 'hierarchies' must be of class 'Hierarchy'")
})


test_that("dots and .list return the same output", {
  expect_equal(hierarchies(hier1, hier2),
               hierarchies(.list = list(hier1, hier2)))
})
