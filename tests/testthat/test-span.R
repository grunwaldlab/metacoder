context("span")

test_that("span: ranks", {
  # without piping
  aa1 <- span(ex_hierarchy1, ranks("family", "genus"))
  expect_is(aa1, "Hierarchy")
  expect_equal(names(aa1$ranklist), c('family', 'genus'))
  expect_equal(length(aa1$taxa), 2)

  # with piping
  aa2 <- ex_hierarchy1 %>% span(ranks("family", "genus"))
  expect_is(aa2, "Hierarchy")
  expect_equal(names(aa2$ranklist), c('family', 'genus'))
  expect_equal(length(aa2$taxa), 2)
  expect_equal(aa1, aa2)

  # with piping, comma separated entries
  aa3 <- ex_hierarchy1 %>% span(ranks("family", "genus"))
  expect_is(aa3, "Hierarchy")
  expect_equal(names(aa3$ranklist), c('family', 'genus'))
  expect_equal(length(aa3$taxa), 2)
})

test_that("span: names", {
  # without piping
  aa1 <- span(ex_hierarchy1, nms("Poaceae", "Poa"))
  expect_is(aa1, "Hierarchy")
  expect_equal(names(aa1$ranklist), c('family', 'genus'))
  expect_equal(length(aa1$taxa), 2)

  # with piping
  aa2 <- ex_hierarchy1 %>% span(nms("Poaceae", "Poa"))
  expect_is(aa2, "Hierarchy")
  expect_equal(names(aa2$ranklist), c('family', 'genus'))
  expect_equal(length(aa2$taxa), 2)
  expect_equal(aa1, aa2)

  # with piping, many entries
  aa3 <- ex_hierarchy1 %>% span(nms("Poaceae", "Poa"))
  expect_is(aa3, "Hierarchy")
  expect_equal(names(aa3$ranklist), c('family', 'genus'))
  expect_equal(length(aa3$taxa), 2)
})

test_that("span: ids", {
  # without piping
  aa1 <- span(ex_hierarchy1, ids(4479, 4544))
  expect_is(aa1, "Hierarchy")
  expect_equal(names(aa1$ranklist), c('family', 'genus'))
  expect_equal(length(aa1$taxa), 2)

  # with piping
  aa2 <- ex_hierarchy1 %>% span(ids(4479, 4544))
  expect_is(aa2, "Hierarchy")
  expect_equal(names(aa2$ranklist), c('family', 'genus'))
  expect_equal(length(aa2$taxa), 2)
  expect_equal(aa1, aa2)

  # with piping, comma separated entries
  aa3 <- ex_hierarchy1 %>% span(ids(4479, 4544))
  expect_is(aa3, "Hierarchy")
  expect_equal(names(aa3$ranklist), c('family', 'genus'))
  expect_equal(length(aa3$taxa), 2)
})


test_that("no variables given to pick - errors", {
  expect_error(span(ex_hierarchy1), "no acceptable selectors passed in")
})

test_that("span fails well", {
  expect_error(span(),
               "argument \".data\" is missing")
  expect_error(span(5),
               "no 'span' method for numeric")
  expect_error(span("adf", 5),
               "no 'span' method for character")
  expect_error(span(ex_hierarchy1, ranks("family")),
               "if no operator, must pass in 2 names")
})

# FIXME: make bigger hierarchy first
# test_that("span: mixed ranks, names and ids", {
#   aa1 <- ex_hierarchy1 %>% span(ranks(family, species), ids(4544))
#
#   expect_is(aa1, "Hierarchy")
#   expect_equal(names(aa1$ranklist), 'species')
#   expect_equal(length(aa1$taxa), 1)
# })
