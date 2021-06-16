## Testing `taxonomy` class

library(taxa)
library(testthat)
context("taxonomy")


## Creating test data
notoryctidae <- taxon(
  name = taxon_name("Notoryctidae"),
  rank = taxon_rank("family"),
  id = taxon_id(4479)
)
notoryctes <- taxon(
  name = taxon_name("Notoryctes"),
  rank = taxon_rank("genus"),
  id = taxon_id(4544)
)
typhlops <- taxon(
  name = taxon_name("typhlops"),
  rank = taxon_rank("species"),
  id = taxon_id(93036)
)
mammalia <- taxon(
  name = taxon_name("Mammalia"),
  rank = taxon_rank("class"),
  id = taxon_id(9681)
)
felidae <- taxon(
  name = taxon_name("Felidae"),
  rank = taxon_rank("family"),
  id = taxon_id(9681)
)
puma <- taxon(
  name = taxon_name("Puma"),
  rank = taxon_rank("genus"),
  id = taxon_id(146712)
)
concolor <- taxon(
  name = taxon_name("concolor"),
  rank = taxon_rank("species"),
  id = taxon_id(9696)
)
panthera <- taxon(
  name = taxon_name("Panthera"),
  rank = taxon_rank("genus"),
  id = taxon_id(146712)
)
tigris <- taxon(
  name = taxon_name("tigris"),
  rank = taxon_rank("species"),
  id = taxon_id(9696)
)
plantae <- taxon(
  name = taxon_name("Plantae"),
  rank = taxon_rank("kingdom"),
  id = taxon_id(33090)
)
solanaceae <- taxon(
  name = taxon_name("Solanaceae"),
  rank = taxon_rank("family"),
  id = taxon_id(4070)
)
solanum <- taxon(
  name = taxon_name("Solanum"),
  rank = taxon_rank("genus"),
  id = taxon_id(4107)
)
lycopersicum <- taxon(
  name = taxon_name("lycopersicum"),
  rank = taxon_rank("species"),
  id = taxon_id(49274)
)
tuberosum <- taxon(
  name = taxon_name("tuberosum"),
  rank = taxon_rank("species"),
  id = taxon_id(4113)
)
unidentified <- taxon(name = taxon_name("unidentified"))

tiger <- hierarchy(mammalia, felidae, panthera, tigris)
cougar <- hierarchy(mammalia, felidae, puma, concolor)
mole <- hierarchy(mammalia, notoryctidae, notoryctes, typhlops)
tomato <- hierarchy(plantae, solanaceae, solanum, lycopersicum)
potato <- hierarchy(plantae, solanaceae, solanum, tuberosum)
potato_partial <- hierarchy(solanaceae, solanum, tuberosum)
unidentified_animal <- hierarchy(mammalia, unidentified)
unidentified_plant <- hierarchy(plantae, unidentified)

test_that("NSE", {
  x <- taxonomy(tiger, cougar, mole)
  expect_equal(all_names(x), x$all_names())
})

test_that("Printing taxonomy", {
  x <- taxonomy(tiger, cougar, mole)
  expect_output(print(x), "Taxonomy")
  expect_output(print(x), "9 taxa")
})

test_that("Simple usage", {
  x <- taxonomy(tiger, cougar, mole)
  expect_length(x$taxa, 9)
  expect_equal(dim(x$edge_list), c(9, 2))
  expect_length(x$roots(), 1)
})


test_that("Multiple roots", {
  x <- taxonomy(tiger, cougar, mole, tomato, potato)
  expect_length(x$taxa, 14)
  expect_equal(dim(x$edge_list), c(14, 2))
  expect_length(x$roots(), 2)
})


test_that("Hierarchies of different lengths", {
  x <- taxonomy(tiger, unidentified_animal)
  expect_length(x$taxa, 5)
  expect_equal(dim(x$edge_list), c(5, 2))
  expect_length(x$roots(), 1)
})


test_that("Same taxon name, different lineage", {
  x <- taxonomy(unidentified_plant, unidentified_animal)
  expect_length(x$taxa, 4)
  expect_equal(dim(x$edge_list), c(4, 2))
  expect_length(x$roots(), 2)
  expect_equal(sum(sapply(x$taxa, function(x) x$name$name) == "unidentified"), 2)
})


test_that("Edge cases", {
  x <- taxonomy()
  expect_length(x$taxa, 0)
  expect_equal(dim(x$edge_list), c(0, 2))
  expect_is(taxonomy(hierarchy()), "Taxonomy")
  expect_equal(length(taxonomy(hierarchy())$taxa), 0)
  expect_length(x$taxa, 0)
  expect_equal(dim(x$edge_list), c(0, 2))
})


test_that("Characters as inputs", {
  x <- taxonomy(c("a", "b", "c"), c("a", "d"))
  expect_length(x$taxa, 4)
  expect_equal(dim(x$edge_list), c(4, 2))
  expect_length(x$roots(), 1)

  # x <- taxonomy(list(c("a", "b", "c"), c("a", "d"))) # does not work yet
})

test_that("Accessing basic info", {
  x <- taxonomy(tiger, cougar, mole, tomato, potato,
                unidentified_plant, unidentified_animal)
  expect_equal(x$taxon_indexes(), taxon_indexes(x))
  expect_equivalent(taxon_indexes(x), seq_along(x$taxa))
})


test_that("Finding roots", {
  x <- taxonomy(tiger, cougar, mole, tomato, potato,
                unidentified_plant, unidentified_animal)
  expect_equal(x$roots(), roots(x))

  # Index return type
  expect_type(roots(x, value = "taxon_indexes"), "integer")

  # Taxon ID return type
  expect_type(roots(x, value = "taxon_ids"), "character")
})


test_that("Finding internodes", {
  x <- taxonomy(tiger, cougar, mole, tomato, potato,
                unidentified_plant, unidentified_animal)
  expect_equal(x$internodes(), internodes(x))
  expect_equal(x$is_internode(), is_internode(x))

  # Index return type
  expect_type(internodes(x, value = "taxon_indexes"), "integer")

  # Taxon ID return type
  expect_type(internodes(x, value = "taxon_ids"), "character")
})


test_that("Finding id_classifications", {
  x <- taxonomy(tiger, cougar, mole, tomato, potato,
                unidentified_plant, unidentified_animal)
  expect_equal(x$id_classifications(), id_classifications(x))
})


test_that("Finding id_classifications", {
  x <- taxonomy(tiger, cougar, mole, tomato, potato,
                unidentified_plant, unidentified_animal)
  expect_equal(x$classifications(), classifications(x))
})


test_that("Finding n_supertaxa", {
  x <- taxonomy(tiger, cougar, mole, tomato, potato,
                unidentified_plant, unidentified_animal)
  expect_equal(x$n_supertaxa(), n_supertaxa(x))
})


test_that("Finding n_supertaxa_1", {
  x <- taxonomy(tiger, cougar, mole, tomato, potato,
                unidentified_plant, unidentified_animal)
  expect_equal(x$n_supertaxa_1(), n_supertaxa_1(x))
})


test_that("Finding n_subtaxa", {
  x <- taxonomy(tiger, cougar, mole, tomato, potato,
                unidentified_plant, unidentified_animal)
  expect_equal(x$n_subtaxa(), n_subtaxa(x))
})


test_that("Finding n_subtaxa_1", {
  x <- taxonomy(tiger, cougar, mole, tomato, potato,
                unidentified_plant, unidentified_animal)
  expect_equal(x$n_subtaxa_1(), n_subtaxa_1(x))
})


test_that("Finding branches", {
  x <- taxonomy(tiger, cougar, mole, tomato, potato,
                unidentified_plant, unidentified_animal)
  expect_equal(x$branches(), branches(x))
  expect_equal(x$is_branch(), is_branch(x))

  # Index return type
  expect_type(branches(x, value = "taxon_indexes"), "integer")

  # Taxon ID return type
  expect_type(branches(x, value = "taxon_ids"), "character")

  # Expected output
  expect_equal(which(! is_root(x) & ! is_leaf(x)), branches(x))

})


test_that("Finding supertaxa", {
  x <- taxonomy(tiger, cougar, mole, tomato, potato,
                unidentified_plant, unidentified_animal)
  expect_equal(x$supertaxa(), supertaxa(x))

  # Index return type
  expect_type(supertaxa(x, value = "taxon_indexes")[[1]], "integer")
  expect_type(supertaxa(x, value = "taxon_indexes", simplify = TRUE), "integer")

  # Taxon ID return type
  expect_type(supertaxa(x, value = "taxon_ids")[[1]], "character")
  expect_type(supertaxa(x, value = "taxon_ids", simplify = TRUE), "character")

  # Recursion settings
  expect_equal(supertaxa(x, recursive = TRUE), supertaxa(x, recursive = -1))
  expect_equal(supertaxa(x, recursive = FALSE), supertaxa(x, recursive = 1))
  expect_equal(max(vapply(supertaxa(x, recursive = 2), length, numeric(1))), 2)

  # Duplicated inputs
  expect_equal(names(x$supertaxa(c(1, 2, 1, 1))), c("b", "c", "b", "b"))
})


test_that("Finding subtaxa", {
  x <- taxonomy(tiger, cougar, mole, tomato, potato,
                unidentified_plant, unidentified_animal)
  expect_equal(x$subtaxa(), subtaxa(x))

  # Index return type
  expect_type(subtaxa(x, value = "taxon_indexes")[[1]], "integer")
  expect_type(subtaxa(x, value = "taxon_indexes", simplify = TRUE), "integer")

  # Taxon ID return type
  expect_type(subtaxa(x, value = "taxon_ids")[[1]], "character")
  expect_type(subtaxa(x, value = "taxon_ids", simplify = TRUE), "character")

  # Subsets and NSE
  my_var <- 2
  expect_equivalent(eval(substitute(sapply(subtaxa(x, subset = n_subtaxa == my_var), length))), c(2, 2))

  # Recursion settings
  expect_equal(subtaxa(x, recursive = TRUE), subtaxa(x, recursive = -1))
  expect_equal(subtaxa(x, recursive = FALSE), subtaxa(x, recursive = 1))
  expect_equivalent(names(subtaxa(x, subset = "e", recursive = 2)$e), c("k", "o"))

  # Edge cases
  expect_equal(subtaxa(x, subset = rep(FALSE, 16)), list())
  expect_equal(subtaxa(x, subset = rep(FALSE, 16), simplify = TRUE),
               integer(0))
})


test_that("Finding stems", {
  x <- taxonomy(tiger, cougar, mole, tomato, potato,
                unidentified_plant, unidentified_animal)
  expect_equal(x$stems(), stems(x))
  expect_equal(x$is_stem(), is_stem(x))

  # Index return type
  expect_type(stems(x, value = "taxon_indexes")[[1]], "integer")
  expect_type(stems(x, value = "taxon_indexes", simplify = TRUE), "integer")

  # Taxon ID return type
  expect_type(stems(x, value = "taxon_ids")[[1]], "character")
  expect_type(stems(x, value = "taxon_ids", simplify = TRUE), "character")
})


test_that("Finding leaves", {
  x <- taxonomy(tiger, cougar, mole, tomato, potato,
                unidentified_plant, unidentified_animal)
  expect_equal(x$leaves(), leaves(x))
  expect_equal(x$is_leaf(), is_leaf(x))

  # Index return type
  expect_type(leaves(x, value = "taxon_indexes")[[1]], "integer")

  # Taxon ID return type
  expect_type(leaves(x, value = "taxon_ids")[[1]], "character")

  # leaves_apply
  expect_equal(sum(leaves_apply(x, length, subset = c(1, 2), simplify = TRUE)), 7)

  # n_leaves
  expect_equal(n_leaves(x), unlist(leaves_apply(x, length)))

  # n_leaves_1
  expect_equivalent(n_leaves_1(x)["l"], 2)
})

test_that("Filtering taxa", {
  x <- taxonomy(tiger, cougar, mole, tomato, potato,
                unidentified_plant, unidentified_animal)

  result <- filter_taxa(x, taxon_names == "Solanum")
  expect_equal(result$taxon_names(), c("l" = "Solanum"))
  expect_warning(filter_taxa(x, taxon_names == "Solanum", drop_obs = FALSE))
  expect_warning(filter_taxa(x, taxon_names == "Solanum", reassign_obs = TRUE))

  # Check that filtering does not change order of taxa
  result <- filter_taxa(x, taxon_names != "tuberosum")
  expected_names <- taxon_names(x)
  expected_names <- expected_names[expected_names != "tuberosum"]
  expect_true(all(expected_names == taxon_names(result)))

  result <- filter_taxa(x, taxon_names == "Solanum", subtaxa = TRUE, invert = TRUE)
  expected_names <- taxon_names(x)
  expected_names <- expected_names[! expected_names %in% c("Solanum", "lycopersicum", "tuberosum")]
  expect_true(all(expected_names == taxon_names(result)))

  # Errors for invalid indexes
  expect_error(filter_taxa(x, 100), "The following taxon indexes are invalid:")

  # Errors for invalid IDs
  expect_error(filter_taxa(x, "zzz"), "The following taxon IDs do not exist:")

  # Errors for invalid logical
  expect_error(filter_taxa(x, TRUE), "must be the same length as the number of taxa")

  # Edge case: filtering everything out
  result <- filter_taxa(x, numeric(0))
  expect_equal(length(result$taxa), 0)
  expect_equal(result, filter_taxa(x, NULL, numeric(0)))
  expect_equal(result, filter_taxa(x, "c", numeric(0)))

  # Edge case: NULL input (shou)
  expect_equal(filter_taxa(x, NULL), x)
  expect_equal(filter_taxa(x, NULL, NULL), x)
  expect_equal(filter_taxa(x, NULL, "c"), filter_taxa(x, "c"))
})


test_that("Sampling taxa",  {
  x <- taxonomy(tiger, cougar, mole, tomato, potato,
                unidentified_plant, unidentified_animal)

  result <- sample_n_taxa(x, size = 3)
  expect_equal(length(taxon_ids(result)), 3)
  expect_warning(sample_n_taxa(x, size = 3, obs_weight = 1))
  expect_warning(sample_n_taxa(x, size = 3, obs_target = 1))

  result <- sample_frac_taxa(x, size = 0.5)
  expect_equal(length(taxon_ids(result)), 8)
})

test_that("Mapping vairables",  {
  x <- taxonomy(tiger, cougar, mole, tomato, potato,
                unidentified_plant, unidentified_animal)
  result <- map_data(x, taxon_names, taxon_ranks)
  expect_equal(result, map_data_(x, taxon_names(x), taxon_ranks(x)))
  expect_equivalent(result, taxon_ranks(x))
  expect_equivalent(names(result), taxon_names(x))
  expect_warning(map_data(x, taxon_names, c("e" = "A", "e" = "B")))
  expect_silent(map_data(x, taxon_names, c("e" = "A", "e" = "B"), warn = FALSE))
  expect_error(map_data(x, taxon_names, 1:10))
})


test_that("dots and .list return the same output", {
  expect_equal(taxonomy(tiger, cougar, mole, tomato, potato),
               taxonomy(.list = list(tiger, cougar, mole, tomato, potato)))
})

test_that("get data frame", {
  x <- taxonomy(tiger, cougar, mole, tomato, potato,
                unidentified_plant, unidentified_animal)
  expect_equal(x$get_data_frame(), get_data_frame(x))

  df <- x$get_data_frame()
  expect_is(df, "data.frame")
  expect_is(df$taxon_names, "character")

  # select columns to return
  expect_named(x$get_data_frame("taxon_ids"), "taxon_ids")
})


test_that("supertaxa_apply function", {
  x <- taxonomy(tiger, cougar, mole, tomato, potato,
                unidentified_plant, unidentified_animal)
  expect_equal(lapply(supertaxa(x), length),
               supertaxa_apply(x, length))
})


test_that("subtaxa_apply function", {
  x <- taxonomy(tiger, cougar, mole, tomato, potato,
                unidentified_plant, unidentified_animal)
  expect_equal(lapply(subtaxa(x), length),
               subtaxa_apply(x, length))
})


test_that("replacing taxon IDs", {
  x <- taxonomy(tiger, cougar, mole, tomato, potato,
                unidentified_plant, unidentified_animal)
  result <- replace_taxon_ids(x, 1:16)
  expect_equivalent(taxon_ids(result), 1:16)
  expect_error(replace_taxon_ids(x, 1:10), "different than the current number of taxa")
  expect_error(replace_taxon_ids(x, rep(1, 16)), "New taxon IDs must be unique")
})


test_that("removing redundant names", {
  lycopersicum <- taxon(
    name = taxon_name("Solanum lycopersicum"),
    rank = taxon_rank("species"),
    id = taxon_id(49274)
  )
  tuberosum <- taxon(
    name = taxon_name("Solanum tuberosum"),
    rank = taxon_rank("species"),
    id = taxon_id(4113)
  )
  x <- taxonomy(tiger, cougar, mole, tomato, potato,
                unidentified_plant, unidentified_animal)
  result <- remove_redundant_names(x)
  expect_true(all(c("tuberosum", "lycopersicum") %in% taxon_names(result)))
})


test_that("taxonomy can be converted to tables", {
  x <- taxonomy(tiger, cougar, mole, tomato, potato,
                unidentified_plant, unidentified_animal)
  expect_message(result <- taxonomy_table(x),
                 "The following ranks will not be included")
  expect_equal(colnames(result), c("class", "family", "genus", "species"))
  result <- taxonomy_table(x, use_ranks = FALSE)
  expect_equal(colnames(result), c("rank_1", "rank_2", "rank_3", "rank_4"))
})


test_that("print_tree works", {
  x <- taxonomy(tiger, cougar, mole, tomato, potato,
                unidentified_plant, unidentified_animal)
  expect_equal(print_tree(x)[1], "Mammalia")
})

