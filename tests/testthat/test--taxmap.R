## Testing `taxmap` class

library(taxa)
library(testthat)
context("taxmap")



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
felis <- taxon(
  name = taxon_name("Felis"),
  rank = taxon_rank("genus"),
  id = taxon_id(9682)
)
catus <- taxon(
  name = taxon_name("catus"),
  rank = taxon_rank("species"),
  id = taxon_id(9685)
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
homo <- taxon(
  name = taxon_name("homo"),
  rank = taxon_rank("genus"),
  id = taxon_id(9605)
)
sapiens <- taxon(
  name = taxon_name("sapiens"),
  rank = taxon_rank("species"),
  id = taxon_id(9606)
)
hominidae <- taxon(
  name = taxon_name("Hominidae"),
  rank = taxon_rank("family"),
  id = taxon_id(9604)
)
unidentified <- taxon(
  name = taxon_name("unidentified")
)

tiger <- hierarchy(mammalia, felidae, panthera, tigris)
cat <- hierarchy(mammalia, felidae, felis, catus)
human <- hierarchy(mammalia, hominidae, homo, sapiens)
mole <- hierarchy(mammalia, notoryctidae, notoryctes, typhlops)
tomato <- hierarchy(plantae, solanaceae, solanum, lycopersicum)
potato <- hierarchy(plantae, solanaceae, solanum, tuberosum)
potato_partial <- hierarchy(solanaceae, solanum, tuberosum)
unidentified_animal <- hierarchy(mammalia, unidentified)
unidentified_plant <- hierarchy(plantae, unidentified)

info <- data.frame(name = c("tiger", "cat", "mole", "human", "tomato", "potato"),
                   n_legs = c(4, 4, 4, 2, 0, 0),
                   dangerous = c(TRUE, FALSE, FALSE, TRUE, FALSE, FALSE))

abund <- data.frame(code = rep(c("T", "C", "M", "H"), 2),
                    sample_id = rep(c("A", "B"), each = 2),
                    count = c(1,2,5,2,6,2,4,0),
                    taxon_index = rep(1:4, 2))

phylopic_ids <- c("e148eabb-f138-43c6-b1e4-5cda2180485a",
                  "12899ba0-9923-4feb-a7f9-758c3c7d5e13",
                  "11b783d5-af1c-4f4e-8ab5-a51470652b47",
                  "9fae30cd-fb59-4a81-a39c-e1826a35f612",
                  "b6400f39-345a-4711-ab4f-92fd4e22cb1a",
                  "63604565-0406-460b-8cb8-1abe954b3f3a")

foods <- list(c("mammals", "birds"),
              c("cat food", "mice"),
              c("insects"),
              c("Most things, but especially anything rare or expensive"),
              c("light", "dirt"),
              c("light", "dirt"))

reaction <- function(x) {
  ifelse(x$data$info$dangerous,
         paste0("Watch out! That ", x$data$info$name, " might attack!"),
         paste0("No worries; its just a ", x$data$info$name, "."))
}

### Make test object

test_obj <- taxmap(tiger, cat, mole, human, tomato, potato,
                   data = list(info = info,
                               phylopic_ids = phylopic_ids,
                               foods = foods,
                               abund = abund),
                   funcs = list(reaction = reaction))

### Manual class construction

test_that("Existing taxon_id column in table data", {
  expect_message(taxmap(tiger, cat, mole, human, tomato, potato,
                        data = list(x = data.frame(taxon_id = "b", x = 2))),
                 'Using existing "taxon_id" column for table')
})

test_that("Existing taxon_id with invalid IDs", {
  expect_error(taxmap(tiger, cat, mole, human, tomato, potato,
                      data = list(x = data.frame(taxon_id = "xxx", x = 2))),
               'The table "x" has a "taxon_id" column, but the values do not appear to be taxon IDs')
})

test_that("Existing taxon_index column in table data", {
  expect_message(taxmap(tiger, cat, mole, human, tomato, potato,
                        data = list(x = data.frame(taxon_index = 6, x = 2))),
                 'Using "taxon_index" column')
})

test_that("No taxon_id or taxon_index column in table data", {
  expect_warning(taxmap(tiger, cat, mole, human, tomato, potato,
                        data = list(x = data.frame(x = 2))),
                 'The table "x" does not have a "taxon_index" column')
})

test_that("Same length table data", {
  expect_message(taxmap(tiger, cat, mole, human, tomato, potato,
                        data = list(x = data.frame(x = rep(2, 6)))),
                 'Assuming that the elements of table "x" are in the same order')
})


test_that("vector/list named by taxon IDs", {
  expect_message(taxmap(tiger, cat, mole, human, tomato, potato,
                        data = list(x = c("b" = 2))),
                 'Using existing names of list/vector "x" as taxon IDs.')
})

test_that("vector/list named, but not by taxon IDs", {
  expect_warning(taxmap(tiger, cat, mole, human, tomato, potato,
                      data = list(x = c("xxx" = 2))),
               'The list/vector "x" is named, but the names do not appear to be taxon IDs.')
})


test_that("No names for vector/list data set", {
  expect_warning(taxmap(tiger, cat, mole, human, tomato, potato,
                        data = list(x = 2)),
                 'The list/vector "x" is unnamed so has no taxon ID information.')
})

test_that("Same length vector/list data set", {
  expect_message(taxmap(tiger, cat, mole, human, tomato, potato,
                        data = list(x = rep(2, 6))),
                 'Assuming that the elements of list/vector "x" are in the same order')
})


### Print methods

test_that("Print methods works", {
  x = test_obj$clone(deep = TRUE)
  x$data <- list()
  x$data$more_data <- list(a = 1, b = 2, c = 3)
  x$data$frame <- data.frame(x = 1:10)
  x$data$mat <- matrix(1:9, nrow = 3)
  x$data$fac <- factor(1:10)
  x$data$tib <- dplyr::as_tibble(data.frame(x = 1:10))
  expect_output(print(x),
                "<Taxmap>.+17 taxa.+17 edges.+1 functions.+reaction")
  x$data$new_vec <- rep(paste0(rep(c("l", "o", "n", "g"), each = 30), collapse = ""), 5)
  expect_output(print(x),
                "\\[truncated\\] \\.\\.\\. ")
  x$data <- c(x$data, list(1:100, a = 1:10, b = 1:100))
  expect_output(print(x), "more data sets")

  # No taxa
  x <- taxmap()
  expect_output(print(x), "No taxa")

  # Empty list
  x = test_obj$clone(deep = TRUE)
  x$data <- list()
  x$data$more_data <- list()
  x$data[[2]] <- list()
  expect_output(print(x), "empty list")

  # List named by taxa
  x = test_obj$clone(deep = TRUE)
  x$data <- list()
  x$data$more_data <- list(c = 3, d = 4)
  expect_output(print(x), "named by taxa")

  # Named vectors
  x = test_obj$clone(deep = TRUE)
  x$data <- list()
  x$data$more_data <- c(sss = 3, dddd = 4)
  x$data$even_more <- c(c = 3, d = 4)
  expect_output(print(x), "named vector")

  # Vector types
  x = test_obj$clone(deep = TRUE)
  x$data <- list()
  x$data$int <- as.integer(1:10)
  x$data$char <- as.character(1:10)
  x$data$fac <- as.factor(1:10)
  x$data$ord <- as.ordered(1:10)
  x$data$log <- as.logical(1:10)
  expect_output(print(x), "integer")
  expect_output(print(x), "character")
  expect_output(print(x), "factor")
  expect_output(print(x), "ordered")
  expect_output(print(x), "logical")


})

### NSE helpers

#### all_names

test_that("Names of table col names are accessible by NSE", {
  expected_col_names <- colnames(test_obj$data$info)
  expected_col_names <- expected_col_names[expected_col_names != "taxon_id"]
  expect_true(all(expected_col_names %in% test_obj$all_names()))
  expect_false(any(expected_col_names %in% test_obj$all_names(tables = FALSE)))
})

test_that("Names of functions are accessible by NSE", {
  expected_funcs <- names(test_obj$funcs)
  expect_true(all(expected_funcs %in% test_obj$all_names()))
  expect_false(any(expected_funcs %in% test_obj$all_names(funcs = FALSE)))
})

test_that("Names of functions are accessible by NSE", {
  expected_others <-
    names(test_obj$data)[sapply(test_obj$data,
                                function(x) ! "data.frame" %in% class(x))]
  expect_true(all(expected_others %in% test_obj$all_names()))
  expect_false(any(expected_others %in% test_obj$all_names(others = FALSE)))
})

test_that("Names of built-in functions are accessible by NSE", {
  expected_builtin <- c("taxon_names", "taxon_ids", "n_supertaxa", "n_subtaxa",
                        "n_subtaxa_1")
  expect_true(all(expected_builtin %in% test_obj$all_names()))
  expect_false(any(expected_builtin %in%
                     test_obj$all_names(builtin_funcs = FALSE)))
})

test_that("Duplicate names give a warning", {
  x = test_obj$clone(deep = TRUE)
  x$data$n_legs = 1:10
  expect_warning(x$all_names(warn = T),
                 "The following names are used more than once: n_legs")
  expect_silent(x$all_names(warn = F))
})

#### names_used

test_that("Names in basic expressions can be found by NSE", {
  expect_true(all(c("n_subtaxa", "n_legs", "taxon_ids")
                  %in% test_obj$names_used(n_legs + n_subtaxa, taxon_ids + 19)))
})

test_that("Names in complex expressions can be found by NSE", {
  expect_true(all(c("n_subtaxa", "n_legs", "taxon_ids", "dangerous", "reaction") %in%
                    test_obj$names_used((((((n_legs))))),
                                        function(x) length(taxon_ids) + x,
                                        {{n_subtaxa}},
                                        taxon_ids[n_subtaxa[dangerous]],
                                        reaction(n_subtaxa))))
})

test_that("Names in invalid expressions can be found by NSE", {
  expect_true(all(c("n_subtaxa")
                  %in% test_obj$names_used(not_a_variable == n_subtaxa,
                                           aslkadsldsa)))
})

test_that("Names of varaibles referred to by full $ path are not returned", {
  expect_equal(length(test_obj$names_used(data$abund$count)), 0)
  expect_equal(length(test_obj$names_used(data$count)), 0)
  expect_equal(length(test_obj$names_used(count)), 1)
})

#### get_data

test_that("NSE values can be found", {
  expect_equal(test_obj$get_data(c("n_subtaxa", "n_legs", "reaction")),
               list(n_subtaxa = test_obj$n_subtaxa(),
                    n_legs = stats::setNames(test_obj$data$info$n_legs,
                                             test_obj$data$info$taxon_id),
                    reaction = test_obj$funcs$reaction(test_obj)))
  expect_error(test_obj$get_data(c("n_subtaxa", "not_valid")),
               "Cannot find the following data: not_valid")
})

test_that("All valid NSE values can be found", {
  expect_equal(names(get_data(test_obj)), unname(all_names(test_obj)))
})


test_that("Using ambiguous names in NSE generates a warning", {
  expect_equal(names(get_data(test_obj)), unname(all_names(test_obj)))
})


#### Get datasets

test_that("datasets can be accessed", {
  # Works right with valid input
  expect_identical(get_dataset(test_obj, "info"), test_obj$data$info)
  expect_identical(get_dataset(test_obj, 1), test_obj$data$info)
  expect_identical(get_dataset(test_obj,  names(test_obj$data) == "info"), test_obj$data$info)

  # Fails with invalid input
  expect_error(get_dataset(test_obj, "not valid"), 'The dataset "not valid" cannot be found')
  expect_error(get_dataset(test_obj, 123), 'The dataset "123" cannot be found')
  expect_error(get_dataset(test_obj, TRUE), 'must be the same length')
})


#### get_data_frame

test_that("get data frame - for now doesn't work on example data", {
  x <- test_obj$clone(deep = TRUE)
  x$data$abund_2 <- x$data$abund
  expect_warning(x$get_data(), "Ambiguous names used")
  expect_warning(x$get_data(c("count", "code")), "Ambiguous names used")
  expect_warning(x$filter_obs("abund", code == "T"), "Ambiguous names used")

})



### Mapping functions

#### obs

test_that("Mapping between table observations and the edge list works", {
  result <- obs(test_obj, "info")
  expect_equal(obs(test_obj, "info"), test_obj$obs("info"))
  expect_true(all(sapply(result, class) == "integer"))
  expect_identical(names(result), unname(test_obj$taxon_ids()))
  expect_identical(result[["b"]], 1:4)
  expect_identical(result, test_obj$obs("phylopic_ids"))
  expect_identical(result, test_obj$obs("foods"))
})

test_that("Returning values for observations", {
  expect_true(all(c(4, 4, 4, 2) %in%
                    obs(test_obj, "info", subset = "b", value = "n_legs", simplify = T)))
})

test_that("Mapping between a subset of observations and the edge list works", {
  expect_identical(obs(test_obj, "info", subset = "b"), list("b" = 1:4))
  expect_identical(obs(test_obj, "info", subset = 1), list("b" = 1:4))
})

test_that("Mapping non-recursivly between observations and the edge list works", {
  result <- obs(test_obj, "info", recursive = FALSE)
  expect_true(all(sapply(result[roots(test_obj)], length) == 0))
  expect_equal(result[["r"]], 6)
})

test_that("Mapping simplification between observations and the edge list works", {
  expect_equal(obs(test_obj, "info", simplify = TRUE), 1:6)
})

test_that("Mapping observations in external variables", {
  external_table <- data.frame(taxon_id = c("p", "n"),
                               my_name = c("Joe", "Fluffy"))
  expect_equal(eval(substitute(obs(test_obj, external_table)$`b`)), c(2, 1))

  external_table <- data.frame(my_name = c("Joe", "Fluffy"))
  expect_error(eval(substitute(obs(test_obj, external_table))), 'no "taxon_id" column')

  extern_vec <- c(p = "Joe", n = "Fluffy")
  expect_equal(eval(substitute(obs(test_obj, extern_vec)$`b`)), c(2, 1))

  extern_vec <- c("Joe", "Fluffy")
  expect_error(eval(substitute(obs(test_obj, extern_vec))), 'no taxon ids')
})

test_that("Mapping observations when there are multiple obs per taxon", {
  result <- obs(test_obj, "abund")
  expect_equal(result$m, which(test_obj$data$abund$taxon_id == "m"))
  expect_equal(result$p, which(test_obj$data$abund$taxon_id == "p"))
  expect_true(all(result$b %in% 1:nrow(test_obj$data$abund)))
})

test_that("Applying a function the observations of each taxon", {
  expect_equal(obs_apply(test_obj, "abund", length),
               lapply(obs(test_obj, "abund"), length))

})

test_that("Counting the observations of each taxon", {
  expect_equal(n_obs(test_obj, "abund"),
               sapply(obs(test_obj, "abund"), length))
  expect_equal(n_obs_1(test_obj, "abund"),
               vapply(taxon_ids(test_obj), function(id) sum(id == test_obj$data$abund$taxon_id), numeric(1)))

})



### Dplyr analogs

#### filter_taxa

test_that("Default taxon filtering works", {
  result <- filter_taxa(test_obj, taxon_names == "Solanum")
  expect_equal(result$taxon_names(), c("l" = "Solanum"))
  expect_equal(as.character(result$data$info$name), c("tomato", "potato"))
  expect_true(length(result$data$phylopic_ids) == 2)
  expect_true(length(result$data$foods) == 2)
})

test_that("Subtaxa can be included when filtering taxa", {
  result <- filter_taxa(test_obj, taxon_names == "Solanum", subtaxa = TRUE)
  expect_equivalent(result$taxon_names(),
                    c("Solanum", "lycopersicum", "tuberosum"))
  expect_equal(filter_taxa(test_obj, 1, subtaxa = TRUE),
               filter_taxa(test_obj, 1, subtaxa = -1))
  expect_equal(filter_taxa(test_obj, 1, subtaxa = FALSE),
               filter_taxa(test_obj, 1, subtaxa = 0))
})

test_that("Supertaxa can be included when filtering taxa", {
  result <- filter_taxa(test_obj, taxon_names == "Solanum", supertaxa = TRUE)
  expect_equivalent(sort(result$taxon_names()),
                    sort(c("Solanum", "Solanaceae", "Plantae")))
  expect_equal(filter_taxa(test_obj, 16, supertaxa = TRUE),
               filter_taxa(test_obj, 16, supertaxa = -1))
  expect_equal(filter_taxa(test_obj, 16, supertaxa = FALSE),
               filter_taxa(test_obj, 16, supertaxa = 0))
})

test_that("Observations can be preserved when filtering taxa", {
  result <- filter_taxa(test_obj, taxon_names == "Solanum", reassign_obs = FALSE)
  expect_equal(nrow(result$data$info), 0)
  expect_equal(length(result$data$phylopic_ids), 0)
  expect_equal(length(result$data$foods), 0)
  result <- filter_taxa(test_obj, taxon_names == "tuberosum", reassign_obs = FALSE)
  expect_equivalent(result$taxon_names(), "tuberosum")
  result <- filter_taxa(test_obj, taxon_names == "Solanum", drop_obs = FALSE)
  expect_equal(result$data$info$taxon_id, c(NA, NA, NA, NA, "l", "l"))
  expect_equal(names(result$data$phylopic_ids), c(NA, NA, NA, NA, "l", "l"))
  expect_equal(names(result$data$foods), c(NA, NA, NA, NA, "l", "l"))
  result <- filter_taxa(test_obj, taxon_names == "Solanum", drop_obs = FALSE,
                        reassign_obs = FALSE)
  expect_true(all(is.na(result$data$info$taxon_id)))
  expect_equal(nrow(result$data$info), 6)
  expect_equal(filter_taxa(test_obj, 2:4, drop_obs = TRUE),
               filter_taxa(test_obj, 2:4, drop_obs = c(TRUE, TRUE, TRUE, TRUE)))
  expect_equal(filter_taxa(test_obj, 2:4, drop_obs = TRUE),
               filter_taxa(test_obj, 2:4, drop_obs = c(info = TRUE,
                                                       phylopic_ids = TRUE,
                                                       foods = TRUE,
                                                       abund = TRUE)))
})

test_that("Taxon ids can be preserved when filtering taxa", {
  result <- filter_taxa(test_obj, taxon_names != "Solanum", reassign_taxa = FALSE)
  expect_true(all(c("lycopersicum", "tuberosum") %in% result$roots(value = "taxon_names")))
})

test_that("The selection of taxa to be filtered can be inverted", {
  result <- filter_taxa(test_obj, taxon_names == "Solanum", subtaxa = TRUE, invert = TRUE)
  expect_true(all(! c("tuberosum", "lycopersicum", "Solanum") %in% result$taxon_names()))
  expect_true(all(c("Mammalia", "Plantae", "sapiens") %in% result$taxon_names()))
})

test_that("Edge cases return reasonable outputs", {
  expect_equal(filter_taxa(test_obj), test_obj)
  expect_error(filter_taxa(test_obj, drop_obs = c(TRUE, TRUE)),
               'Invalid input for logical vector selecting')
  expect_error(filter_taxa(test_obj, drop_obs = c(not_valid = TRUE)),
               'Invalid input for logical vector selecting')
  expect_error(filter_taxa(test_obj, drop_obs = c(not_valid = TRUE, TRUE, TRUE, FALSE)),
               'Invalid input for logical vector selecting')
})

test_that("Filtering taxa when there are multiple obs per taxon", {
  result <- filter_taxa(test_obj, taxon_names == "Solanum")
  expect_equal(nrow(result$data$abund), 0) # There were no plants in that set

  result <- filter_taxa(test_obj, taxon_names == "Felidae", subtaxa = TRUE)
  expect_equal(nrow(result$data$abund), 4) # There were 2 cats and 2 tigers
})

#### filter_obs

test_that("Default observation filtering works", {
  result <- filter_obs(test_obj, "info", n_legs == 2, dangerous == TRUE)
  expect_equivalent(as.character(result$data$info$name), "human")
  result <- filter_obs(test_obj, "phylopic_ids", n_legs == 2, dangerous == TRUE)
  expect_equal(length(result$data$phylopic_ids), 1)
  result <- filter_obs(test_obj, "foods", n_legs == 2, dangerous == TRUE)
  expect_equal(result$data$foods[[1]],
               "Most things, but especially anything rare or expensive")
})

test_that("Datasets can be specified using names, numbers, and characters", {
  result <- filter_obs(test_obj, c("info", "foods"), n_legs == 2)
  expect_equal(result, filter_obs(test_obj, c(1, 3), n_legs == 2))
  expect_equal(result, filter_obs(test_obj, 1:4 %in% c(1, 3), n_legs == 2))
})


test_that("Filtering observations with external variables work", {
  my_logical <- test_obj$data$info$n_legs == 2 & test_obj$data$info$dangerous == TRUE
  expect_equal(filter_obs(test_obj, "info", n_legs == 2, dangerous == TRUE),
               filter_obs(test_obj, "info", my_logical))
  expect_equal(filter_obs(test_obj, "info", n_legs == 2, dangerous == TRUE, drop_taxa = TRUE),
               filter_obs(test_obj, "info", my_logical, drop_taxa = TRUE))

  result <- filter_obs(test_obj, "info", n_legs == 2, dangerous == TRUE)
  expect_equivalent(as.character(result$data$info$name), "human")
  result <- filter_obs(test_obj, "phylopic_ids", n_legs == 2, dangerous == TRUE)
  expect_equal(length(result$data$phylopic_ids), 1)
  result <- filter_obs(test_obj, "foods", n_legs == 2, dangerous == TRUE)
  expect_equal(result$data$foods[[1]],
               "Most things, but especially anything rare or expensive")
})

test_that("Removing taxa when filtering observations work", {

  result <- filter_obs(test_obj, "info", n_legs == 2, drop_taxa = TRUE)
  expect_equivalent(as.character(result$data$info$name), "human")
  expect_equivalent(sort(result$taxon_names()),
                    sort(c("Mammalia", "Hominidae", "homo", "sapiens")))
  expect_equal(names(result$data$phylopic_ids), result$data$info$taxon_id)
  expect_equal(names(result$data$foods), result$data$info$taxon_id)
  expect_equal(unique(result$data$abund$taxon_id), result$data$info$taxon_id)

  # Removing taxa that appear in some datasets
  result <- filter_obs(test_obj, "info", n_legs == 2, drop_taxa = TRUE,
                       drop_obs = c(abund = FALSE))
  expect_equal(result$data$abund$taxon_id, test_obj$data$abund$taxon_id)
  expect_equivalent(result$roots(value = "taxon_names"), "Mammalia")
  expect_true(length(test_obj$taxa) > length(result$taxa))

})

test_that("Edge cases return reasonable outputs", {
  expect_equal(filter_obs(test_obj, "info"), test_obj)
  expect_error(filter_obs(test_obj, "not_valid",
                          "not the name of a data set. Valid targets "))
  expect_error(filter_obs(test_obj, "info", "11"),
               "observation filtering with taxon IDs is not currently")
  expect_error(filter_obs(test_obj, character(0)),
               "At least one dataset must be specified.")

})

test_that("Filtering obs when there are multiple obs per taxon", {
  result <- filter_obs(test_obj, "abund", code == "C", drop_taxa = TRUE)
  expect_equal(nrow(result$data$abund), 2)
})


test_that("Filtering multiple datasets at once", {
  result <- filter_obs(test_obj, c("phylopic_ids", "info"), n_legs < 4, drop_taxa = TRUE)
  expect_equal(length(result$data$phylopic_ids), 3)
  expect_equal(nrow(result$data$info), 3)

  # Multiple datasets with different taxon IDs
  test_obj_2 <- test_obj$clone(deep = TRUE)
  test_obj_2$data$abund_2 <- test_obj_2$data$abund
  test_obj_2$data$abund_2$taxon_id <- rep("r", 8)
  result <- expect_warning(filter_obs(test_obj_2, c("abund", "abund_2"), code == "C", drop_taxa = TRUE))
  expect_true("tuberosum" %in% taxon_names(result))
  expect_true(! "sapiens" %in% taxon_names(result))
  expect_equal(nrow(result$data$abund), 2)
  expect_equal(nrow(result$data$abund_2), 2)

  # Datasets of different length cannot be filtered
  expect_error(filter_obs(test_obj, c("phylopic_ids", "abund"), n_legs < 4, drop_taxa = TRUE),
               "If multiple datasets are filtered at once, then they must the same length")
})


#### select_obs

test_that("Default observation column subsetting works",  {
  result <- select_obs(test_obj, "info", dangerous)
  expect_equal(colnames(result$data$info), c("taxon_id", "dangerous"))
})

test_that("Edge cases return reasonable outputs during observation column subsetting", {
  result <- select_obs(test_obj, "info")
  expect_equal(colnames(result$data$info), c("taxon_id"))
  expect_error(select_obs(test_obj, "not_valid"),
               "The input does not correspond to a valid dataset")
  expect_error(select_obs(test_obj), " missing, with no default")
  expect_error(select_obs(test_obj, "foods"), 'not a table, so columns cannot be selected')
})

test_that("The columns of multiple datasets can be subset at once", {
  test_obj_2 <- test_obj$clone(deep = TRUE)
  test_obj_2$data$abund_2 <- test_obj_2$data$abund
  result <- select_obs(test_obj_2, c("abund", "abund_2"), count, code)
  expect_equal(colnames(result$data$abund), c("taxon_id", "count", "code"))
  expect_equal(colnames(result$data$abund_2), c("taxon_id", "count", "code"))
})


#### mutate_obs

test_that("Observation column addition works",  {
  result <- mutate_obs(test_obj, "info",
                       new_col = "new",
                       newer_col = paste0(new_col, "er"))
  expect_true(all(c("new_col", "newer_col") %in% colnames(result$data$info)))
})

test_that("Observation column replacement works",  {
  result <- mutate_obs(test_obj, "info", name = "replacement")
  expect_true(all(result$data$info$name == "replacement"))
})

test_that("Edge cases for observation column addition",  {
  expect_equal(mutate_obs(test_obj, "info"), test_obj)
  expect_error(select_obs(test_obj, "foods"), 'not a table, so columns cannot be selected.')
  expect_error(select_obs(test_obj, "phylopic_ids"),
               'is not a table, so columns cannot be selected.')
})

test_that("New tables and vectors can be made",  {
 # New tables
  result <- mutate_obs(test_obj, "new_table",
                       ranks = taxon_ranks,
                       new_col = "new",
                       newer_col = paste0(new_col, "er"))
  expect_equal(dim(result$data$new_table), c(17, 3))
  result <- mutate_obs(test_obj, "new_table", a = 1:10)
  expect_equal(dim(result$data$new_table), c(10, 1))
  result <- mutate_obs(test_obj, "new_table", a = numeric(0), b = character(0))
  expect_equal(dim(result$data$new_table), c(0, 2))

  # New vectors
  result <- mutate_obs(test_obj, "new_table", 1:10)
  expect_equal(length(result$data$new_table), 10)
  result <- mutate_obs(test_obj, "new_table", numeric(0))
  expect_equal(length(result$data$new_table), 0)

 # Invlaid: inputs of mixed lengths
  expect_error(mutate_obs(test_obj, "new_table", a = 1:3, b = 2:8),
               "Cannot make a new table out of multiple values of unequal length")

 # Invlaid: unnamed inputs
  expect_error(mutate_obs(test_obj, "new_table", 1:10, 1:10),
               "Cannot add a new dataset with")

 # invalid: not a table
  expect_error(mutate_obs(test_obj, "foods", 1:10),
               "is not a table")

})

#### transmute_obs

test_that("Observation column addition (transmute) works",  {
  result <- transmute_obs(test_obj, "info",
                          new_col = paste("new", name),
                          newer_col = paste0(new_col, "!!"))
  expect_equal(c("taxon_id", "new_col", "newer_col"),
               colnames(result$data$info))
  x <- test_obj$clone(deep = TRUE)
  x$data$new_table <- data.frame(y = 1:4)
  result <- transmute_obs(x, "new_table", # no taxon ids in new_table
                          new_col = paste("new", name),
                          newer_col = paste0(new_col, "!!"))
  expect_equal(c("new_col", "newer_col"),
               colnames(result$data$new_table))
})

test_that("Edge cases for observation column addition (transmute) ",  {
  result <- transmute_obs(test_obj, "info")
  expect_equal("taxon_id", colnames(result$data$info))
  expect_error(transmute_obs(test_obj, "not_valid"),
               "The input does not correspond to a valid dataset")
  expect_error(select_obs(test_obj, "foods"), 'not a table, so columns cannot be selected')
  expect_error(select_obs(test_obj, "phylopic_ids"),
               'not a table, so columns cannot be selected')
})


#### arrange_obs

test_that("Sorting observations work",  {
  result <- arrange_obs(test_obj, "info", dangerous, name)
  expect_equal(test_obj$data$info$taxon_id[order(test_obj$data$info$dangerous,
                                                 test_obj$data$info$name)],
               result$data$info$taxon_id)
  result <- arrange_obs(test_obj, "info", desc(dangerous), desc(name))
  expect_equal(test_obj$data$info$taxon_id[order(test_obj$data$info$dangerous,
                                                 test_obj$data$info$name,
                                                 decreasing = TRUE)],
               result$data$info$taxon_id)
  list_results <- arrange_obs(test_obj, "foods", desc(dangerous), desc(name))
  expect_equal(result$data$info$taxon_id, names(list_results$data$foods))
  list_results <- arrange_obs(test_obj, "phylopic_ids", desc(dangerous), desc(name))
  expect_equal(result$data$info$taxon_id, names(list_results$data$phylopic_ids))
})

test_that("Sorting observations with non-target NSE values",  {
  result <- arrange_obs(test_obj, "info", phylopic_ids)
  expect_equal(test_obj$data$info$taxon_id[order(test_obj$data$phylopic_ids)],
               result$data$info$taxon_id)
})

test_that("Edge cases during observation sorting works",  {
  expect_equal(arrange_obs(test_obj, "info"), test_obj)
  expect_error(arrange_obs(test_obj, "not_valid"),
               "The input does not correspond to a valid dataset")
})

test_that("Sorting multiple datasets works",  {
  result <- arrange_obs(test_obj, c("info", "phylopic_ids", "foods"), n_legs)
  expect_equal(test_obj$data$info$taxon_id[order(test_obj$data$info$n_legs)],
               result$data$info$taxon_id)
  expect_equal(result$data$info$taxon_id, names(result$data$phylopic_ids))
  expect_equal(result$data$info$taxon_id, names(result$data$foods))
})


#### arrange_taxa

test_that("Sorting taxa work",  {
  expect_equal(arrange_taxa(test_obj, taxon_ids)$taxon_ids(),
               sort(test_obj$taxon_ids()))
  expect_equal(arrange_taxa(test_obj, desc(taxon_ids))$taxon_ids(),
               sort(test_obj$taxon_ids(), decreasing = TRUE))
})

test_that("Edge cases during observation sorting works",  {
  expect_equal(arrange_taxa(test_obj), test_obj)
})


#### sample_n_obs

test_that("Sampling observations works",  {
  result <- sample_n_obs(test_obj, "info", size = 3)
  expect_equal(nrow(result$data$info), 3)
  result <- sample_n_obs(test_obj, "info", size = 30, replace = TRUE)
  expect_equal(nrow(result$data$info), 30)
  result <- sample_n_obs(test_obj, "foods", size = 3)
  expect_equal(length(result$data$foods), 3)
  result <- sample_n_obs(test_obj, "phylopic_ids", size = 3)
  expect_equal(length(result$data$phylopic_ids), 3)

  result <- sample_frac_obs(test_obj, "info", size = 0.5)
  expect_equal(nrow(result$data$info), 3)

  result <- sample_frac_obs(test_obj, "info", size = 0.5, taxon_weight = 1 / n_obs)
  expect_equal(nrow(result$data$info), 3)


  result <- sample_frac_obs(test_obj, "info", size = 0.5, taxon_weight = 1 / n_obs)
  expect_equal(nrow(result$data$info), 3)
})

test_that("Sampling using data from supertaxa works",  { # Not complete
  expect_equal({
    set.seed(1)
    sample_n_obs(test_obj, "info", size = 3, use_supertaxa = 0)
  },
  {
    set.seed(1)
    sample_n_obs(test_obj, "info", size = 3, use_supertaxa = FALSE)
  })
  expect_equal({
    set.seed(1)
    sample_n_obs(test_obj, "info", size = 3, use_supertaxa = -1)
  },
  {
    set.seed(1)
    sample_n_obs(test_obj, "info", size = 3, use_supertaxa = TRUE)
  })
})


test_that("Edge cases during sampling observations",  {
  expect_error(sample_n_obs(test_obj),
               "missing, with no default")
  expect_error(sample_n_obs(test_obj, "not_valid"),
               "The input does not correspond to a valid dataset.")
})


#### sample_n_taxa

test_that("Sampling taxa works",  {
  result <- sample_n_taxa(test_obj, size = 3)
  expect_equal(length(result$taxon_ids()), 3)
  expect_error(sample_n_taxa(test_obj, obs_weight = 1:10),
               "`obs_target` must also be defined.")

  result <- sample_n_taxa(test_obj, 3, obs_target = "info", obs_weight = 1:6)
  expect_equal(length(result$taxon_ids()), 3)
})


test_that("Edge cases during sampling observations",  {
  expect_error(sample_n_taxa(test_obj),
               "missing, with no default")
  expect_error(sample_n_taxa(),
               "missing, with no default")
})


test_that("Sampling observations using data from subtaxa works", { # Not complete
  expect_equal({
    set.seed(1)
    sample_n_taxa(test_obj, size = 3, use_subtaxa = 0)
  },
  {
    set.seed(1)
    sample_n_taxa(test_obj, size = 3, use_subtaxa = FALSE)
  })
  expect_equal({
    set.seed(1)
    sample_n_taxa(test_obj, size = 3, use_subtaxa = -1)
  },
  {
    set.seed(1)
    sample_n_taxa(test_obj, size = 3, use_subtaxa = TRUE)
  })
})


test_that("dots and .list return the same output", {
  expect_equal(taxmap(tiger, cat, mole, human, tomato, potato,
                      data = list(info = info,
                                  phylopic_ids = phylopic_ids,
                                  foods = foods),
                      funcs = list(reaction = reaction)),
               taxmap(.list = list(tiger, cat, mole, human, tomato, potato),
                      data = list(info = info,
                                  phylopic_ids = phylopic_ids,
                                  foods = foods),
                      funcs = list(reaction = reaction)))
})

