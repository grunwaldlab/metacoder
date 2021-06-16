#' Taxmap class
#'
#' A class designed to store a taxonomy and associated information. This class
#' builds on the [taxonomy()] class. User defined data can be stored in the list
#' `obj$data`, where `obj` is a taxmap object. Data that is associated with taxa
#' can be manipulated in a variety of ways using functions like [filter_taxa()]
#' and [filter_obs()]. To associate the items of lists/vectors with taxa, name
#' them by [taxon_ids()]. For tables, add a column named `taxon_id` that stores
#' [taxon_ids()].
#'
#' To initialize a `taxmap` object with associated data sets, use the parsing
#' functions [parse_tax_data()], [lookup_tax_data()], and [extract_tax_data()].
#'
#' @export
#' @param data A list of tables with data associated with the taxa.
#' @param funcs A named list of functions to include in the class. Referring to
#'   the names of these in functions like [filter_taxa()] will execute the
#'   function and return the results. If the function has at least one argument,
#'   the taxmap object is passed to it.
#' @inheritParams taxonomy
#'
#' @family classes
#' @return An `R6Class` object of class [taxmap()]
#'
#' @details on initialize, function sorts the taxon list based on rank (if rank
#'   information is available), see [ranks_ref] for the reference rank names and
#'   orders
#'
#' @template taxmapegs
taxmap <- function(..., .list = NULL, data = NULL, funcs = list(), named_by_rank = FALSE) {
  Taxmap$new(..., .list = .list, data = data, funcs = funcs, named_by_rank = named_by_rank)
}

Taxmap <- R6::R6Class(
  "Taxmap",
  inherit = Taxonomy,
  public = list(
    data = list(),
    funcs = list(),

    # -------------------------------------------------------------------------
    # Constructor
    initialize = function(..., .list = NULL, data = list(), funcs = list(), named_by_rank = FALSE) {

      # Call `taxonomy` constructor
      super$initialize(..., .list = .list, named_by_rank = named_by_rank)

      # Make sure `data` is in the right format and add to object
      self$data <- init_taxmap_data(self, data, self$input_ids)
      check_taxmap_data(self)

      # Make sure `funcs` is in the right format and add to object
      self$funcs <- validate_taxmap_funcs(funcs)
    },

    # -------------------------------------------------------------------------
    print = function(indent = "", max_rows = 3, max_items = 6,
                     max_width = getOption("width") - 10) {

      # Print taxonomy
      cat(paste0(indent, "<Taxmap>\n"))
      taxon_names <- vapply(self$taxa, function(x) x$name$name, character(1))
      taxon_ids <- names(self$taxa)
      if (length(self$taxa) > 0) {
        limited_print(paste(tid_font(taxon_ids), taxon_names,
                            sep = punc_font(". ")),
                      sep = punc_font(", "),
                      mid = punc_font(" ... "),
                      trunc_char = punc_font("[truncated]"),
                      prefix = paste0(indent, "  ",
                                      length(self$taxa), " taxa:"),
                      type = "cat")
        limited_print(private$make_graph(),
                      sep = punc_font(", "),
                      mid = punc_font(" ... "),
                      trunc_char = punc_font("[truncated]"),
                      prefix = paste0(indent, "  ",
                                      nrow(self$edge_list), " edges:"),
                      type = "cat")
      } else {
        cat("  No taxa\n  No edges\n")
      }

      # Get item names
      if (is.null(names(self$data))) {
        data_names <- paste0("[[", seq_len(length(self$data)), "]]")
      } else {
        data_names <- names(self$data)
        data_names[data_names == ""] <- paste0("[[", which(data_names == ""),
                                               "]]")
      }

      # Print a subset of each item, up to a max number, then just print names
      cat(paste0("  ", length(self$data), " data sets:\n"))
      if (length(self$data) > 0) {
        for (i in 1:min(c(max_items, length(self$data)))) {
          print_item(self, self$data[[i]],
                     name = data_names[i], max_rows = max_rows,
                     max_width = max_width, prefix = "    ")
        }
        if (length(self$data) > max_items) {
          cat(paste0("    And ", length(self$data) - max_items,
                     " more data sets:"))
          limited_print(data_names[(max_items + 1):length(self$data)],
                        type = "cat")
        }
      }

      # Print the names of functions
      cat(paste0("  ", length(self$funcs), " functions:\n"))
      limited_print(prefix = "   ", names(self$funcs), type = "cat")

      invisible(self)
    },

    # -------------------------------------------------------------------------
    # Check that a set of IDs are valid taxon IDs
    is_taxon_id = function(ids) {
      valid_ids <- c(unlist(self$edge_list), names(self$taxa), NA)
      ids %in% valid_ids
    },

    # -------------------------------------------------------------------------
    # Returns the names of things to be accessible using non-standard evaluation
    all_names = function(tables = TRUE, funcs = TRUE, others = TRUE,
                         builtin_funcs = TRUE, warn = FALSE) {
      output <- character(0)

      # Add functions included in the package
      if (builtin_funcs) {
        output <- c(output, private$nse_accessible_funcs)
      }

      # Get column names in each table, removing 'taxon_id'
      is_table <- vapply(self$data, is.data.frame, logical(1))
      if (tables && length(self$data[is_table]) > 0) {
        table_col_names <- unlist(lapply(self$data[is_table], colnames))
        names(table_col_names) <- paste0("data[['",
                                         rep(names(self$data[is_table]),
                                             vapply(self$data[is_table],
                                                    ncol, integer(1))),
                                         "']]")
        table_col_names <- table_col_names[table_col_names != "taxon_id"]
        output <- c(output, table_col_names)
      }

      # Get other object names in data
      is_other <- !is_table
      if (others && length(self$data[is_other]) > 0) {
        other_names <- names(self$data[is_other])
        names(other_names) <- rep("data", length(other_names))
        output <- c(output, other_names)
      }


      # Get function names
      if (funcs && length(self$funcs) > 0) {
        func_names <- names(self$funcs)
        names(func_names) <- rep("funcs", length(func_names))
        output <- c(output, func_names)
      }

      # Check for duplicates
      if (warn) {
        duplicated_names <- unique(output[duplicated(output)])
        if (length(duplicated_names) > 0) {
          warning(paste0("The following names are used more than once: ",
                         paste0(duplicated_names, collapse = ", ")))
        }
      }


      # Add the name to the name of the name and return
      names(output) <- paste0(names(output),
                              ifelse(names(output) == "", "", "[['"),
                              output,
                              ifelse(names(output) == "", "", "']]"))
      return(output)
    },


    # -------------------------------------------------------------------------
    # Get data indexes or other values associated with taxa
    obs = function(data, value = NULL, subset = NULL, recursive = TRUE,
                   simplify = FALSE) {

      # non-standard argument evaluation
      data_used <- eval(substitute(self$data_used(subset)))
      subset <- lazyeval::lazy_eval(lazyeval::lazy(subset), data = data_used)
      subset <- private$parse_nse_taxon_subset(subset)
      obs_taxon_ids <- self$get_data_taxon_ids(data, require = TRUE)

      # Get observations of taxa
      if (is.logical(recursive) && recursive == FALSE) {
        recursive = 0
      }
      if (recursive || is.numeric(recursive)) {
        my_subtaxa <- self$subtaxa(subset = unname(subset),
                                   recursive = recursive,
                                   include_input = TRUE,
                                   value = "taxon_indexes")
        #unname is needed for some reason.. something to look into...
      } else {
        my_subtaxa <- subset
      }
      obs_taxon_index <- match(obs_taxon_ids, self$taxon_ids())
      obs_key <- split(seq_along(obs_taxon_ids), obs_taxon_index)
      output <- stats::setNames(
        lapply(my_subtaxa,function(x) unname(unlist(obs_key[as.character(x)]))),
        names(subset)
      )
      is_null <- vapply(output, is.null, logical(1))
      output[is_null] <- lapply(1:sum(is_null), function(x) numeric(0))

      # Look up values
      if (!is.null(value)) {
        possible_values <- self$get_data(value)[[1]]
        if (length(possible_values) != length(obs_taxon_ids)) {
          stop(call. = FALSE,
               'The value "', value, '" is not the same length as the data set "', data, '".')
        }
        if (! is.null(names(possible_values)) && any(names(possible_values) != names(obs_taxon_ids))) {
          stop(call. = FALSE,
               'The value "', value, '" is in a different order than the data set "', data,
               '" according to taxon IDs.')
        }
        output <- lapply(output, function(i) possible_values[i])
      }

      # Reduce dimensionality
      if (simplify) {
        output <- simplify(output)
      }

      return(output)
    },


    # -------------------------------------------------------------------------
    # Apply a function to data for the observations for each taxon.
    # This is similar to using obs() with lapply() or sapply().
    obs_apply = function(data, func, simplify = FALSE, value = NULL,
                         subset = NULL, recursive = TRUE, ...) {
      my_obs <- self$obs(data, simplify = FALSE, value = value,
                         subset = eval(substitute(subset)),
                         recursive = recursive)
      output <- lapply(my_obs, func, ...)
      if (simplify) {
        output <- simplify(output)
      }
      return(output)
    },

    # -------------------------------------------------------------------------
    # Filter data in a taxmap() object (in obj$data) with a set of conditions.
    filter_obs = function(data, ..., drop_taxa = FALSE, drop_obs = TRUE,
                          subtaxa = FALSE, supertaxa = TRUE,
                          reassign_obs = FALSE, target = NULL) {

      # Check for use of "target"
      if (! is.null(target)) {
        warning(call. = FALSE,
                'Use of "target" is depreciated. Use "data" instead.')
        data <- target
      }

      # Parse data option
      data <- parse_dataset(self, data)

      # Check that multiple datasets are the same length
      dataset_length <- vapply(self$data[data], FUN.VALUE = numeric(1),
                               function(one) {
                                 if (is.data.frame(one)) {
                                   return(nrow(one))
                                 } else {
                                   return(length(one))
                                 }
                               })
      dataset_length <- unique(dataset_length)
      if (length(dataset_length) > 1) {
        stop(call. = FALSE,
             "If multiple datasets are filtered at once, then they must the same length. ",
             "The following lengths were found for the specified datasets:\n",
             limited_print(type = "silent", prefix = "  ", dataset_length))
      }

      # non-standard argument evaluation
      selection <- lazyeval::lazy_eval(lazyeval::lazy_dots(...),
                                       data = self$data_used(...))

      # Parse drop_obs option
      drop_obs <- parse_possibly_named_logical(
        drop_obs,
        self$data,
        formals(self$filter_obs)$drop_obs
      )

      # If no selection is supplied, match all rows
      if (length(selection) == 0) {
        selection <- list(seq_len(dataset_length))
      }

      # convert taxon_ids to indexes
      is_char <- vapply(selection, is.character, logical(1))
      if (sum(is_char) > 0) {
        stop(paste("observation filtering with taxon IDs is not currently",
                   "supported. If you want to filter observation by taxon IDs,",
                   "use something like: `obj$data$my_target$taxon_ids %in%",
                   "my_subset`"))
      }

      # convert logical to indexes
      is_logical <- vapply(selection, is.logical, logical(1))
      for (one in selection[is_logical]) {
        if (length(one) != dataset_length) {
          stop(call. = FALSE,
               "All logical filtering criteria must be the same length as the data sets filtered.")
        }
      }
      selection[is_logical] <- lapply(selection[is_logical], which)

      # combine filters
      intersect_with_dups <- function(a, b) {
        rep(sort(intersect(a, b)), pmin(table(a[a %in% b]), table(b[b %in% a])))
      }
      selection <- Reduce(intersect_with_dups, selection)

      # Remove observations
      data_taxon_ids <- NULL
      for (one in data) {
        data_taxon_ids <- c(data_taxon_ids, self$get_data_taxon_ids(one, require = drop_taxa)[selection])
        private$remove_obs(data = one, indexes = selection)
      }

      # Remove unobserved taxa
      data_taxon_ids <- unique(data_taxon_ids)
      if (drop_taxa & ! is.null(data_taxon_ids)) {

        # dont remove taxa that appear in other data sets if they are not also filtered
        sets_to_keep_ids_from <- names(drop_obs[! drop_obs])
        other_ids_to_keep <- unique(unlist(lapply(sets_to_keep_ids_from,
                                                  self$get_data_taxon_ids)))
        taxon_ids_to_keep <- unique(c(data_taxon_ids, other_ids_to_keep))

        # Remove taxa that are not in the filtered data set
        self$filter_taxa(taxon_ids_to_keep, drop_obs = drop_obs,
                         subtaxa = subtaxa, supertaxa = supertaxa,
                         reassign_obs = reassign_obs)
      }

      return(self)
    },


    # -------------------------------------------------------------------------
    # Subsets columns in a data set
    select_obs = function(data, ..., target = NULL) {

      # Check for use of "target"
      if (! is.null(target)) {
        warning(call. = FALSE,
                'Use of "target" is depreciated. Use "data" instead.')
        data <- target
      }

      # Parse data option
      data <- parse_dataset(self, data)

      # Check that the datasets are tables
      for (one in data) {
        if (! is.data.frame(self$data[[one]])) {
          stop(paste0('"data" ', one, ' is not a table, so columns cannot be selected.'))
        }
      }

      # Subset columns
      for (one in data) {
        self$data[[one]] <-
          dplyr::bind_cols(self$data[[one]][ , c("taxon_id"), drop = FALSE],
                           dplyr::select(self$data[[one]], ...))

      }

      return(self)
    },


    # -------------------------------------------------------------------------
    # Add columns to tables in obj$data
    mutate_obs = function(data, ..., target = NULL) {

      # Check for use of "target"
      if (! is.null(target)) {
        warning(call. = FALSE,
                'Use of "target" is depreciated. Use "data" instead.')
        data <- target
      }

      # Parse data option
      dataset_index <- parse_dataset(self, data, must_be_valid = FALSE, needed = FALSE)

      # Check that only one data is specified
      if (length(data) > 1) {
        stop(call. = FALSE,
             'Only one data can be mutated at a time.')
      }

      # Get data used in expressions to add
      data_used <- self$data_used(...)
      unevaluated <- lazyeval::lazy_dots(...)

      # add columns
      if (length(dataset_index) > 0) {
        # Check that the data is a table
        if (! is.data.frame(self$data[[dataset_index]])) {
          stop(paste0('Dataset "', data, '" is not a table, so columns cannot be added'))
        } else {
          for (index in seq_along(unevaluated)) {
            new_col <- lazyeval::lazy_eval(unevaluated[index], data = data_used)
            data_used <- c(data_used, new_col) # Allows this col to be used in next cols
            self$data[[dataset_index]][[names(new_col)]] <- new_col[[1]]
          }
        }
      } else { # not a current data
        new_dataset <- list()
        for (index in seq_along(unevaluated)) {
          new_col <- lazyeval::lazy_eval(unevaluated[index], data = data_used)
          data_used <- c(data_used, new_col) # Allows this col to be used in next cols
          new_dataset <- c(new_dataset, new_col)
        }
        if (any(names(unevaluated) == "")) { # unnammed inputs cant be put in tables
          if (length(unevaluated) == 1) { # Add as a vector
            message('Adding a new "', class(new_dataset[[1]]),'" vector of length ', length(new_dataset[[1]]), '.')
            self$data[[data]] <- new_dataset[[1]]
          } else {
            stop(call. = FALSE,
                 "Cannot add a new dataset with multiple values if any are unnamed.",
                 " The following input indexes are unnamed:\n",
                 limited_print(which(names(unevaluated) == ""), type = "silent", prefix = "  "))
          }
        } else { # Try to put in new table
          part_lengths <- vapply(new_dataset, length, numeric(1))
          if (length(unique(part_lengths[part_lengths != 1])) == 1) { # All inputs are same length or 1
            new_dataset <- dplyr::as_tibble(new_dataset)
            message('Adding a new ', nrow(new_dataset), ' x ', ncol(new_dataset),
                    ' table called "', data, '"')
            self$data[[data]] <- dplyr::as_tibble(new_dataset)
          } else {
            stop(call. = FALSE,
                 "Cannot make a new table out of multiple values of unequal length.",
                 " The inputs have the following lengths:\n",
                 limited_print(part_lengths, type = "silent", prefix = "  "))

          }
        }
      }

      return(self)
    },


    # -------------------------------------------------------------------------
    # Replace columns of tables in obj$data
    transmute_obs = function(data, ..., target = NULL) {

      # Check for use of "target"
      if (! is.null(target)) {
        warning(call. = FALSE,
                'Use of "target" is depreciated. Use "data" instead.')
        data <- target
      }

      # Parse data option
      data <- parse_dataset(self, data)

      # Check that only one data is specified
      if (length(data) > 1) {
        stop(call. = FALSE,
             'Only one dataset can be transmuted at a time.')
      }

      # Check that the dataset is a table
      if (! is.data.frame(self$data[[data]])) {
        stop(paste0('Dataset "', data, '" is not a table, so columns cannot be selected.'))
      }

      if ("taxon_id" %in% colnames(self$data[[data]])) {
        result <- list(taxon_id = self$data[[data]]$taxon_id)
      } else {
        result <- list()
      }
      data_used <- self$data_used(...)
      unevaluated <- lazyeval::lazy_dots(...)
      for (index in seq_along(unevaluated)) {
        new_col <- lazyeval::lazy_eval(unevaluated[index], data = data_used)
        # Allow this col to be used in evaluating the next cols
        data_used <- c(data_used, new_col)
        result[[names(new_col)]] <- new_col[[1]]
      }
      self$data[[data]] <- tibble::as_tibble(result)
      return(self)
    },


    # -------------------------------------------------------------------------
    # Sort columns of tables in obj$data
    arrange_obs = function(data, ..., target = NULL) {

      # Check for use of "target"
      if (! is.null(target)) {
        stop('Use of "target" is depreciated. Use "data" instead.')
      }

      # Parse data option
      data <- parse_dataset(self, data)

      # Check that multiple datasets are the same length
      dataset_length <- vapply(self$data[data], FUN.VALUE = numeric(1),
                               function(one) {
                                 if (is.data.frame(one)) {
                                   return(nrow(one))
                                 } else {
                                   return(length(one))
                                 }
                               })
      dataset_length <- unique(dataset_length)
      if (length(dataset_length) > 1) {
        stop(call. = FALSE,
             "If multiple datasets are filtered at once, then they must the same length. ",
             "The following lengths were found for the specified datasets:\n",
             limited_print(type = "silent", prefix = "  ", dataset_length))
      }

      # Sort observations
      data_used <- self$data_used(...)
      for (one in data) {
        if (is.data.frame(self$data[[one]])) { # if it is a table
          sort_cols <- data_used[! names(data_used) %in% names(self$data[[one]])]
          if (length(sort_cols) == 0) {
            self$data[[one]] <- dplyr::arrange(self$data[[one]], ...)
          } else {
            target_with_extra_cols <-
              dplyr::bind_cols(sort_cols, self$data[[one]])
            self$data[[one]] <-
              dplyr::arrange(target_with_extra_cols, ...)[, -seq_along(sort_cols)]
          }
        } else { # if it is a list or vector
          dummy_table <- data.frame(index = seq_along(self$data[[one]]))
          if (length(data_used)!= 0) {
            dummy_table <- dplyr::bind_cols(data_used, dummy_table)
          }
          dummy_table <- dplyr::arrange(dummy_table, ...)
          self$data[[one]] <- self$data[[one]][dummy_table$index]
        }
      }

      return(self)
    },


    # -------------------------------------------------------------------------
    # Randomly sample some number of observations from a table
    sample_n_obs = function(data, size, replace = FALSE, taxon_weight = NULL,
                            obs_weight = NULL, use_supertaxa = TRUE,
                            collapse_func = mean, ..., target = NULL) {

      # Check for use of "target"
      if (! is.null(target)) {
        warning(call. = FALSE,
                'Use of "target" is depreciated. Use "data" instead.')
        data <- target
      }

      # Parse data option
      data <- parse_dataset(self, data)

      # non-standard argument evaluation
      data_used <- eval(substitute(self$data_used(taxon_weight, obs_weight)))
      taxon_weight <- lazyeval::lazy_eval(lazyeval::lazy(taxon_weight),
                                          data = data_used)
      obs_weight <- lazyeval::lazy_eval(lazyeval::lazy(obs_weight),
                                        data = data_used)

      # Check that multiple datasets are the same length
      dataset_length <- vapply(self$data[data], FUN.VALUE = numeric(1),
                               function(one) {
                                 if (is.data.frame(one)) {
                                   return(nrow(one))
                                 } else {
                                   return(length(one))
                                 }
                               })
      dataset_length <- unique(dataset_length)
      if (length(dataset_length) > 1) {
        stop(call. = FALSE,
             "If multiple datasets are sampled at once, then they must the same length. ",
             "The following lengths were found for the specified datasets:\n",
             limited_print(type = "silent", prefix = "  ", dataset_length))
      }

      # Calculate taxon component of taxon weights
      if (is.null(taxon_weight)) {
        obs_taxon_weight <- rep(1, dataset_length)
      } else {
        obs_index <- match(self$get_data_taxon_ids(data, require = TRUE),
                           self$taxon_ids())
        my_supertaxa <- self$supertaxa(recursive = use_supertaxa,
                                       simplify = FALSE, include_input = TRUE,
                                       na = FALSE,
                                       value = "taxon_indexes")
        taxon_weight_product <- vapply(
          my_supertaxa,
          function(x) collapse_func(taxon_weight[x]),
          numeric(1)
        )
        obs_taxon_weight <- taxon_weight_product[obs_index]
      }
      obs_taxon_weight <- obs_taxon_weight / sum(obs_taxon_weight)

      # Calculate observation component of observation weights
      if (is.null(obs_weight)) {
        obs_weight <- rep(1, dataset_length)
      }
      obs_weight <- obs_weight / sum(obs_weight)

      # Combine observation and taxon weight components
      combine_func <- prod
      weight <- mapply(obs_taxon_weight, obs_weight,
                       FUN = function(x, y) combine_func(c(x,y)))
      weight <- weight / sum(weight)

      # Sample observations
      sampled_rows <- sample.int(dataset_length, size = size,
                                 replace = replace, prob = weight)
      self$filter_obs(data, sampled_rows, ...)
    },

    # -------------------------------------------------------------------------
    # Randomly sample some proportion of observations from a table
    sample_frac_obs = function(data, size, replace = FALSE,
                               taxon_weight = NULL, obs_weight = NULL,
                               use_supertaxa = TRUE,
                               collapse_func = mean, ..., target = NULL) {

      # Check for use of "target"
      if (! is.null(target)) {
        warning(call. = FALSE,
                'Use of "target" is depreciated. Use "data" instead.')
        data <- target
      }

      # Parse data option
      data <- parse_dataset(self, data)

      # Check that multiple datasets are the same length
      dataset_length <- vapply(self$data[data], FUN.VALUE = numeric(1),
                               function(one) {
                                 if (is.data.frame(one)) {
                                   return(nrow(one))
                                 } else {
                                   return(length(one))
                                 }
                               })
      dataset_length <- unique(dataset_length)
      if (length(dataset_length) > 1) {
        stop(call. = FALSE,
             "If multiple datasets are sampled at once, then they must the same length. ",
             "The following lengths were found for the specified datasets:\n",
             limited_print(type = "silent", prefix = "  ", dataset_length))
      }

      # Call sample_n_obs
      eval(substitute(self$sample_n_obs(data = data,
                                        size = size * dataset_length,
                                        replace = replace,
                                        taxon_weight = taxon_weight,
                                        obs_weight = obs_weight,
                                        use_supertaxa = use_supertaxa,
                                        collapse_func = collapse_func,
                                        ...)))
    },


    # -------------------------------------------------------------------------
    # Count observations for each taxon in a data set
    n_obs = function(data = NULL, target = NULL) {

      # Check for use of "target"
      if (! is.null(target)) {
        warning(call. = FALSE,
                'Use of "target" is depreciated. Use "data" instead.')
        data <- target
      }

      # If no data is specified, use the first dataset
      if (is.null(data)) {
        if (length(self$data) > 0) {
          data <- names(self$data)[1]
        } else {
          stop(paste0('There are no data sets to get observation info from.'))
        }
      }

      # Parse data option
      data <- self$get_data(data)

      # Count observations
      vapply(self$obs(data, recursive = TRUE, simplify = FALSE),
             length, numeric(1))
    },

    # -------------------------------------------------------------------------
    # Count observations for each taxon in a data set, including observations
    # for the specific taxon but NOT the observations of its subtaxa.
    n_obs_1 = function(data = NULL, target = NULL) {

      # Check for use of "target"
      if (! is.null(target)) {
        warning(call. = FALSE,
                'Use of "target" is depreciated. Use "data" instead.')
        data <- target
      }

      if (is.null(data)) {
        if (length(self$data) > 0) {
          data <- names(self$data)[1]
        } else {
          stop(paste0('There are no data sets to get observation info from.'))
        }
      }
      vapply(self$obs(data, recursive = FALSE, simplify = FALSE),
             length, numeric(1))
    },

    # Find taxon ids for datasets by dataset name
    #
    # require: if TRUE, require that taxon ids be present, or make an error
    get_data_taxon_ids = function(dataset_name, require = FALSE, warn = FALSE, message = FALSE) {

      stop_or_warn <- function(text) {
        if (require) {
          stop(call. = FALSE, text)
        }
        if (warn) {
          warning(call. = FALSE, text)
        } else if (message) {
          message(text)
        }
      }

      # Get the dataset
      if (length(dataset_name) == 1 && # data is name/index of dataset in object
          (dataset_name %in% names(self$data) || is.numeric(dataset_name))) {
        data <- self$data[[dataset_name]]
      } else { # it is an external data set, not in the object
        data <- dataset_name
        dataset_name <- deparse(substitute(dataset_name))
      }

      # Extract taxon ids if they exist
      if (is.data.frame(data)) {
        if ("taxon_id" %in% colnames(data)) {
          is_valid <- private$ids_are_valid(data$taxon_id)
          if (all(is_valid)) {
            return(data$taxon_id)
          } else  {
            stop_or_warn(paste0('There is a "taxon_id" column in the data set "',
                                dataset_name, '", but the following invalid IDs were found:\n  ',
                                limited_print(data$taxon_id[! is_valid], type = "silent")))
            return(NULL)
          }
        } else {
          stop_or_warn(paste0('There is no "taxon_id" column in the data set "',
                              dataset_name, '", so there are no taxon IDs.'))
          return(NULL)
        }
      } else if (class(data) == "list" || is.vector(data) || can_be_used_in_taxmap(data)) {
        if (! is.null(names(data))) {
          is_valid <- private$ids_are_valid(names(data))
          if (all(is_valid)) {
            return(names(data))
          } else if (all(! is_valid)) {
            stop_or_warn(paste0('The data set "', dataset_name,
                                '" is named, but not named by taxon ids.'))
            return(NULL)
          } else { # some are valid, but not all
            stop_or_warn(paste0('Dataset "', dataset_name, '" appears to be named by taxon IDs, but contains ', sum(! is_valid), ' invalid IDs:\n  ',
                                limited_print(names(data)[! is_valid], type = "silent")))
            return(NULL)
          }
        } else {
          stop_or_warn(paste0('The data set "', dataset_name,
                              '" is unnamed, ',
                              'so there are no taxon ids.'))
          return(NULL)
        }
      } else {
        stop_or_warn(paste0('I dont know how to extract taxon ids from dataset "', dataset_name,
                            '" of type "', class(data)[1], '".'))
        return(NULL)
      }
    },

    # Get a data set from a taxmap selfect
    get_dataset = function(data) {

      # Convert logicals to numerics
      if (is.logical(data)) {
        if (length(data) != length(self$data)) {
          stop("When using a TRUE/FALSE vector to specify the data set, it must be the same length as the number of data sets",
               call. = FALSE)
        } else {
          data <- which(data)
        }
      }

      # Check for multiple/no values
      if (length(data) == 0) {
        stop('No data specified.', call. = FALSE)
      }
      if (length(data) > 1) {
        stop('Only one dataset can be used.', call. = FALSE)
      }

      # Check that dataset exists
      error_msg <- paste0('The dataset "', data,
                          '" cannot be found. Datasets found include:\n  ',
                          limited_print(paste0("[", seq_along(self$data), "] ", names(self$data)),
                                        type = "silent"))
      if (is.character(data)) {
        if (! data %in% names(self$data)) {
          stop(error_msg, call. = FALSE)
        }
      } else if (is.numeric(data)) {
        if (! data %in% seq_along(self$data)) {
          stop(error_msg, call. = FALSE)
        }
      }

      # Return without printing
      return(self$data[[data]])
    }

  ),

  private = list(
    nse_accessible_funcs = c(
      "taxon_names",
      "taxon_ids",
      "taxon_indexes",
      "classifications",
      "n_supertaxa",
      "n_supertaxa_1",
      "n_subtaxa",
      "n_subtaxa_1",
      "n_leaves",
      "n_leaves_1",
      "taxon_ranks",
      "is_root",
      "is_stem",
      "is_branch",
      "is_leaf",
      "is_internode",
      "n_obs",
      "n_obs_1"
    ),

    # Remove observations from a particular dataset or just remove the taxon ids
    # NOTE: indexes = what is NOT removed
    remove_obs = function(data, indexes, unname_only = FALSE) {
      if (unname_only) {
        if (is.data.frame(self$data[[data]])) {
          self$data[[data]][! indexes, "taxon_id"] <- as.character(NA)
        } else {
          names(self$data[[data]])[! indexes] <- as.character(NA)
        }
      } else {
        if (is.data.frame(self$data[[data]])) {
          self$data[[data]] <-
            self$data[[data]][indexes, , drop = FALSE]
        } else {
          self$data[[data]] <- self$data[[data]][indexes]
        }
      }
    },

    # Checks if a character vector contains only taxon IDs.
    # Returns logical vector same length as input
    ids_are_valid = function(ids_to_check) {
      ids_to_check %in% c(self$taxon_ids(), NA_character_)
    }
  )
)
