#' Convert `data` input for Taxamp
#'
#' Make sure `data` is in the right format and complain if it is not. Then, add
#' a `taxon_id` column to data with the same length as the input
#'
#' @param self The newly created [taxmap()] object
#' @param data The `data` variable passed to the `Taxmap` constructor
#' @param input_ids The taxon IDs for the inputs that made the taxonomy
#' @param assume_equal If `TRUE`, and a data set length is the same as the
#'   `input_ids` length, then assume that `input_ids` applies to the data set as
#'   well.
#'
#' @return A `data` variable with the right format
#'
#' @keywords internal
init_taxmap_data <- function(self, data, input_ids, assume_equal = TRUE) {

  process_one <- function(x, name) {

    # Check for zero-length inputs
    if (length_of_thing(x) <= 0) {
      return(x)
    }

    if (is.data.frame(x)) {
      # Convert all data.frames to tibbles
      if  (! tibble::is_tibble(x)) {
        x <- tibble::as_tibble(x)
      }

      # Add the `taxon_id` column if it is not already there
      if ("taxon_id" %in% colnames(x)) {
        is_valid <- self$is_taxon_id(x$taxon_id)
        if (all(is_valid)) {
          message(paste0('Using existing "taxon_id" column for table "',
                         name, '"'))
        } else { # has a "taxon_id" column but invalid IDs
          stop(call. = FALSE,
               paste0('The table "', name,
                      '" has a "taxon_id" column, but the values do not appear to be taxon IDs.',
                      'The following ', sum(! is_valid), ' of ',
                      length(is_valid), ' values are not valid indexes:\n',
                      limited_print(x$taxon_id[! is_valid], type = "silent")))
        }
      } else if ("taxon_index" %in% colnames(x)) {
        message(paste0('Using "taxon_index" column to assign taxon IDs for table "',
                       name, '".'))
        x <- dplyr::bind_cols(taxon_id = unname(input_ids[as.integer(x$taxon_index)]), x)
      } else if (assume_equal && nrow(x) == length(input_ids)) {
        message(paste0('Assuming that the elements of table "', name,
                       '" are in the same order as taxon information.'))
        x <- dplyr::bind_cols(taxon_id = unname(input_ids), x)
      } else {
        warning(call. = FALSE,
                paste0('The table "', name,
                       '" does not have a "taxon_index" column, "taxon_id" column, or a number of',
                       'rows equal to the number of inputs, so no "taxon_id"',
                       ' can be assigned.'))
      }
    } else if (is.vector(x) || is.list(x)) {
      # Add the `taxon_id` column if it is not already there
      if (is.null(names(x))) { # data set is unnamed
        if (assume_equal && length(x) == length(input_ids)) { # if data length is same as input assume they match up
          names(x) <- input_ids
          message(paste0('Assuming that the elements of list/vector "', name,
                         '" are in the same order as taxon information.'))
        } else { # No taxonomy information
          warning(call. = FALSE,
                  paste0('The list/vector "', name,
                         '" is unnamed so has no taxon ID information.'))
        }
      } else { # data set has names
        is_valid <- self$is_taxon_id(names(x))
        if (all(is_valid)) { # All are valid taxon ids
          message(paste0('Using existing names of list/vector "', name,
                         '" as taxon IDs.'))
        } else { # data set has names, but they are not valid IDs
          warning(call. = FALSE,
                  paste0('The list/vector "', name,
                         '" is named, but the names do not appear to be taxon IDs.',
                         'The following ', sum(! is_valid), ' of ',
                         length(is_valid), ' names are not valid indexes:\n',
                         limited_print(names(x)[! is_valid], type = "silent")))
        }
      }
    }
    return(x)
  }

  # Get names of data inputs for messages
  data_names <- names(data)
  if (is.null(data_names)) {
    data_names <- rep(NA, length(data))
  }
  data_names <- ifelse(is.na(data_names) | data_names == "",
                       paste0("input_", seq_along(data)),
                       data_names)

  # Process each input
  mapply(process_one, data, data_names, SIMPLIFY = FALSE)
}


#' Validate `funcs` input for Taxamp
#'
#' Make sure `funcs` is in the right format and complain if it is not.
#' NOTE: This currently does nothing.
#'
#' @param funcs The `funcs` variable passed to the `Taxmap` constructor
#'
#' @return A `funcs` variable with the right format
#'
#' @keywords internal
validate_taxmap_funcs <- function(funcs) {
  funcs
}


#' used to parse inputs to `drop_obs` and `reassign_obs`
#'
#' @keywords internal
parse_possibly_named_logical <- function(input, data, default) {
  if (is.null(names(input))) {
    if (length(input) == 1) {
      output <- stats::setNames(rep(input, length(data)),
                                names(data))
    } else if (length(input) == length(data)) {
      output <- stats::setNames(input, names(data))
    } else {
      stop(paste("Invalid input for logical vector selecting which data",
                 "sets to affect. Valid inputs include:\n",
                 "1) a single unnamed logical (e.g. TRUE)\n",
                 "2) one or more named logicals with names matching",
                 "data sets in obj$data (e.g. c(data_1 = TRUE, data_2",
                 "= FALSE)\n  3) an unamed logical vector of the same",
                 "length as obj$data."))
    }
  } else {
    if (length(not_data_names <-
               names(input)[! names(input) %in% names(data)]) > 0) {
      stop(paste0("Invalid input for logical vector selecting which data",
                  " sets to affect. The following names are not in",
                  " data: ",
                  paste0(not_data_names, collapse = ", ")))
    }
    output <- stats::setNames(rep(default, length(data)),
                              names(data))
    output[names(input)] <- input
  }
  return(output)
}

#' Get list of usable functions
#'
#' Returns the names of all functions that can be called from any environment
#'
#' @return vector
#'
#' @keywords internal
all_functions <- function() {
  objects <- unlist(lapply(seq_len(length(search())), function(i) ls(pos = i)))
  is_func <- vapply(objects, function(obj) is.function(get(obj)), logical(1))
  return(objects[is_func])
}


#' Check dataset format
#'
#' Check that the datasets in a [taxmap()] object are in the correct format.
#' * Checks that column names are not the names of functions
#'
#' @param obj A [taxmap()] object
#'
#' @return NULL
#'
#' @keywords internal
check_taxmap_data <- function(obj) {
  #  Check that column names are not the names of functions
  # data_names <- all_names(obj, funcs = FALSE, builtin_funcs = FALSE)
  # suspect_names <- data_names[data_names %in% all_functions()]
  # if (length(suspect_names) > 0) {
  #   warning(paste0("Naming table columns/vectors/lists the same name as ",
  #                  "functions can sometimes interfere with non-standard ",
  #                  "evaluation. The following data shares names with ",
  #                  "functions:\n", limited_print(names(suspect_names),
  #                                                type = "silent")),
  #           call. = FALSE)
  # }

  return(invisible(NULL))
}


#' Convert a vector to database IDs
#'
#' This is a convenience function to convert to identifiers of various  data
#' sources. It wraps the \code{as.*id} functions in \code{\link[taxize]{taxize}}
#'
#' @param ids The character or numeric vector of raw taxon IDs.
#' @param database The database format to convert the IDs to. Either ncbi,
#'   itis, eol, col, tropicos, gbif, nbn, worms, natserv, bold, or wiki
#' @param ... Passed to \code{as.*id} function.
#'
#' @keywords internal
as_id <- function(ids, database, ...) {
  id_constructors <- list(ncbi = taxize::as.uid,
                          itis = taxize::as.tsn,
                          iucn = taxize::as.iucn,
                          eol = taxize::as.eolid,
                          col = taxize::as.colid,
                          tropicos = taxize::as.tpsid,
                          gbif = taxize::as.gbifid,
                          nbn = taxize::as.nbnid,
                          worms = taxize::as.wormsid,
                          natserv = taxize::as.natservid,
                          bold = taxize::as.boldid,
                          wiki = taxize::as.wiki)

  if (! database %in% names(id_constructors)) {
    stop(call. = FALSE,
         paste0('"', database,
                '" is not a recognized database name. Must be one of the following:\n',
                limited_print(names(id_constructors), type = "silent")))
  }

  id_constructors[[database]](ids, ...)
}


#' lappy with progress bars
#'
#' Immitates lapply with optional progress bars
#'
#' @param X The thing to iterate over
#' @param FUN The function to apply to each element
#' @param progress (logical of length 1) Whether or not to print a progress bar. Default is to only print a progress bar during interactive use.
#' @param ... Passed to function
#'
#' @return list
#' @keywords internal
progress_lapply <- function(X, FUN, progress = interactive(), ...) {
  if (progress) {
    progress_bar <- utils::txtProgressBar(min = 0, max = length(X), style = 3)
    one_iteration <- function(index) {
      output <- FUN(X[[index]], ...)
      utils::setTxtProgressBar(progress_bar, index)
      return(output)
    }
    output <- lapply(seq_len(length(X)), one_iteration)
    close(progress_bar)
  } else {
    output <- lapply(X, FUN, ...)
  }

  return(output)
}


#' Check that a unknown object can be used with taxmap
#'
#' Check that a unknown object can be assigned taxon IDs and filtered.
#'
#' @param obj
#'
#' @return TRUE/FALSE
#'
#' @keywords internal
can_be_used_in_taxmap <- function(obj) {
  # Check if is one-dimension, with a known length (i.e. works with length function)
  obj_length <- tryCatch({length(obj)}, error = function(e) return(NA))
  if (is.na(obj_length)) {
    return(FALSE)
  }

  # Can be subset (i.e. [ changes the value of length)
  if (obj_length > 0) {
    new_length <- tryCatch({length(obj[numeric(0)])}, error = function(e) return(NA))
    if (is.na(new_length) || new_length != 0) {
      return(FALSE)
    }
    new_length <- tryCatch({length(obj[c(1, 1)])}, error = function(e) return(NA))
    if (is.na(new_length) || new_length != 2) {
      return(FALSE)
    }
  }

  # Has per-element names that can be set (i.e. works with names and names<- function)
  my_names <- tryCatch({names(obj)}, error = function(e) return(NA))
  if (length(my_names) == 1 && is.na(my_names)) {
    return(FALSE)
  }
  if (obj_length > 0) {
    new_names <- tryCatch({
      names(obj)[1] <- "name_test"
      names(obj)[1]
    },
    error = function(e) return(NA))
    if (new_names != "name_test") {
      return(FALSE)
    }
  }

  return(TRUE)
}
