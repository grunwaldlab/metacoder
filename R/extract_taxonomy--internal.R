#' Report a error/warning if needed
#' 
#' Report a error/warning if needed
#' 
#' @param text The error to report
#' @param vigilance (\code{character} of length 1) Controls the reporting of possible problems.
#' The following values are possible: 
#'  \describe{
#'    \item{\code{"none"}}{No warnings or errors are generated if the function can complete.}
#'    \item{\code{"message"}}{A message is generated when atypical events occur.}
#'    \item{\code{"warning"}}{Warnings are generated when atypical events occur.}
#'    \item{\code{"error"}}{Errors are generated when atypical events occur, stopping the completion of the function.}
#'  } 
#' 
#' @return \code{NULL}
#' 
#' @keywords internal
vigilant_report <- function(text, vigilance = c("error", "none", "message", "warning")) {
  vigilance <- match.arg(vigilance)
  text <- paste0(text, "\n",
                "To avoid this ", vigilance, ", change the setting of the `vigilance` option")
  response <- list("error" = stop, "warning" = warning, "message" = message,
                   "none" = function(text) invisible(NULL))
  response[[vigilance]](text)
}


#' Check that all match input
#' 
#' Ensure that all of a character vector matches a regex.
#' Inputs that do not match are excluded.
#' 
#' @param input (\code{character})
#' @param regex (\code{character} of length 1)
#' @param max_print  (\code{numeric} of length 1)
#' The maximum number of lines to print in error/warning messages.
#' @param ... passed to \code{\link{vigilant_report}}
#'  
#' @return \code{character} Parts of \code{input} matching \code{regex}
#' 
#' @keywords internal
validate_regex_match <- function(input, regex, max_print = 10, ...) {
  # check which input values match
  input <- as.character(input)
  not_matching <- ! grepl(pattern = regex, x = input)
  # complain about those that dont
  if (sum(not_matching) > 0) {
    invalid_list <- paste("   ", which(not_matching), ": ", input[not_matching], "\n")
    if (length(invalid_list) > max_print) {
      invalid_list <- c(invalid_list[1:max_print], "    ...")
    }
    vigilant_report(...,
                    paste0(collapse = "",
                           c("The following ", sum(not_matching), " of ", length(input),
                             " input(s) could not be matched by the regex supplied:\n",
                             invalid_list)))
  }
  # return matching inputs
  return(input[! not_matching])
}




#' Check a regex-key pair
#' 
#' Checks that the number of capture groups in the regex matches the length of the key.
#' Checks that only certian values of \code{key} can appear more that once.
#' Adds names to keys that will be used for column names in the output of \code{extract_taxonomy}.
#' Uses non-standard evaluation to get the name of input variables.
#'
#' @param regex (\code{character})
#' A regex with capture groups
#' @param key (\code{character})
#' A key corresponding to \code{regex}
#' @param multiple_allowed (\code{character})
#' Values of \code{key_options} that can appear more than once.
#' 
#' @return Returns the result of \code{\link{match.arg}} on the key.
#' 
#' @keywords internal
#' @rdname validate_regex_key_pair
validate_regex_key_pair_ <- function(regex, key, multiple_allowed) {
  key_exp <- deparse(key$expr)
  key_value <- lazyeval::lazy_eval(key)
  regex_exp <- deparse(regex$expr)
  regex_value <- lazyeval::lazy_eval(regex)
  
  # Check key length
  regex_capture_group_count <- count_capture_groups(regex_value)
  key_length <- length(key_value)
  if (key_length != regex_capture_group_count) {
    stop(paste0(collapse = "",
                'The number of capture groups in `', regex_exp, '` and the length of "', 
                key_exp, '" do not match.\n',
                'The key has ', key_length, ' term(s) and the regex has ', regex_capture_group_count))
  }
  
  # Check that only names in `multiple_allowed` appear more than once
  counts <- table(key_value)
  duplicated_keys <- names(counts[counts > 1])
  if (! all(duplicated_keys %in% multiple_allowed)) {
    stop(paste0(collapse = "",
                'The following values in `', key_exp, '` have been used more than once:\n',
                paste0(collapse ="", "    ", duplicated_keys, "\n"),
                '  Only the following keys can be duplicated:\n',
                paste0(collapse ="", "    ", multiple_allowed, "\n")))
  }
  
  # Apply key name defaults
  key_names <- names(key_value)
  if (is.null(key_names)) { key_names <- rep("", length(key_value)) }
  key_names[key_names == ""] <- key_value[key_names == ""]
  
  # Number key names that appear more than once
  for (a_key in duplicated_keys) {
    key_names[key_names == a_key] <- paste0(key_names[key_names == a_key], "_",
                                            seq_along(key_names[key_names == a_key]))
  }
  
  names(key_value) <- key_names
  return(key_value)
}

#' @keywords internal
#' @rdname validate_regex_key_pair
validate_regex_key_pair <- function(regex, key, multiple_allowed) {
  validate_regex_key_pair_(lazyeval::lazy(regex), lazyeval::lazy(key), multiple_allowed)
}


#' Count capture groups
#' 
#' Count the number of capture groups in a regular expression.
#' 
#' @param regex (\code{character} of length 1)
#' 
#' @return \code{numeric} of length 1
#' 
#' @source http://stackoverflow.com/questions/16046620/regex-to-count-the-number-of-capturing-groups-in-a-regex
#' 
#' @keywords internal
count_capture_groups <- function(regex) {
  new_regex <- paste0(regex, "|") # Allow the regex to match nothing
  ncol(stringr::str_match(string = "", pattern = new_regex)) - 1
}



#' List of classifications to taxonomy tree
#' 
#' Convert a list of classifications to a taxonomic tree structure.
#' Columns in the input not used to create the taxonomy are preserved in the output.
#' Unique IDs can be assigned if the column used to create the taxonomy does not
#' have unique values for place in the hierarchy (e.g. same species epithet in different phlya.)
#' Missing data does not effect the taxonomy, but its positon is preserved in the item output. 
#' 
#' @param classifications (\code{list} of \code{data.frame}) 
#' The classifications of a set of items.
#' Not necessarily unique.
#' Rows should correspond to taxa and columns to information associated with those taxa.
#' @param id_column (\code{character} of length 1)
#' The name of the column present in each \code{data.frame} in \code{classifications} that will
#' be used to infer the taxonomic tree structure.
#' @param make_ids (\code{logical} of length 1) If \code{TRUE}, make unique IDs based off of the
#' contents of \code{id_column}.
#' 
#' @return An object of type \code{classified}
#' 
#' @keywords internal
class_to_taxonomy <- function(classifications, id_column, make_ids = FALSE) {
  
  recursive_part <- function(group, level = 1, parent = NA) {
    # Split list of classifications based on level
    split_group <- split_class_list(group, level, id_column)
    # Make rows for current taxon-parent relationships
    taxon_rows <- lapply(split_group, function(x) c(x[[1]][level, ], my_parent_ = parent))
    # Add taxon index to output
    taxon_index <- row_count + seq_along(taxon_rows)
    taxon_rows <- mapply(function(my_row, index) c(my_row, my_taxon_id_ = index),
                         SIMPLIFY = FALSE, taxon_rows, taxon_index)
    # Increment the number of rows in the output
    row_count <<- row_count + length(taxon_index)
    # Run this function on each of the subgroups
    child_taxa <- mapply(recursive_part, SIMPLIFY = FALSE,
                         group = child_groups, level = level + 1, parent = taxon_index)
    # Return the result of this instance of the function and the ones it makes
    return(c(taxon_row, child_taxa))
  }
  
  # make index counter to be used inside the recursive part
  row_count <- 0
  
  # Run recursive part of the function
  taxonomy <- do.call(rbind, recursive_part(classifications))
  
  # Format the output into a classified object
  
}



#' Split a list of classifications by a row/column value
#' 
#' Split a list of classifications into a list of list by unique values of a
#' specific row/column.
#' 
#' @param classifications (\code{list} of \code{data.frame}) 
#' The classifications of a set of items.
#' Not necessarily unique.
#' Rows should correspond to taxa and columns to information associated with those taxa.
#' @param row_index
#' The row in each \code{data.frame}
#' @param col_index
#' The column in each \code{data.frame}
#' 
#' @return \code{list} of \code{list} of \code{data.frame}
#' 
#' @keywords internal 
split_class_list <- function(classifications, row_index, col_index)  {
  split_by_values <- unlist(lapply(classifications, function(x) x[row_index, col_index]))
  split(classifications, split_by_values)
}