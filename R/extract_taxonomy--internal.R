#' Report a error/warning if needed
#' 
#' Report a error/warning if needed
#' NOTE: This function is unusual in that it looks for a varaible names `vigilance` in a parent namespace.
#' 
#' @param text The error to report
#' 
#' @return \code{NULL}
#' 
#' @keywords internal
vigilant_report <- function(text) {
  vigilance <- dynGet("vigilance", ifnotfound = "error")
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
#'  
#' @return \code{character} Parts of \code{input} matching \code{regex}
#' 
#' @keywords internal
validate_regex_match <- function(input, regex, max_print = 10) {
  # check which input values match
  input <- as.character(input)
  not_matching <- ! grepl(pattern = regex, x = input)
  # complain about those that dont
  if (sum(not_matching) > 0) {
    invalid_list <- paste("   ", which(not_matching), ": ", input[not_matching], "\n")
    if (length(invalid_list) > max_print) {
      invalid_list <- c(invalid_list[1:max_print], "    ...")
    }
    vigilant_report(paste0(collapse = "",
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
validate_regex_key_pair <- function(regex, key, multiple_allowed) {
  
  # Non-standard evaluation
  regex_var_name <- deparse(substitute(regex))
  key_var_name <- deparse(substitute(key))
  
  # Check that the keys used are valid
  allowed <- eval(formals(extract_taxonomy.default)[[key_var_name]])
  invalid_keys <- key[! key %in% allowed]
  if (length(invalid_keys) > 0) {
    stop(paste0('Invalid key value "', invalid_keys[1], '" given.\n',
                'Valid options are: ', paste0(collapse = ", ", allowed)))
  }
  
  # Check key length
  regex_capture_group_count <- count_capture_groups(regex)
  key_length <- length(key)
  if (key_length != regex_capture_group_count) {
    stop(paste0(collapse = "",
                'The number of capture groups in "', regex_var_name, '" and the length of "', 
                key_var_name, '" do not match.\n',
                'The key has ', key_length, ' term(s) and the regex has ', regex_capture_group_count))
  }
  
  # Check that only keys in `multiple_allowed` appear more than once
  counts <- table(key)
  duplicated_keys <- names(counts[counts > 1])
  invalid_duplicates <- duplicated_keys[! duplicated_keys %in% multiple_allowed]
  if (length(invalid_duplicates) > 0) {
    stop(paste0(collapse = "",
                'The following values in `', key_var_name, '` have been used more than once: ',
                paste0(collapse =", ", invalid_duplicates), '\n',
                '  Only the following keys can be duplicated: ',
                paste0(collapse =", ", multiple_allowed)))
  }
  
  # Apply key name defaults
  key_names <- names(key)
  if (is.null(key_names)) { key_names <- rep("", length(key)) }
  key_names[key_names == ""] <- key[key_names == ""]
  
  names(key) <- key_names
  return(key)
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
#' @param item_data (\code{data.frame}) 
#' Data associated with the \code{classifications}.
#' 
#' @return An object of type \code{classified}
#' 
#' @keywords internal
class_to_taxonomy <- function(classifications, id_column, item_data = NULL) {
  
  recursive_part <- function(group, level = 1, parent = NA) {
    # Remove taxa that do not have inforamtion for the current level
    finished <- vapply(group, nrow, numeric(1)) < level
    group <- group[! finished]
    # Assign items to tip taxa
    item_taxa <- setNames(rep(parent, sum(finished)), names(finished[finished]))
    # Split list of classifications based on level
    split_group <- split_class_list(group, level, id_column)
    # Make rows for current taxon-parent relationships
    taxon_rows <- lapply(split_group, function(x) cbind(x[[1]][level, , drop = FALSE], my_parent_ = parent))
    # Get taxon index of parent
    taxon_index <- row_count + seq_along(taxon_rows)
    row_count <<- row_count + length(taxon_index)
    # Run this function on each of the subgroups
    if (length(split_group) > 0) {
      child_results <- mapply(recursive_part, SIMPLIFY = FALSE,
                           group = split_group, level = level + 1, parent = taxon_index)
      child_taxa <- unlist(lapply(child_results, function(x) x$taxon_data), recursive = FALSE, use.names = FALSE)
      child_items <- unlist(setNames(lapply(child_results, function(x) x$item_data), NULL))
    } else {
      child_taxa <- NULL
      child_items <- NULL
    }
    # Return the result of this instance of the function and the ones it makes
    return(list(taxon_data = c(taxon_rows, child_taxa), item_data = c(item_taxa, child_items)))
  }
  
  # make index counter to be used inside the recursive part
  row_count <- 0
  
  # Remove invalid classifications
  classifications <- classifications[!is.na(classifications)]
  
  # Run recursive part of the function
  names(classifications) <- seq_along(classifications)
  data <- recursive_part(classifications)
  taxonomy <- do.call(rbind, data$taxon_data)
  row.names(taxonomy) <- NULL
  item_index <- data$item_data[order(as.numeric(names(data$item_data)))]
  
  # Format the output into a classified object
  if (is.null(taxonomy))  {
    taxon_id <- character(0)
    parent_id <- integer(0)
    item_id <- rep(NA, nrow(item_data))
  } else {
    taxon_id <- 1:nrow(taxonomy) 
    parent_id <- taxonomy$my_parent_
    item_id <- unname(item_index)
  }
  classified(taxon_id = taxon_id,
             parent_id = parent_id,
             item_taxon_id = item_id,
             taxon_data = taxonomy[ , ! colnames(taxonomy) %in% c("my_parent_"), drop = FALSE],
             item_data = item_data)
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


#' Number duplicated names
#' 
#' Add a number to any character values that appear more that once.
#' 
#' @param my_names (\code{character})
#' @param replace_auto (\code{logical} of length 1)
#' If \code{TRUE}, undo the automated 'name.#' dataframe numbering
#' 
#' @return \code{character}
#' 
#' @keywords internal 
rename_duplicated <- function(my_names, replace_auto = TRUE) {
  # Undo automated numbering
  if (replace_auto) {
    my_names <- gsub(my_names, pattern = "\\.[0-9]+$", replacement = "")
  }
  
  # Number duplicates
  counts <- table(my_names)
  duplicated <- names(counts[counts > 1])
  for (x in duplicated) {
    my_names[my_names == x] <- paste0(my_names[my_names == x], "_",
                                              seq_along(my_names[my_names == x]))
  }
  return(my_names)
}


