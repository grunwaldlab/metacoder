#' Check that all match input
#' 
#' Ensure that all of a character vector matches a regex.
#' Inputs that do not match are excluded.
#' 
#' @param input (\code{character})
#' @param regex (\code{character} of length 1)
#' @param mismatch_action (\code{character} of length 1)
#' What to do if the regex does not match one or more of input.
#'  \describe{
#'    \item{\code{"allow"}}{Silently exclude mismatches.}
#'    \item{\code{"warn"}}{Like \code{"allow"} but issue a warning.}
#'    \item{\code{"error"}}{Throw an error if mismatches are found.}
#'  }
#' @param max_print  (\code{numeric} of length 1)
#' The maximum number of lines to print in error/warning messages.
#'  
#' @return \code{character} Parts of \code{input} matching \code{regex}
#' 
#' @keywords internal
validate_regex_match <- function(input, regex, mismatch_action = "error", max_print = 10) {
  mismatch_action <- match.arg(mismatch_action, c("allow", "warn", "error"))
  input <- as.character(input)
  not_matching <- ! grepl(pattern = regex, x = input)
  if (sum(not_matching) > 0) {
    invalid_list <- paste("   ", which(not_matching), ": ", input[not_matching], "\n")
    if (length(invalid_list) > max_print) {
      invalid_list <- c(invalid_list[1:max_print], "    ...")
    }
    error_text <- paste0(collapse = "",
                         c("The following ", sum(not_matching), " of ", length(input),
                           " input(s) could not be matched by the regex supplied:\n",
                           invalid_list))
    if (mismatch_action == "error") { stop(error_text) }
    if (mismatch_action == "warn") { warning(paste(sep = "\n", 
                                                   error_text, 
                                                   "They will be excluded.")) }
  }
  return(input[! not_matching])
}




#' Check a regex-key pair
#' 
#' Check that the number of capture groups in the regex matches the length of the key
#' and that the key consists of valid options.
#' Uses non-standard evaluation to get the name of input variables.
#'
#' @param regex (\code{character})
#' A regex with capture groups
#' @param key (\code{character})
#' A key corresponding to \code{regex}
#' @param key_options (\code{character})
#' Values \code{key} can take on.
#' @param multiple_allowed (\code{character})
#' Values of \code{key_options} that can appear more than once.
#' 
#' @return Returns the result of \code{\link{match.arg}} on the key.
#' 
#' @keywords internal
#' @rdname validate_regex_key_pair
validate_regex_key_pair_ <- function(regex, key, key_options, multiple_allowed = key_options) {
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
  
  # Check that key values are correct
  output_key <- match.arg(key_value, key_options, several.ok = TRUE)
  
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
  key_names[key_names == ""] <- output_key[key_names == ""]
  
  # Number key names that appear more than once
  for (a_key in duplicated_keys) {
    key_names[key_names == a_key] <- paste0(key_names[key_names == a_key], "_",
                                            seq_along(key_names[key_names == a_key]))
  }
  
  names(output_key) <- key_names
  return(output_key)
}

#' @keywords internal
#' @rdname validate_regex_key_pair
validate_regex_key_pair <- function(regex, key, key_options, multiple_allowed = key_options) {
  validate_regex_key_pair_(lazyeval::lazy(regex), lazyeval::lazy(key), key_options, multiple_allowed)
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
