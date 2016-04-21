#' Check a regex-key pair
#' 
#' Check that the number of capture groups in the regex matches the length of the key
#' and that the key consists of valid options.
#' Uses non-standard evaluation to get the name of input variables.
#'
#' @param regex A regex with capture groups
#' @param key A key corresponding to \code{regex}
#' @param key_options Values \code{key} can take on.
#' 
#' @return Returns the result of \code{\link{match.arg}} on the key.
#' 
#' @keywords internal
#' @rdname validate_regex_key_pair
validate_regex_key_pair_ <- function(regex, key, key_options) {
  regex_capture_group_count <- count_capture_groups(lazy_eval(regex))
  key_length <- length(lazy_eval(key))
  if (key_length != regex_capture_group_count) {
    stop(paste0(collapse = "",
                'The number of capture groups in "', regex$expr, '"and the length of "', 
                key$expr, '" do not match.\n',
                'The key has ', key_length, ' term(s) and the regex has ', regex_capture_group_count))
  }
  match.arg(lazy_eval(key), key_options, several.ok = TRUE)
}

#' @keywords internal
#' @rdname validate_regex_key_pair
validate_regex_key_pair <- function(regex, key, key_options) {
  validate_regex_key_pair_(lazy(regex), lazy(key), key_options)
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