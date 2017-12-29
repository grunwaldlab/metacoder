#' Get patterns for ambiguous taxa
#'
#' This function stores the regex patterns for ambiguous taxa.
#'
#' @param unknown If \code{TRUE}, Remove taxa with names the suggest they are
#'   placeholders for unknown taxa (e.g. "unknown ...").
#' @param uncultured If \code{TRUE}, Remove taxa with names the suggest they are
#'   assinged to uncultured organisms (e.g. "uncultured ...").
#' @param unknown If \code{TRUE}, add "^" to front and "$" to the back of each
#'   pattern to indicate they are to match whole words.
#' @param name_regex The regex code to match a valid character in a taxon name.
#'   For example, "[a-z]" would mean taxon names can only be lower case letters.
#'
#' @keywords internal
ambiguous_patterns <- function(unknown = TRUE, uncultured = TRUE,
                               whole_match = FALSE, name_regex = ".") {
  # Initialize output vector
  output <- c()
  
  # Add patterns for unknown taxa
  if (unknown) {
    output <- c(output, 
                paste0("unknown", name_regex, "+"),
                paste0("Unknown", name_regex, "+"),
                paste0("UNKNOWN", name_regex, "+"),
                paste0("unidentified", name_regex, "+"),
                paste0("Unidentified", name_regex, "+"),
                paste0("UNIDENTIFIED", name_regex, "+")
    )
  }
  
  # Add patterns for uncultured taxa
  if (uncultured) {
    output <- c(output, 
                paste0("uncultured", name_regex, "+"),
                paste0("Uncultured", name_regex, "+"),
                paste0("UNCULTURED", name_regex, "+"))
  }
  
  # Add regex code for full matches
  if (whole_match) {
    output <- paste0("^", output, "$")
  }
  
  return(output)
}



#' Find ambiguous taxon names
#'
#' Find taxa with ambiguous names, such as "unknown" or "uncultured". If
#' you encounter a taxon name that represents an ambiguous taxon that is not
#' filtered out by this function, let us know and we will add it.
#' 
#' @param taxon_names A taxmap object
#' @inheritParams ambiguous_patterns 
#' 
#' @return TRUE/FALSE vector corresponding to \code{taxon_names}
#' 
#' @export
is_ambiguous <- function(taxon_names, unknown = TRUE, uncultured = TRUE, name_regex = ".") {
  # Get patterns to filter out
  patterns <- ambiguous_patterns(unknown = unknown, uncultured = uncultured,
                                 name_regex = name_regex)
  
  # Find which taxa to filter out 
  Reduce(`|`, lapply(patterns, function(x) {
    grepl(taxon_names, pattern = x)
  }))
}


#' Filter ambiguous taxon names
#'
#' Filter out taxa with ambiguous names, such as "unknown" or "uncultured".
#' NOTE: some parameters of this function are passed to
#' \code{\link[taxa]{filter_taxa}} with the "invert" option set to \code{TRUE}.
#'
#' If you encounter a taxon name that represents an ambiguous taxon that is not
#' filtered out by this function, let us know and we will add it.
#'
#' @param obj A taxmap object
#' @inheritParams is_ambiguous
#' @inheritParams taxa::filter_taxa
#'
#' @return TRUE/FALSE vector corresponding to \code{taxon_names}
#'
#' @export
filter_ambiguous_taxa <- function(obj, unknown = TRUE, uncultured = TRUE, name_regex = ".", 
                                  subtaxa = FALSE, drop_obs = TRUE,
                                  reassign_obs = TRUE, reassign_taxa = TRUE) {
  # Identify taxa to filter out 
  to_remove <- is_ambiguous(obj$taxon_names(), unknown = unknown,
                            uncultured = uncultured, name_regex = name_regex)
  
  taxa::filter_taxa(obj, to_remove, invert = TRUE, 
                    subtaxa = subtaxa, drop_obs = drop_obs,
                    reassign_obs = reassign_obs, reassign_taxa = reassign_taxa)
}