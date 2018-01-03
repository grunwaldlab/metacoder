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
                paste0(name_regex, "*", "unknown", name_regex, "*"),
                paste0(name_regex, "*", "Unknown", name_regex, "*"),
                paste0(name_regex, "*", "UNKNOWN", name_regex, "*"),
                paste0(name_regex, "*", "unidentified", name_regex, "*"),
                paste0(name_regex, "*", "Unidentified", name_regex, "*"),
                paste0(name_regex, "*", "UNIDENTIFIED", name_regex, "*"),
                paste0(name_regex, "*", "Incertae Sedis", name_regex, "*"),
                paste0(name_regex, "*", "incertae sedis", name_regex, "*"),
                paste0(name_regex, "*", "INCERTAE SEDIS", name_regex, "*"),
                paste0(name_regex, "*", "ambiguous", name_regex, "*"),
                paste0(name_regex, "*", "Ambiguous", name_regex, "*"),
                paste0(name_regex, "*", "AMBIGUOUS", name_regex, "*"),
                paste0(name_regex, "*", "unassigned", name_regex, "*"),
                paste0(name_regex, "*", "Unassigned", name_regex, "*"),
                paste0(name_regex, "*", "UNASSIGNED", name_regex, "*")
    )
  }
  
  # Add patterns for uncultured taxa Incertae Sedis
  if (uncultured) {
    output <- c(output, 
                paste0(name_regex, "*", "uncultured", name_regex, "*"),
                paste0(name_regex, "*", "Uncultured", name_regex, "*"),
                paste0(name_regex, "*", "UNCULTURED", name_regex, "*"), 
                paste0(name_regex, "*", "candidatus", name_regex, "*"),
                paste0(name_regex, "*", "Candidatus", name_regex, "*"),
                paste0(name_regex, "*", "CANDIDATUS", name_regex, "*")
    )
  }
  
  # Add regex code for full matches
  if (whole_match) {
    output <- paste0("^", output, "$")
  }
  
  return(output)
}



#' Find ambiguous taxon names
#'
#' Find taxa with ambiguous names, such as "unknown" or "uncultured". If you
#' encounter a taxon name that represents an ambiguous taxon that is not
#' filtered out by this function, let us know and we will add it.
#'
#' @param taxon_names A taxmap object
#' @inheritParams ambiguous_patterns
#' @param ignore_case If \code{TRUE}, dont consider the case of the text when
#'   determining a match.
#'
#' @return TRUE/FALSE vector corresponding to \code{taxon_names}
#'
#' @export
is_ambiguous <- function(taxon_names, unknown = TRUE, uncultured = TRUE,
                         name_regex = ".", ignore_case = TRUE) {
  # Get patterns to filter out
  patterns <- ambiguous_patterns(unknown = unknown, uncultured = uncultured,
                                 name_regex = name_regex)
  
  # Find which taxa to filter out 
  Reduce(`|`, lapply(patterns, function(x) {
    grepl(taxon_names, pattern = x, ignore.case = ignore_case)
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
                                  ignore_case = TRUE, subtaxa = FALSE, drop_obs = TRUE,
                                  reassign_obs = TRUE, reassign_taxa = TRUE) {
  # Identify taxa to filter out 
  to_remove <- is_ambiguous(obj$taxon_names(), unknown = unknown, 
                            uncultured = uncultured, name_regex = name_regex,
                            ignore_case = ignore_case)
  
  taxa::filter_taxa(obj, to_remove, invert = TRUE, 
                    subtaxa = subtaxa, drop_obs = drop_obs,
                    reassign_obs = reassign_obs, reassign_taxa = reassign_taxa)
}