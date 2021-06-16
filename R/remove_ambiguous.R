#' Get patterns for ambiguous taxa
#'
#' This function stores the regex patterns for ambiguous taxa.
#'
#' @param unknown If \code{TRUE}, include names that suggest they are
#'   placeholders for unknown taxa (e.g. "unknown ...").
#' @param uncultured If \code{TRUE}, include names that suggest they are
#'   assigned to uncultured organisms (e.g. "uncultured ...").
#' @param regex If \code{TRUE}, includes regex syntax to make matching things like spaces more robust.
#' @param case_variations If \code{TRUE}, include variations of letter case. 
#'
#' @export
ambiguous_synonyms <- function(unknown = TRUE, uncultured = TRUE, regex = TRUE, case_variations = FALSE) {
  unknown_syns <- c(
    'unknown',
    'unidentified',
    'incertae sedis',
    'ambiguous',
    'ambiguous taxa',
    'unassigned',
    'possible',
    'putative'
  )
  uncultured_syns <- c(
    'uncultured',
    'candidatus',
    'metagenome'
  )
  output <- c()
  if (unknown) {
    output <- c(output, unknown_syns)
  }
  if (uncultured) {
    output <- c(output, uncultured_syns)
  }
  if (case_variations) {
    output <- c(output,
                capitalize(output),
                toupper(output))
  }
  if (regex) {
    output <- gsub(output, pattern = ' ', replacement = '[_ -]+')
  }
  return(output)
}



#' Get patterns for ambiguous taxa
#'
#' This function stores the regex patterns for ambiguous taxa.
#'
#' @param unknown If \code{TRUE}, Remove taxa with names the suggest they are
#'   placeholders for unknown taxa (e.g. "unknown ...").
#' @param uncultured If \code{TRUE}, Remove taxa with names the suggest they are
#'   assigned to uncultured organisms (e.g. "uncultured ...").
#' @param case_variations If \code{TRUE}, include variations of letter case. 
#' @param whole_match If \code{TRUE}, add "^" to front and "$" to the back of each
#'   pattern to indicate they are to match whole words.
#' @param name_regex The regex code to match a valid character in a taxon name.
#'   For example, "[a-z]" would mean taxon names can only be lower case letters.
#'
#' @keywords internal
ambiguous_patterns <- function(unknown = TRUE, uncultured = TRUE, case_variations = FALSE,
                               whole_match = FALSE, name_regex = ".") {
  # Initialize output vector
  output <-  paste0(name_regex, "*",
                    ambiguous_synonyms(unknown = unknown,
                                       uncultured = uncultured,
                                       case_variations = case_variations),
                    name_regex, "*")
  
  # Add regex code for full matches
  if (whole_match) {
    output <- paste0("^", output, "$")
  }
  
  return(output)
}



#' Find ambiguous taxon names
#'
#' Find taxa with ambiguous names, such as "unknown" or "uncultured".
#' 
#' If you encounter a taxon name that represents an ambiguous taxon that is not
#' filtered out by this function, let us know and we will add it.
#'
#' @param taxon_names A \code{\link{taxmap}} object
#' @inheritParams ambiguous_patterns
#' @param ignore_case If \code{TRUE}, dont consider the case of the text when
#'   determining a match.
#'
#' @return TRUE/FALSE vector corresponding to \code{taxon_names}
#' 
#' @examples 
#' is_ambiguous(c("unknown", "uncultured", "homo sapiens", "kfdsjfdljsdf"))
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
#' \code{\link{filter_taxa}} with the "invert" option set to \code{TRUE}.
#' Works the same way as \code{\link{filter_taxa}} for the most part.
#' 
#' If you encounter a taxon name that represents an ambiguous taxon that is not
#' filtered out by this function, let us know and we will add it.
#'
#' @param obj A \code{\link{taxmap}} object
#' @inheritParams is_ambiguous
#' @inheritParams filter_taxa
#'
#' @return A \code{\link{taxmap}} object
#'
#' @examples 
#' obj <- parse_tax_data(c("Plantae;Solanaceae;Solanum;lycopersicum",
#'                         "Plantae;Solanaceae;Solanum;tuberosum",
#'                         "Plantae;Solanaceae;Solanum;unknown",
#'                         "Plantae;Solanaceae;Solanum;uncultured",
#'                         "Plantae;UNIDENTIFIED"))
#' filter_ambiguous_taxa(obj)
#'
#' @export
filter_ambiguous_taxa <- function(obj, unknown = TRUE, uncultured = TRUE,
                                  name_regex = ".", ignore_case = TRUE,
                                  subtaxa = FALSE, drop_obs = TRUE,
                                  reassign_obs = TRUE, reassign_taxa = TRUE) {
  # Identify taxa to filter out 
  to_remove <- is_ambiguous(obj$taxon_names(), unknown = unknown, 
                            uncultured = uncultured, name_regex = name_regex,
                            ignore_case = ignore_case)
  
  filter_taxa(obj, to_remove, invert = TRUE, 
                    subtaxa = subtaxa, drop_obs = drop_obs,
                    reassign_obs = reassign_obs, reassign_taxa = reassign_taxa)
}
