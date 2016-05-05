#===================================================================================================
#' Extract taxonomy information from sequence headers
#' 
#' Extracts the taxonomy from metadata (e.g. sequence headers) or parsed sequence data. 
#' The location and identity of important information in the input is specified using a regular expression
#' with capture groups and an corresponding key.
#' An object of type \code{\link{classified}} is returned containing the specifed information.
#' Taxa are translated into unique codes if they are not already encoded this way.
#' 
#' @param input A vector from which to extract taxonomy information or an object of class
#' \code{\link{ape}{DNAbin}}. 
#' @param key (\code{character}) The identity of the capturing groups defined using \code{regex}.
#'  The length of \code{key} must be equal to the number of capturing groups specified in \code{regex}.
#'  Any names added to the terms will be used as column names in the output.
#'  Only \code{"taxon_info"} and \code{"item_info"} can be used multiple times.
#'  Each term must be one of those decribed below:
#'  \describe{
#'    \item{\code{taxon_id}}{A unique numeric id for a taxon for a particular \code{database} (e.g. ncbi accession number).
#'          Requires an internet connection.}
#'    \item{\code{name}}{The name of a taxon. Not necessarily unique, but are interpretable
#'          by a particular \code{database}. Requires an internet connection.}
#'    \item{\code{taxon_info}}{Arbitrary taxon info you want included in the output. Can be used more than once.}
#'    \item{\code{class}}{A list of taxa information that consitutes the full taxonomic classification
#'          from broad to specific (see \code{class_rev}) for a particular \code{database}. Individual taxa
#'          are separated by the \code{class_sep} argument and the information is parsed by the
#'          \code{class_regex} and \code{class_key} arguments.}
#'    \item{\code{item_id}}{An unique item (e.g. sequence) identifier for a particular \code{database}.
#'          Requires an internet connection.}
#'    \item{\code{item_info}}{Arbitrary item info you want included in the output. Can be used more than once.}
#'  }
#' @param regex (\code{character; length == 1}) A regular expression with capturing groups
#'  indicating the locations of relevant information. The identity of the information must
#'  be specified using the \code{key} argument.
#'  
#' @param class_key (\code{character} of length 1)
#' The identity of the capturing groups defined using \code{class_iregex}.
#' The length of \code{class_key} must be equal to the number of capturing groups specified in \code{class_regex}.
#' Any names added to the terms will be used as column names in the output.
#' At least \code{"taxon_id"} or \code{"name"} must be specified.
#' Only \code{"taxon_info"} can be used multiple times.
#' Each term must be one of those decribed below:
#'  \describe{
#'    \item{\code{taxon_id}}{A unique numeric id for a taxon for a particular \code{database} (e.g. ncbi accession number).
#'          Requires an internet connection.}
#'    \item{\code{name}}{The name of a taxon. Not necessarily unique, but are interpretable
#'          by a particular \code{database}. Requires an internet connection.}
#'    \item{\code{taxon_info}}{Arbitrary taxon info you want included in the output. Can be used more than once.}
#'  }
#' @param class_regex (\code{character} of length 1)
#' A regular expression with capturing groups indicating the locations of data for each taxon in the \code{class} term in the \code{key} argument.
#' The identity of the information must be specified using the \code{class_key} argument.
#' The \code{class_sep} option can be used to split the classification into data for each taxon before matching.
#' If \code{class_sep} is \code{NULL}, each match of \code{class_regex} defines a taxon in the classification. 
#' @param class_sep (\code{character} of length 1)
#' Used with the \code{class} term in the \code{key} argument.
#' The character(s) used to separate individual taxa within a classification.
#' After the string defined by the \code{class} capture group in \code{regex} is split by \code{class_sep},
#' its capture groups are extracted by \code{class_regex} and defined by \code{class_key}.
#' If \code{NULL}, every match of \code{class_regex} is used instead with first splitting by \code{class_sep}.
#' @param class_rev (\code{logical} of length 1)
#' Used with the \code{class} term in the \code{key} argument.
#' If \code{TRUE}, the order of taxon data in a classfication is reversed to be specific to broad.
#' 
#' @param database (\code{character} of length 1) The name of the database that patterns given in 
#'  \code{parser} will apply to. Valid databases include "ncbi", "itis", "eol", "col", "tropicos",
#'  "nbn", and "none". \code{"none"} will cause no database to be quired; use this if you want to not use the
#'  internet. NOTE: Only \code{"ncbi"} has been tested so far.
#' @param allow_na (\code{logical} of length 1) If \code{TRUE}, any missing data will be represented as \code{NA}s
#' in the output. This preserves the correspondance between the input and output values.
#' Missing data can be generated if the regex does not match the input or online queries fail.
#' @param vigilance (\code{character} of length 1) Controls the reporting of possible problems, such
#' as missing data and failed online queries (see \code{allow_na}).
#' The following values are possible: 
#'  \describe{
#'    \item{\code{"none"}}{No warnings or errors are generated if the function can complete.}
#'    \item{\code{"message"}}{A message is generated when atypical events occur.}
#'    \item{\code{"warning"}}{Warnings are generated when atypical events occur.}
#'    \item{\code{"error"}}{Errors are generated when atypical events occur, stopping the completion of the function.}
#'  } 
#' @param return_match (\code{logical} of length 1)
#' If \code{TRUE}, include the part of the input matched by \code{regex} in the output object.
#' @param return_input (\code{logical} of length 1) If \code{TRUE}, include the input in the output object.
#' @param verbosity (\code{character} of length 1) Controls the printing of progress updates.
#' The following values are possible: 
#'  \describe{
#'    \item{\code{"none"}}{No progress reports are printed}
#'    \item{\code{"low"}}{Minimal progress reports of a fixed length are printed.}
#'    \item{\code{"high"}}{Lots of information is printed depending on the amount of the input.}
#'  } 
#' @param ... Not used.
#'  
#' @return Returns an object of type \code{classified}
#'
#' @examples
#' \dontrun{
#' # Extract embedded classifications from UNITE FASTA file offline
#' file_path <- system.file("extdata", "unite_general_release.fasta", package = "metacoder")
#' sequences <- ape::read.FASTA(file_path)
#' unite_ex_data_3 <- extract_taxonomy(sequences,
#'                                     regex = "^(.*)\\|(.*)\\|(.*)\\|.*\\|(.*)$",
#'                                     key = c(name = "item_info", seq_id = "item_info",
#'                                             other_id = "item_info", "class_name"),
#'                                     database = "none")
#' # Look up taxonomic data online using sequence ID
#' unite_ex_data <- extract_taxonomy(sequences,
#'                                   regex = "^(.*)\\|(.*)\\|(.*)\\|.*\\|(.*)$",
#'                                 key = c(name = "name", seq_id = "item_id",
#'                                        other_id = "item_info", tax_string = "item_info"))
#' }
#' 
#' @export
#' @rdname extract_taxonomy
extract_taxonomy <- function(input, ...) {
  UseMethod("extract_taxonomy")
}


#' @method extract_taxonomy default
#' @export
#' @rdname extract_taxonomy
extract_taxonomy.default <- function(input,
                                     key = c("class", "taxon_id", "name", "taxon_info", "item_id", "item_info"),
                                     regex = "(.*)",
                                     class_key = c("name", "taxon_id", "taxon_info"),
                                     class_regex = "(.*)",
                                     class_sep = NULL,
                                     class_rev = FALSE,
                                     database = c("none", "ncbi", "itis", "eol", "col", "tropicos", "nbn"),
                                     allow_na = TRUE,
                                     vigilance = c("error", "warning", "message", "none"),
                                     return_match = FALSE,
                                     return_input = TRUE,
                                     verbosity = c("low", "none", "high"),
                                     ...) {
  my_print <- function(text, level = "low") {
    options <- c("none", "low", "high")
    level <- factor(level, ordered = TRUE, levels = options)
    if (level <= verbosity[1]) { message(text) }
  }
  # Validate and standardize input ----------------------------------------------------------------
  my_print(level = "high", "Validating input ----------------------------------")
  # vigilance 
  vigilance <- match.arg(vigilance)
  # verbosity
  verbosity <- match.arg(verbosity)
  # input
  input <- validate_regex_match(input, regex)
  # regex and key
  if (missing(key)) { key <- key[1] }
  key <- validate_regex_key_pair(regex, key, multiple_allowed = c("taxon_info", "item_info"))
  # classification regex and key
  if (missing(class_key)) { class_key <- class_key[1] }
  class_key <- validate_regex_key_pair(class_regex, class_key, multiple_allowed = c("taxon_info"))
  # classification sep
  if (!is.null(class_sep) && (class(class_sep) != "character" | length(class_sep) != 1)) {
    stop('"class_sep" must be a character vector of length 1 or NULL')
  }
  # classification order
  if (class(class_rev) != "logical" | length(class_rev) != 1) {
    stop('"class_rev" must be a logical (aka boolean) vector of length 1')
  }
  # database
  database <- match.arg(database)
  if (database == "none" && ! "class" %in% key) {
    stop("Cannot look up data without a `database` specified.")
  }

  # Parse input -----------------------------------------------------------------------------------
  my_print(level = "high", "Parsing input -------------------------------------")
  parsed_input <- data.frame(stringr::str_match(input, regex), stringsAsFactors = FALSE)
  colnames(parsed_input) <- c("match", key)
  # Consolidate item data
  item_data <- parsed_input
  colnames(item_data) <- c("match", names(key))
  item_data <- item_data[ , c(TRUE, key %in% c("item_id", "item_info")), drop = FALSE]
  if (! return_match) { item_data <- item_data[, -1, drop = FALSE] }
  if (return_input) { item_data <- cbind(data.frame(input = input), item_data) }

  # Determine item classifications ----------------------------------------------------------------
  # This step produces a list of dataframes corresponding the in input values.
  # The rows of each data.frame in the list correspond to taxa in a classification.
  # Columns correspond to information for each taxon.
  my_print(level = "high", "Getting item classifications ----------------------")
  precedence <- c("class", "taxon_id", "item_id", "name")
  classification_data <- parsed_input[ , precedence[precedence %in% key], drop = FALSE] # extract and order data that can be used to get classifications
  classification_func <- get(paste0("class_from_", colnames(classification_data)[1]))
  current_arg_values <- mget(names(formals(extract_taxonomy.default)))
  current_arg_values <- current_arg_values[! names(current_arg_values) %in% c( "input", "...")]
  item_classifications <- do.call(classification_func, c(classification_data, current_arg_values))

  # Infer taxonomy structure ----------------------------------------------------------------------
  my_print(level = "high", "Inferring taxonomic structure ---------------------")
  class_precedence <- c("taxon_id", "name")
  class_source <- class_precedence[class_precedence %in% class_key][1]
  taxonomy <- class_to_taxonomy(item_classifications, id_column = class_source, item_data = item_data) # returns an `classified` object with no item data

  # Add taxon info columns to taxon_data ----------------------------------------------------------
  my_print(level = "high", "Formatting output ---------------------------------")
  taxon_info_column <- function(content, col_name) {
    taxon_values <- lapply(split(content, taxonomy$item_taxon_id), unique)
    if (any(lapply(taxon_values, length) > 1)) {
      stop(paste0('Values for "', col_name, '" are not consistent with the inferred taxonomy (More than one unique value found for at least one taxon). Perhaps a "item_info" key value would be more appropriate?'))
    }
    unlist(taxon_values)[taxonomy$taxon_id]
  }
  if ("taxon_info" %in% key) {
    taxon_info_col_names <- names(key)[key == "taxon_info"]
    taxon_info_source_cols <- setNames(parsed_input[ , colnames(parsed_input) == "taxon_info", drop = FALSE],
                                       taxon_info_col_names)
    new_columns <- mapply(taxon_info_column, taxon_info_source_cols, taxon_info_col_names,
                          SIMPLIFY = FALSE)
    taxonomy$taxon_data <- cbind(taxonomy$taxon_data, new_columns, stringsAsFactors = FALSE) 
  }
  
  # Rename duplicated column names
  colnames(taxonomy$taxon_data) <- rename_duplicated(colnames(taxonomy$taxon_data))
  colnames(taxonomy$item_data) <- rename_duplicated(colnames(taxonomy$item_data))

  # Return output
  my_print(level = "low",
           paste0(length(input), " inputs used to classify ", nrow(taxonomy$item_data),
                  " items by ", length(taxonomy$taxon_id), " taxa."))
  return(taxonomy)
  }




#' @method extract_taxonomy DNAbin
#' @export
#' @rdname extract_taxonomy
extract_taxonomy.DNAbin <- function(input, ...) {
  output <- extract_taxonomy(names(input), ...)
  output$item_data$sequence <- unlist(lapply(as.character(input),
                                             function(x) paste0(x, collapse = "")))
  return(output)
}





