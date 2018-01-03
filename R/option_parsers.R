#' Check that an object is a taxmap
#' 
#' Check that an object is a taxmap
#' This is intended to be used to parse options in other functions.
#' 
#' @param obj A taxmap object
#' 
#' @family option parsers
#' 
#' @keywords internal
verify_taxmap <- function(obj) {
  if (! "Taxmap" %in% class(obj)) {
    stop(paste0('The object supplied is not a taxmap object. ',
                'It appears to be of type "', class(obj)[1], '"'), 
         call. = FALSE)
  }
}


#' Get a data set from a taxmap object
#' 
#' Get a data set from a taxmap object and complain if it does not exist.
#' This is intended to be used to parse options in other functions.
#' 
#' @param obj A taxmap object
#' @param dataset Which data set to use. Can be any of the following:
#'   \describe{
#'     \item{Name}{The name of the data set to use.}
#'     \item{Index}{The index of the data set to use.}
#'     \item{TRUE/FALSE vector}{A TRUE/FALSE vector the same length as the
#'     number of datasets, with exactly one TRUE corresponding to the
#'     selected data set.}
#'   }
#' 
#' @family option parsers
#' 
#' @keywords internal
#' 
#' @examples
#' \dontrun{
#' # Parse dataset
#' x = parse_tax_data(hmp_otus, class_cols = "lineage", class_sep = ";",
#'                    class_key = c(tax_rank = "info", tax_name = "taxon_name"),
#'                    class_regex = "^(.+)__(.+)$")
#'                    
#' # Get data set by name
#' print(metacoder:::get_taxmap_table(x, "tax_data"))
#' print(metacoder:::get_taxmap_table(x, "invalid"))
#' 
#' # Get data set by index
#' print(metacoder:::get_taxmap_table(x, 1))
#' print(metacoder:::get_taxmap_table(x, 3)) # invalid
#' 
#' # Get data set by T/F vector
#' print(metacoder:::get_taxmap_table(x, c(T, F)))
#' print(metacoder:::get_taxmap_table(x, c(T, T))) # invalid
#' print(metacoder:::get_taxmap_table(x, c(T, F, F))) # invalid
#'                    
#' }
get_taxmap_data <- function(obj, dataset) {
  # Check that obj is a taxmap object
  verify_taxmap(obj)
  
  # Convert logicals to numerics 
  if (is.logical(dataset)) {
    if (length(dataset) != length(obj$data)) {
      stop("When using a TRUE/FALSE vector to specify the data set, it must be the same length as the number of data sets",
           call. = FALSE)
    } else {
      dataset <- which(dataset)
    }
  }
  
  # Check for multiple/no values
  if (length(dataset) == 0) {
    stop('No dataset specified.', call. = FALSE)
  }
  if (length(dataset) > 1) {
    stop('Only one dataset can be used.', call. = FALSE)
  }
  
  # Check that dataset exists
  error_msg <- paste0('The dataset "', dataset,
                      '" is not in the object supplied. Datasets found include:\n  ',
                      limited_print(paste0("[", seq_along(obj$data), "] ", names(obj$data)),
                                    type = "silent"))
  if (is.character(dataset)) {
    if (! dataset %in% names(obj$data)) {
      stop(error_msg, call. = FALSE)
    }
  } else if (is.numeric(dataset)) {
    if (! dataset %in% seq_along(obj$data)) {
      stop(error_msg, call. = FALSE)
    }
  }
  
  # Return without printing
  return(invisible(obj$data[[dataset]]))
}


#' Get a table from a taxmap object
#' 
#' Get a table from a taxmap object and complain if it does not exist.
#' This is intended to be used to parse options in other functions.
#' 
#' @inheritParams get_taxmap_data
#' 
#' @return A table
#' 
#' @family option parsers
#' 
#' @keywords internal
get_taxmap_table <- function(obj, dataset) {
  # Get the data set and do checks
  table <- get_taxmap_data(obj, dataset)
  
  # Check that the dataset is a table
  if (! is.data.frame(table)) {
    stop(paste0('The dataset "', dataset,  '" is not a table.'), call. = FALSE)
  }
  
  # Return without printing
  return(invisible(table))
}


#' Get a column subset 
#' 
#' Convert logical, names, or indexes to column names and check that they exist.
#' 
#' @param obj A taxmap object
#' @param dataset The name of a table in \code{obj} that contains counts.
#' @param cols The columns in the data set to use. Takes one of
#'   the following inputs:
#'   \describe{
#'     \item{TRUE/FALSE:}{All non-target columns will be preserved or not.}
#'     \item{Vector of TRUE/FALSE of length equal to the number of columns:}{Preserve the columns
#'   corresponding to \code{TRUE} values.}
#'     \item{Character vector:}{The names of columns to preserve}
#'     \item{Numeric vector:}{The indexes of columns to preserve}
#'   }
#' 
#' @keywords internal
#' 
#' @family option parsers
#' 
#' @examples
#' \dontrun{
#' # Parse dataset
#' x = parse_tax_data(hmp_otus, class_cols = "lineage", class_sep = ";",
#'                    class_key = c(tax_rank = "info", tax_name = "taxon_name"),
#'                    class_regex = "^(.+)__(.+)$")
#'                    
#' # Get all col names
#' metacoder:::parse_taxmap_cols(x, "tax_data")
#' 
#' # Get col names by index
#' metacoder:::parse_taxmap_cols(x, "tax_data", 2:4)
#' 
#' # Get col names by name (i.e. verify)
#' metacoder:::parse_taxmap_cols(x, "tax_data", c("taxon_id", "lineage"))
#' metacoder:::parse_taxmap_cols(x, "tax_data", c("taxon_id", "not_valid"))
#' 
#' # Get colnames by TRUE/FALSE vector
#' metacoder:::parse_taxmap_cols(x, "tax_data", startsWith(colnames(x$data$tax_data), "7"))
#'                    
#' }
get_taxmap_cols <- function(obj, dataset, cols = NULL) {
  # Get table used. This checks the obj as well
  my_table <- get_taxmap_table(obj, dataset)
  
  # If NULL, return all cols
  if (is.null(cols)) {
    cols <- TRUE
  }
  
  # Convert logical/numeric to column names
  if (is.logical(cols)) {
    if (length(cols) == ncol(my_table)) { # Is a TRUE/FALSE vector
      result <- colnames(my_table)[cols]
    } else if (length(cols) == 1) { # Is a single TRUE/FALSE
      if (cols) {
        result <- colnames(my_table)
      } else {
        result <- character(0)
      }
    } else { # Incorrect length of TRUE/FALSE vector
      stop(paste0("When specifying columns with a TRUE/FALSE vector, it must be", 
                  "either a single value or have a length equal to the number of", 
                  "columns in '", dataset, "'."),
           call. = FALSE)
    }
  } else if (is.character(cols)) { # Is already column names
    invalid_cols <- cols[! cols %in% colnames(my_table)]
    if (length(invalid_cols) == 0) {
      result <- cols
    } else {
      stop(paste0('The following ', length(invalid_cols), ' column(s) are not in "', dataset, '":\n  ',
                  limited_print(invalid_cols, type = "silent")),
           call. = FALSE)
    }
  } else if (is.numeric(cols)) { # If column indexes
    invalid_cols <- cols[cols > ncol(my_table) | cols < 1]
    if (length(invalid_cols) == 0) {
      result <- colnames(my_table)[cols]
    } else {
      stop(paste0('The following ', length(invalid_cols), ' column indexes are not valid for "', dataset, '":\n  ',
                  limited_print(invalid_cols, type = "silent")),
           call. = FALSE)
    }
  } else {
    stop(paste0("`cols` is of the invalid type: ", class(cols), ".\n", 
                'The "cols" option must either be TRUE/FALSE or a vector of valid column names/indexes.',
                call. = FALSE))
  }
  
  # Retrun result
  return(result)
}


#' Parse the other_cols option
#'
#' Parse the other_cols option used in many calculation functions.
#'
#' @param obj A taxmap object
#' @param dataset The name of a table in \code{obj} that contains counts.
#' @param cols The names/indexes of columns in \code{dataset} to use. Takes one
#'   of the following inputs:
#'   \describe{
#'     \item{TRUE/FALSE:}{All columns will used.}
#'     \item{Vector of TRUE/FALSE of length equal to the number of columns:}{Use the columns
#'   corresponding to \code{TRUE} values.}
#'     \item{Character vector:}{The names of columns to use}
#'     \item{Numeric vector:}{The indexes of columns to use}
#'   }
#' @param other_cols Preserve in the output non-target columns present in the
#'   input data. The "taxon_id" column will always be preserved. Takes one of
#'   the following inputs:
#'   \describe{
#'     \item{TRUE/FALSE:}{All non-target columns will be preserved or not.}
#'     \item{Vector of TRUE/FALSE of length equal to the number of columns:}{Preserve the columns
#'   corresponding to \code{TRUE} values.}
#'     \item{Character vector:}{The names of columns to preserve}
#'     \item{Numeric vector:}{The indexes of columns to preserve}
#'   }
#'   
#' @keywords internal
#' 
#' @family option parsers
#' 
#' @examples
#' \dontrun{
#' # Parse dataset for examples
#' x = parse_tax_data(hmp_otus, class_cols = "lineage", class_sep = ";",
#'                    class_key = c(tax_rank = "info", tax_name = "taxon_name"),
#'                    class_regex = "^(.+)__(.+)$")
#' 
#' # If all cols are used, there are no other cols, only "taxon_id"
#' metacoder:::get_taxmap_other_cols(x, dataset = "tax_data", cols = TRUE)
#' 
#' # If a subset of target columns is specified, the rest are returned 
#' metacoder:::get_taxmap_other_cols(x, dataset = "tax_data", cols = 2:3)
#' 
#' # Additionally, a subset of other columns can be specified
#' metacoder:::get_taxmap_other_cols(x, dataset = "tax_data", cols = 2:3,
#'                                   other_cols = 4:5)
#'                    
#' }
get_taxmap_other_cols <- function(obj, dataset, cols, other_cols = NULL) {
  # Get table used 
  my_table <- get_taxmap_table(obj, dataset)
  
  # Get target cols
  cols <- get_taxmap_cols(obj, dataset, cols)
  
  # Get other cols
  other_cols <- get_taxmap_cols(obj, dataset, other_cols)
  
  # Remove target cols if present
  in_both <- other_cols %in% cols
  if (sum(in_both) > 0) {
    warning(paste0("The following columns will be replaced in the output:\n  ",
                   limited_print(other_cols[in_both], type = "silent")),
            call. = FALSE)
  }
  result <- other_cols[! in_both]
  
  # Add taxon id column regardless
  if (! "taxon_id" %in% result) {
    result <- c("taxon_id", result)
  }
  
  return(result)
}


#' Read sequences in an unknown format
#'
#' Read sequences in an unknown format. This is meant to parse the sequence
#' input arguments of functions like \code{\link{primersearch}}.
#' 
#' @param input (\code{character}) One of the following: 
#' \describe{
#'   \item{A character vector of sequences}{See the example below for what this
#'   looks like. The parser \code{\link{read_fasta}} produces output like this.}
#'   \item{A list of character vectors}{Each vector should have one base per element.}
#'   \item{A "DNAbin" object}{This is the result of parsers like
#'   \code{\link[ape]{read.FASTA}}.}
#'   \item{A list of "SeqFastadna" objects}{This is the result of parsers like
#'   \code{\link[seqinr]{read.fasta}}.}
#'   Either "input" or "file" must be supplied but not both.
#' }
#' @param file The path to a FASTA file containing sequences to use. Either
#'   "input" or "file" must be supplied but not both.
#' 
#' @return A named character vector of sequences
#' 
#' @keywords internal
parse_seq_input <- function(input = NULL, file = NULL) {
  # Check parameters
  if (sum(! c(is.null(file), is.null(input))) != 1) {
    stop(call. = FALSE,
         "Either `file` or `input` must be supplied, but not both.")
  }
  
  if (! is.null(file) && (! is.character(file) || length(file) != 1)) {
    stop(call. = FALSE,
         "`file` must be a character vector of length 1 that is a valid path to a file.")
  }
  
  # Convert to common format
  if (! is.null(file)) {
    result <- read_fasta(file)
  } else if (length(input) == 0 || class(input) == "character") {
    result <- input
  } else if (class(input) == "DNAbin") {
    result <- toupper(vapply(as.character(input), paste, character(1), collapse = ""))
  } else if (class(input[[1]]) == "SeqFastadna" || class(input) == "list") {
    result <- vapply(input, paste, character(1), collapse = "")
  } else {
    stop(paste0('Could not parse sequence information of class "', class(input), '".'),
         call. = FALSE)
  }
  
  return(result)
}