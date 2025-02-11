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
#' NOTE: This will be replaced by the function `get_dataset` in the `taxa`
#' package. Get a data set from a taxmap object and complain if it does not
#' exist. This is intended to be used to parse options in other functions.
#' 
#' @param obj A taxmap object
#' @param data Which data set to use. Can be any of the following:
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
get_taxmap_data <- function(obj, data) {
  # Check that obj is a taxmap object
  verify_taxmap(obj)
  
  # Convert logicals to numerics 
  if (is.logical(data)) {
    if (length(data) != length(obj$data)) {
      stop("When using a TRUE/FALSE vector to specify the data set, it must be the same length as the number of data sets",
           call. = FALSE)
    } else {
      data <- which(data)
    }
  }
  
  # Check for multiple/no values
  if (length(data) == 0) {
    stop('No dataset specified.', call. = FALSE)
  }
  if (length(data) > 1) {
    stop('Only one dataset can be used.', call. = FALSE)
  }
  
  # Check that data exists
  error_msg <- paste0('The dataset "', data,
                      '" is not in the object supplied. Datasets found include:\n  ',
                      limited_print(paste0("[", seq_along(obj$data), "] ", names(obj$data)),
                                    type = "silent"))
  if (is.character(data)) {
    if (! data %in% names(obj$data)) {
      stop(error_msg, call. = FALSE)
    }
  } else if (is.numeric(data)) {
    if (! data %in% seq_along(obj$data)) {
      stop(error_msg, call. = FALSE)
    }
  }
  
  # Return without printing
  return(invisible(obj$data[[data]]))
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
get_taxmap_table <- function(obj, data) {
  # Get the data set and do checks
  table <- get_taxmap_data(obj, data)
  
  # Check that the data is a table
  if (! is.data.frame(table)) {
    stop(paste0('The dataset "', data,  '" is not a table.'), call. = FALSE)
  }
  
  # Return without printing
  return(invisible(table))
}


#' Get a column subset 
#' 
#' Convert logical, names, or indexes to column names and check that they exist.
#' 
#' @param obj A taxmap object
#' @param data The name of a table in \code{obj} that contains counts.
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
get_taxmap_cols <- function(obj, data, cols = NULL) {
  # Get table used. This checks the obj as well
  my_table <- get_taxmap_table(obj, data)
  
  # If NULL, return all cols
  if (is.null(cols)) {
    cols <- TRUE
  }
  
  # Convert factors to characters
  if (is.factor(cols)) {
    cols <- as.character(cols)
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
                  "columns in '", data, "'."),
           call. = FALSE)
    }
  } else if (is.character(cols)) { # Is already column names
    invalid_cols <- cols[! cols %in% colnames(my_table)]
    if (length(invalid_cols) == 0) {
      result <- cols
    } else {
      stop(paste0('The following ', length(invalid_cols), ' column(s) are not in "', data, '":\n  ',
                  limited_print(invalid_cols, type = "silent")),
           call. = FALSE)
    }
  } else if (is.numeric(cols)) { # If column indexes
    invalid_cols <- cols[cols > ncol(my_table) | cols < 1]
    if (length(invalid_cols) == 0) {
      result <- colnames(my_table)[cols]
    } else {
      stop(paste0('The following ', length(invalid_cols), ' column indexes are not valid for "', data, '":\n  ',
                  limited_print(invalid_cols, type = "silent")),
           call. = FALSE)
    }
  } else {
    stop(paste0("`cols` is of the invalid type: ", class(cols), ".\n", 
                'The "cols" option must either be TRUE/FALSE or a vector of valid column names/indexes.'),
         call. = FALSE)
  }
  
  # Retrun result
  return(result)
}


#' Parse the other_cols option
#'
#' Parse the other_cols option used in many calculation functions.
#'
#' @param obj A taxmap object
#' @param data The name of a table in \code{obj} that contains counts.
#' @param cols The names/indexes of columns in \code{data} to use. Takes one
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
get_taxmap_other_cols <- function(obj, data, cols, other_cols = NULL) {
  # Get table used 
  my_table <- get_taxmap_table(obj, data)
  
  # Get target cols
  cols <- get_taxmap_cols(obj, data, cols)
  
  # Get other cols
  other_cols <- get_taxmap_cols(obj, data, other_cols)
  
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
#' @param output_format The format of the sequences returned. Either "character" or "DNAbin".
#' @param u_to_t If `TRUE`, then "U" in the sequence will be converted to "T".
#' 
#' @return A named character vector of sequences
#' 
#' @keywords internal
parse_seq_input <- function(input = NULL, file = NULL, output_format = "character", u_to_t = FALSE) {
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
  if (output_format == "character") {
    if (! is.null(file)) {
      result <- read_fasta(file)
    } else if (length(input) == 0 || inherits(input, "character")) {
      result <- input
    } else if (inherits(input,"DNAbin")) {
      result <- toupper(vapply(as.character(input), paste, character(1), collapse = ""))
    } else if (inherits(input[[1]], "SeqFastadna") || inherits(input, "list")) {
      result <- vapply(input, paste, character(1), collapse = "")
    } else {
      stop(paste0('Could not parse sequence information of class "', class(input), '".'),
           call. = FALSE)
    }
    
    if (u_to_t) {
      result <- vapply(result, FUN = gsub, FUN.VALUE = character(1),
                       pattern = "U", replacement = "T", fixed = TRUE)
      result <- vapply(result, FUN = gsub, FUN.VALUE = character(1),
                       pattern = "u", replacement = "t", fixed = TRUE)
    }
    
  } else if (output_format == "DNAbin") {
    if (! is.null(file)) {
      if (u_to_t) {
        file <- make_fasta_with_u_replaced(file)
      }
      result <- ape::read.FASTA(file)
    } else if (length(input) == 0 || inherits(input, "character")) {
      if (u_to_t) {
        input <- vapply(input, FUN = gsub, FUN.VALUE = character(1),
                        pattern = "U", replacement = "T", fixed = TRUE)
        input <- vapply(input, FUN = gsub, FUN.VALUE = character(1),
                        pattern = "u", replacement = "t", fixed = TRUE)
      }
      result <- ape::as.DNAbin(strsplit(input, split = ""))
    } else if (inherits(input,  "DNAbin")) {
      result <- input
    } else if (inherits(input[[1]], "SeqFastadna") || inherits(input, "list")) {
      input <- lapply(input, function(x) {
        attributes(x) <- NULL
        return(x)
      })
      if (u_to_t) {
        input <- lapply(input, gsub, pattern = "U", replacement = "T", fixed = TRUE)
        input <- lapply(input, gsub, pattern = "u", replacement = "t", fixed = TRUE)
      }
      result <- ape::as.DNAbin(input)
    } else {
      stop(paste0('Could not parse sequence information of class "', class(input), '".'),
           call. = FALSE)
    }
  } else {
    stop(paste0('Invalid output format "', output_format, '".'))
  }
  
  return(result)
}
